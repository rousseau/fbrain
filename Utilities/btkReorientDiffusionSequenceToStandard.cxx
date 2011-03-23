/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/11/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr)

  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.

==========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "CmdLine.h"
#include "vnl/vnl_cross.h"

#include "itkEuler3DTransform.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkJoinSeriesImageFilter.h"

#include "btkDiffusionGradientTable.h"

#include "vnl/vnl_inverse.h"

// Local includes
#include "btkLandmarksFileReader.h"


int main( int argc, char *argv[] )
{

  try {

  // Parse arguments
  const char *inputImageFile = NULL, *outputImageFile = NULL, *lmksFile = NULL;
  const char *gTableFile = NULL, *cTableFile = NULL;

  TCLAP::CmdLine cmd("Reorient a DWI sequence to standard orientation", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> gTableArg("g","gtable","Gradient table",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Reoriented image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> cTableArg("c","ctable","Corrected table",true,"","string",cmd);
  TCLAP::ValueArg<std::string> lmksArg("l","landmarks","Landmarks file",true,"","string",cmd);
  TCLAP::SwitchArg nnSwitch("","nn","Nearest Neighbor interpolation", cmd, false);
  TCLAP::SwitchArg fovSwitch("e","extendFOV","Extends the FOV to avoid image cropping", cmd, false);

  cmd.parse( argc, argv );

  inputImageFile  = inputArg.getValue().c_str();
  gTableFile      = gTableArg.getValue().c_str();
  outputImageFile = outputArg.getValue().c_str();
  cTableFile      = cTableArg.getValue().c_str();
  lmksFile        = lmksArg.getValue().c_str();
  bool nn         = nnSwitch.getValue();

  const    unsigned int    Dimension = 4;

  // Read image

  typedef itk::Image< short, Dimension >  ImageType;
  typedef itk::Image< short, 3 >  Image3DType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
  imageReader -> SetFileName(  inputImageFile );
  imageReader -> Update();
  ImageType::Pointer image = imageReader->GetOutput();

  btk::btkLandmarksFileReader::Pointer landRead = btk::btkLandmarksFileReader::New();
  landRead->SetInputFileName(lmksFile);
  landRead->Update();

  double *lpt_ras = landRead->GetOutputLPT();
  double *rpt_ras = landRead->GetOutputRPT();
  double *apt_ras = landRead->GetOutputAPT();
  double *ppt_ras = landRead->GetOutputPPT();

  // Compute rotation matrix

  vnl_vector<double> left_lps(3), pos_lps(3), sup_lps(3);

  left_lps(0) = rpt_ras[0] - lpt_ras[0];
  left_lps(1) = rpt_ras[1] - lpt_ras[1];
  left_lps(2) = lpt_ras[2] - rpt_ras[2];

  left_lps.normalize();

  pos_lps(0) = apt_ras[0] - ppt_ras[0];
  pos_lps(1) = apt_ras[1] - ppt_ras[1];
  pos_lps(2) = ppt_ras[2] - apt_ras[2];

  pos_lps.normalize();

  sup_lps = vnl_cross_3d(left_lps,pos_lps);

  // TODO The following code is quite repeated throug the library
  // we should include it somewhere else

  vnl_matrix<double> R(3,3);
  R.set_identity();

  R(0,0) = left_lps(0); R(0,1) = left_lps(1); R(0,2) = left_lps(2);
  R(1,0) = pos_lps(0) ; R(1,1) = pos_lps(1);  R(1,2) = pos_lps(2);
  R(2,0) = sup_lps(0) ; R(2,1) = sup_lps(1);  R(2,2) = sup_lps(2);

  // Corrects rotation matrix

  vnl_matrix<double> PQ = R;
  vnl_matrix<double> NQ = R;
  vnl_matrix<double> PQNQDiff;

  const unsigned int maximumIterations = 100;

  for(unsigned int ni = 0; ni < maximumIterations; ni++ )
  {
    // Average current Qi with its inverse transpose
    NQ = ( PQ + vnl_inverse_transpose( PQ ) ) / 2.0;
    PQNQDiff = NQ - PQ;
    if( PQNQDiff.frobenius_norm() < 1e-7 )
    {
//      std::cout << "Polar decomposition used " << ni << " iterations " << std::endl;
      break;
    }
    else
    {
      PQ = NQ;
    }
  }

  // Correct and write image

  typedef itk::Euler3DTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();

  ImageType::SizeType size = image -> GetLargestPossibleRegion().GetSize();
  ImageType::IndexType index = image -> GetLargestPossibleRegion().GetIndex();
  unsigned int numberOfFrames = size[3];

  typedef itk::ContinuousIndex< double, 3 > ContinuousIndexType;
  ContinuousIndexType centerIndex;

  centerIndex[0] = (size[0] - 1.0)/2.0;
  centerIndex[1] = (size[1] - 1.0)/2.0;
  centerIndex[2] = (size[2] - 1.0)/2.0;

  Image3DType::PointType centerPoint;

  typedef itk::ExtractImageFilter< ImageType, Image3DType > ExtractorType;
  ExtractorType::Pointer extractor = ExtractorType::New();

  extractor -> SetInput( image );
  index[3] = 0;
  size[3]  = 0;

  ImageType::RegionType region;
  region.SetIndex(index);
  region.SetSize(size);

  extractor -> SetExtractionRegion( region );
  extractor -> Update();

  Image3DType::Pointer b0 = extractor -> GetOutput();

  b0 -> TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);

  transform -> SetCenter( centerPoint );
  transform -> SetRotationMatrix( NQ );
//  TODO Cannot constrain the rotation because sometimes it is necessary to rotate in xy
//  Something I could do is to round the rotations in xy to n*pi*/2 ...
//  transform -> SetRotation( 0.0, 0.0, transform -> GetAngleZ() );

  // Extend FOV according to the transformation
  Image3DType::SizeType  resamplingSize;
  Image3DType::IndexType lowerIndex; lowerIndex.Fill(0);
  Image3DType::IndexType upperIndex; upperIndex.Fill(0);

  if ( fovSwitch.isSet() )
  {
    typedef itk::ImageRegionIteratorWithIndex< Image3DType > Image3DIteratorType;
    Image3DIteratorType b0It( b0, b0 -> GetLargestPossibleRegion() );

    Image3DType::IndexType b0Index;
    Image3DType::IndexType b0TransformedIndex;

    Image3DType::PointType b0Point;
    Image3DType::PointType b0TransformedPoint;

    // Even when only the corners are necessary, the code is shorter in this way
    // (and it's not very expensive to compute it)

    for(b0It.GoToBegin(); !b0It.IsAtEnd(); ++b0It)
    {
      b0Index = b0It.GetIndex();
      b0 -> TransformIndexToPhysicalPoint(b0Index, b0Point);
      b0TransformedPoint = transform -> TransformPoint(b0Point);
      b0 -> TransformPhysicalPointToIndex(b0TransformedPoint, b0TransformedIndex);

      for (unsigned int i = 0; i<3; i++)
        if ( b0TransformedIndex[i] < lowerIndex[i] ) lowerIndex[i] = b0TransformedIndex[i];
        else if
           ( b0TransformedIndex[i] > upperIndex[i] ) upperIndex[i] = b0TransformedIndex[i];

    }


    resamplingSize[0] = upperIndex[0] - lowerIndex[0] + 1;
    resamplingSize[1] = upperIndex[1] - lowerIndex[1] + 1;
    resamplingSize[2] = upperIndex[2] - lowerIndex[2] + 1;

  } else
  {
    resamplingSize = b0 -> GetLargestPossibleRegion().GetSize();
  }

  Image3DType::PointType resamplingOrigin;
  b0 -> TransformIndexToPhysicalPoint(lowerIndex, resamplingOrigin );

  transform -> SetRotationMatrix( NQ.transpose() );

  // Resampling

  typedef itk::ResampleImageFilter< Image3DType, Image3DType >  ResamplerType;

  typedef itk::NearestNeighborInterpolateImageFunction< Image3DType >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::JoinSeriesImageFilter< Image3DType, ImageType > JoinerType;
  JoinerType::Pointer joiner = JoinerType::New();
  joiner -> SetOrigin( 0.0 );
  // TODO Does it change something the 4th dimension spacing in diffusion sequences?
  joiner -> SetSpacing( 1 );

  for (unsigned int i=0; i<numberOfFrames; i++)
  {
    index[3] = i;
    size[3]  = 0;

    ImageType::RegionType region;
    region.SetIndex(index);
    region.SetSize(size);

    extractor -> SetExtractionRegion( region );
    extractor -> Update();

    ResamplerType::Pointer resampler = ResamplerType::New();

    if (nn)
    {
      resampler -> SetInterpolator(interpolator);
    }

    resampler -> SetInput( extractor -> GetOutput() );
    resampler -> SetTransform(transform);
    resampler -> SetDefaultPixelValue(0);
    resampler -> SetSize( resamplingSize );
    resampler -> SetOutputOrigin( resamplingOrigin );
    resampler -> SetOutputDirection( extractor -> GetOutput() -> GetDirection()  );
    resampler -> SetOutputSpacing( extractor -> GetOutput() -> GetSpacing() );

    resampler -> Update();

    joiner -> SetInput( i, resampler -> GetOutput() );

  }

  joiner -> Update();

  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer  imageWriter  = ImageWriterType::New();
  imageWriter -> SetInput( joiner -> GetOutput() );
  imageWriter -> SetFileName(  outputImageFile );
  imageWriter -> Update();

  // Change gradient table

  typedef btk::DiffusionGradientTable< ImageType > GradientTableType;
  GradientTableType::Pointer gradientTable = GradientTableType::New();

  gradientTable -> SetNumberOfGradients(numberOfFrames);
  gradientTable -> SetImage( image );
  gradientTable -> SetTransform( transform );
  gradientTable -> LoadFromFile( gTableFile);
  gradientTable -> RotateGradientsInWorldCoordinates();
  gradientTable -> SaveToFile( cTableFile);


  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

