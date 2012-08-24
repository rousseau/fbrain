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

/* Standard includes */
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkEuler3DTransform.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkOrientImageFilter.h"

/* Vnl includes */
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_cross.h"

/* Btk includes */
#include "btkLandmarksFileReader.h"
#include "btkDiffusionGradientTable.h"
#include "btkFileHelper.h"


int main( int argc, char *argv[] )
{

  try {

  // Parse arguments
  const char *lmksFile = NULL;
  const char *mgradResampled = NULL;
  const char *mgradName = NULL;

  TCLAP::CmdLine cmd("Reorient a DWI sequence to standard orientation", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Reoriented image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> lmksArg("l","landmarks","Landmarks file",true,"","string",cmd);
  TCLAP::SwitchArg nnSwitch("","nn","Nearest Neighbor interpolation", cmd, false);
  TCLAP::SwitchArg fovSwitch("e","extendFOV","Extends the FOV to avoid image cropping", cmd, false);
  TCLAP::ValueArg<std::string> mgradArg("","mean-gradient","Mean gradient",false,"","string",cmd);
  TCLAP::ValueArg<std::string> mgradResampledArg("","mean-gradient-resampled","Mean gradient resampled",false,"","string",cmd);

  cmd.parse( argc, argv );

  lmksFile        = lmksArg.getValue().c_str();
  bool nn         = nnSwitch.getValue();

  mgradName = mgradArg.getValue().c_str();
  mgradResampled = mgradResampledArg.getValue().c_str();


  const    unsigned int    Dimension = 4;


  std::string inputImageFile = inputArg.getValue();
  std::string inRadix = btk::FileHelper::GetRadixOf(inputImageFile);
  std::string bvec_in = inRadix + ".bvec";
  std::string bval_in = inRadix + ".bval";

  std::string outputImageFile = outputArg.getValue();
  std::string outRadix = btk::FileHelper::GetRadixOf(outputImageFile);
  std::string bvec_out = outRadix + ".bvec";
  std::string bval_out = outRadix + ".bval";

  // Read image

  typedef itk::Image< short, Dimension >  ImageType;
  typedef itk::Image< short, 3 >  Image3DType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
  imageReader -> SetFileName(  inputImageFile.c_str() );
  imageReader -> Update();
  ImageType::Pointer image = imageReader->GetOutput();

  typedef itk::ImageFileReader< Image3DType > Image3DReaderType;

  Image3DType::Pointer mgrad;
  if ( strcmp(mgradName,"") )
  {
    Image3DReaderType::Pointer    mgradReader    = Image3DReaderType::New();
    mgradReader -> SetFileName( mgradName );
    mgradReader -> Update();
    mgrad = mgradReader -> GetOutput();
  }

  // Read landmarks

  btk::btkLandmarksFileReader::Pointer landRead = btk::btkLandmarksFileReader::New();
  landRead->SetInputFileName(lmksFile);
  landRead->Update();

  double *lpt_ras = landRead->GetOutputLPT();
  double *rpt_ras = landRead->GetOutputRPT();
  double *apt_ras = landRead->GetOutputAPT();
  double *ppt_ras = landRead->GetOutputPPT();

  // Convert points to lps coordinates

  Image3DType::PointType leftPoint, rightPoint, anteriorPoint, posteriorPoint;

  leftPoint[0]      = -lpt_ras[0]; leftPoint[1]      = -lpt_ras[1]; leftPoint[2]      = lpt_ras[2];
  rightPoint[0]     = -rpt_ras[0]; rightPoint[1]     = -rpt_ras[1]; rightPoint[2]     = rpt_ras[2];
  anteriorPoint[0]  = -apt_ras[0]; anteriorPoint[1]  = -apt_ras[1]; anteriorPoint[2]  = apt_ras[2];
  posteriorPoint[0] = -ppt_ras[0]; posteriorPoint[1] = -ppt_ras[1]; posteriorPoint[2] = ppt_ras[2];

  // Extract b0 original sequence

  ImageType::SizeType     size    = image -> GetLargestPossibleRegion().GetSize();
  ImageType::IndexType    index   = image -> GetLargestPossibleRegion().GetIndex();
  ImageType::SpacingType  spacing = image -> GetSpacing();
  unsigned int numberOfFrames = size[3];

  typedef itk::ExtractImageFilter< ImageType, Image3DType > ExtractorType;
  ExtractorType::Pointer extractor = ExtractorType::New();
  extractor -> SetInput( image );

  ImageType::RegionType region;
  index[3] = 0;
  size[3]  = 0;
  region.SetIndex(index);
  region.SetSize(size);

  extractor -> SetExtractionRegion( region );
  extractor -> SetDirectionCollapseToSubmatrix();
  extractor -> Update();

  Image3DType::Pointer b0_ori = extractor -> GetOutput();

  // Recover original direction matrix

  vnl_matrix<double> oriDirection(3,3);
  oriDirection = b0_ori -> GetDirection().GetVnlMatrix();

  // Reorient mean gradient to RAI

  typedef itk::OrientImageFilter< Image3DType, Image3DType >  FilterType;

  Image3DType::Pointer mgrad_reo;

  if ( strcmp(mgradName,"") )
  {
    FilterType::Pointer filter = FilterType::New();
    filter -> SetInput( mgrad  );
    filter -> UseImageDirectionOn();
    filter -> SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
    filter -> Update();
    mgrad_reo = filter -> GetOutput();
  }

  // Reorient sequence to RAI

  typedef itk::JoinSeriesImageFilter< Image3DType, ImageType > JoinerType;
  JoinerType::Pointer joiner = JoinerType::New();
  joiner -> SetOrigin( 0.0 );
  joiner -> SetSpacing( spacing[3] );

  for (unsigned int i = 0; i < numberOfFrames; i++)
  {
    index[3] = i;

    ImageType::RegionType desiredRegion;
    desiredRegion.SetSize(  size  );
    desiredRegion.SetIndex( index );

    extractor -> SetExtractionRegion( desiredRegion );

    FilterType::Pointer filter = FilterType::New();
    filter -> SetInput( extractor -> GetOutput()  );
    filter -> UseImageDirectionOn();
    filter -> SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
    filter -> Update();

    joiner -> SetInput( i, filter -> GetOutput() );
  }

  joiner -> Update();

  ImageType::Pointer image_reo = joiner -> GetOutput();

  // Store original Direction and origin.
  ImageType::SizeType size_reo = image_reo -> GetLargestPossibleRegion().GetSize();
  ImageType::IndexType index_reo = image_reo -> GetLargestPossibleRegion().GetIndex();

  size_reo[3] = 0;
  index_reo[3] = 0;

  ImageType::RegionType region_reo;
  region_reo.SetIndex(index_reo);
  region_reo.SetSize(size_reo);

  ExtractorType::Pointer extractor_reo = ExtractorType::New();

  extractor_reo -> SetInput( image_reo );
  extractor_reo -> SetExtractionRegion( region_reo );
  extractor_reo -> SetDirectionCollapseToSubmatrix();
  extractor_reo -> Update();

  Image3DType::Pointer b0_reo = extractor_reo -> GetOutput();

  vnl_matrix<double> reoDirection(3,3);
  reoDirection = b0_reo -> GetDirection().GetVnlMatrix();
  vnl_matrix<double> reoDirectionTrans(3,3);
  reoDirectionTrans = reoDirection.transpose();

  // Reorient gradients accordingly

  typedef btk::DiffusionGradientTable< ImageType > GradientTableType;
  GradientTableType::Pointer gradientTable = GradientTableType::New();

  vnl_matrix<double> rotationMatrix;
  rotationMatrix = reoDirectionTrans*oriDirection;

  gradientTable -> SetNumberOfGradients(numberOfFrames);
  gradientTable -> SetRotationMatrix( rotationMatrix );
  gradientTable -> LoadFromFile( bvec_in.c_str());
  gradientTable -> RotateGradients();

  // Create the image in standard orientation
  // p = I*S*idx + O_center

  ImageType::Pointer stdImage = joiner -> GetOutput();

  ImageType::DirectionType stdImageDirection  = stdImage -> GetDirection();
  stdImageDirection.SetIdentity();
  stdImage -> SetDirection (stdImageDirection);

  ImageType::SpacingType stdSpacing  = stdImage -> GetSpacing();
  ImageType::SizeType 	 stdSize     = stdImage -> GetLargestPossibleRegion().GetSize();
  ImageType::IndexType   stdIndex    = stdImage -> GetLargestPossibleRegion().GetIndex();

  ImageType::PointType stdOrigin;
  Image3DType::PointType stdOrigin3D;

  for (unsigned int i=0; i<3; i++)
  {
    stdOrigin[i] = (1.0-stdSize[i])*0.5*stdSpacing[i];
    stdOrigin3D[i] = stdOrigin[i];
  }

  stdImage -> SetOrigin( stdOrigin );

  // Create the mean gradient in standard orientation

  Image3DType::Pointer mgrad_std = mgrad_reo;

  if ( strcmp(mgradName,"") )
  {
    Image3DType::DirectionType mgrad_std_direction  = mgrad_std -> GetDirection();
    mgrad_std_direction.SetIdentity();
    mgrad_std -> SetDirection (mgrad_std_direction);
    mgrad_std -> SetOrigin(stdOrigin3D);
  }

  // Obtain directions LPS directions in standard image

  vnl_vector<double> left_lps(3), pos_lps(3), sup_lps(3);

  left_lps = reoDirectionTrans*( leftPoint.GetVnlVector() - rightPoint.GetVnlVector() );
  left_lps.normalize();

  pos_lps = reoDirectionTrans*( posteriorPoint.GetVnlVector() - anteriorPoint.GetVnlVector() );
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

  ExtractorType::Pointer stdExtractor = ExtractorType::New();
  stdExtractor -> SetInput( stdImage );

  stdIndex[3] = 0;
  stdSize[3]  = 0;

  ImageType::RegionType stdRegion;
  stdRegion.SetIndex(stdIndex);
  stdRegion.SetSize(stdSize);

  stdExtractor -> SetExtractionRegion( stdRegion );
  stdExtractor -> SetDirectionCollapseToSubmatrix();
  stdExtractor -> Update();

  Image3DType::Pointer stdB0 = stdExtractor -> GetOutput();

  typedef itk::Euler3DTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::ContinuousIndex< double, 3 > ContinuousIndexType;
  ContinuousIndexType centerIndex;

  centerIndex[0] = (stdSize[0] - 1.0) / 2.0;
  centerIndex[1] = (stdSize[1] - 1.0) / 2.0;
  centerIndex[2] = (stdSize[2] - 1.0) / 2.0;

  Image3DType::PointType centerPoint;

  stdB0 -> TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);

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
    Image3DIteratorType b0It( stdB0, stdB0 -> GetLargestPossibleRegion() );

    Image3DType::IndexType b0Index;
    Image3DType::IndexType b0TransformedIndex;

    Image3DType::PointType b0Point;
    Image3DType::PointType b0TransformedPoint;

    // Even when only the corners are necessary, the code is shorter in this way
    // (and it's not very expensive to compute it)

    for(b0It.GoToBegin(); !b0It.IsAtEnd(); ++b0It)
    {
      b0Index = b0It.GetIndex();
      stdB0 -> TransformIndexToPhysicalPoint(b0Index, b0Point);
      b0TransformedPoint = transform -> TransformPoint(b0Point);
      stdB0 -> TransformPhysicalPointToIndex(b0TransformedPoint, b0TransformedIndex);

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
    resamplingSize = stdB0 -> GetLargestPossibleRegion().GetSize();
  }

  Image3DType::PointType resamplingOrigin;
  stdB0 -> TransformIndexToPhysicalPoint(lowerIndex, resamplingOrigin );

  transform -> SetRotationMatrix( NQ.transpose() );

  // Resampling

  typedef itk::ResampleImageFilter< Image3DType, Image3DType >  ResamplerType;

  typedef itk::NearestNeighborInterpolateImageFunction< Image3DType >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  typedef itk::JoinSeriesImageFilter< Image3DType, ImageType > JoinerType;
  JoinerType::Pointer joiner2 = JoinerType::New();
  joiner2 -> SetOrigin( 0.0 );
  // TODO Does it change something the 4th dimension spacing in diffusion sequences?
  joiner2 -> SetSpacing( 1 );

  for (unsigned int i=0; i<numberOfFrames; i++)
  {
    stdIndex[3] = i;
    stdSize[3]  = 0;

    ImageType::RegionType stdRegion;
    stdRegion.SetIndex( stdIndex );
    stdRegion.SetSize( stdSize );

    stdExtractor -> SetExtractionRegion( stdRegion );
    stdExtractor -> Update();

    ResamplerType::Pointer resampler = ResamplerType::New();

    if (nn)
    {
      resampler -> SetInterpolator(interpolator);
    }

    resampler -> SetInput( stdExtractor -> GetOutput() );
    resampler -> SetTransform(transform);
    resampler -> SetDefaultPixelValue(0);
    resampler -> SetSize( resamplingSize );
    resampler -> SetOutputOrigin( resamplingOrigin );
    resampler -> SetOutputDirection( stdExtractor -> GetOutput() -> GetDirection()  );
    resampler -> SetOutputSpacing( stdExtractor -> GetOutput() -> GetSpacing() );

    resampler -> Update();

    joiner2 -> SetInput( i, resampler -> GetOutput() );

  }

  joiner2 -> Update();


  Image3DType::Pointer mgrad_resampled;

  if ( strcmp(mgradName,"") )
  {
    ResamplerType::Pointer resampler = ResamplerType::New();

    if (nn)
    {
      resampler -> SetInterpolator(interpolator);
    }

    resampler -> SetInput( mgrad_std );
    resampler -> SetTransform( transform );
    resampler -> SetDefaultPixelValue(0);
    resampler -> UseReferenceImageOn();
    resampler -> SetReferenceImage( mgrad_std );

    resampler -> Update();
    mgrad_resampled = resampler -> GetOutput();
  }

  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer  imageWriter  = ImageWriterType::New();
  imageWriter -> SetInput( joiner2 -> GetOutput() );
  imageWriter -> SetFileName(  outputImageFile.c_str() );
  imageWriter -> Update();

  typedef itk::ImageFileWriter< Image3DType > Image3DWriterType;
  if ( strcmp(mgradName,"") )
  {
    Image3DWriterType::Pointer  image3DWriter  = Image3DWriterType::New();
    image3DWriter -> SetInput( mgrad_resampled );
    image3DWriter -> SetFileName(  mgradResampled );
    image3DWriter -> Update();
  }

  // Change gradient table

  gradientTable -> SetNumberOfGradients(numberOfFrames);
  gradientTable -> SetImage( stdImage );
  gradientTable -> SetTransform( transform );
  gradientTable -> RotateGradientsInWorldCoordinates();
  gradientTable -> SaveToFile( bvec_out.c_str());

  // Write b-values
  char clcopybval[255];
  sprintf(clcopybval,"cp %s %s",bval_in.c_str(),bval_out.c_str());
  system(clcopybval);


  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

