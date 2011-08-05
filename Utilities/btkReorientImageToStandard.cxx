/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 31/01/2011
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

#include <tclap/CmdLine.h>
#include "vnl/vnl_cross.h"

#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkOrientImageFilter.h"

#include "vnl/vnl_inverse.h"

// Local includes
#include "btkLandmarksFileReader.h"


int main( int argc, char *argv[] )
{

  try {

  // Parse arguments

  TCLAP::CmdLine cmd("Reorient an image to standard orientation", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Reoriented image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> lmksArg("l","landmarks","Landmarks file",true,"","string",cmd);
  TCLAP::SwitchArg nnSwitch("","nn","Nearest Neighbor interpolation", cmd, false);
  TCLAP::SwitchArg fovSwitch("e","extend-fov","Extends the FOV to avoid image cropping", cmd, false);

  cmd.parse( argc, argv );

  const char *inputImageFile  = inputArg.getValue().c_str();
  const char *outputImageFile = outputArg.getValue().c_str();
  const char *lmksFile   = lmksArg.getValue().c_str();

  bool nn         = nnSwitch.getValue();

  const    unsigned int    Dimension = 3;

  // Read image

  typedef itk::Image< short, Dimension >  ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
  imageReader -> SetFileName(  inputImageFile );
  imageReader -> Update();
  ImageType::Pointer image = imageReader->GetOutput();

  // Read landmarks

  btk::btkLandmarksFileReader::Pointer landRead = btk::btkLandmarksFileReader::New();
  landRead->SetInputFileName(lmksFile);
  landRead->Update();

  double *lpt_ras = landRead->GetOutputLPT();
  double *rpt_ras = landRead->GetOutputRPT();
  double *apt_ras = landRead->GetOutputAPT();
  double *ppt_ras = landRead->GetOutputPPT();

  // Convert points to lps coordinates

  ImageType::PointType leftPoint, rightPoint, anteriorPoint, posteriorPoint;

  leftPoint[0]      = -lpt_ras[0]; leftPoint[1]      = -lpt_ras[1]; leftPoint[2]      = lpt_ras[2];
  rightPoint[0]     = -rpt_ras[0]; rightPoint[1]     = -rpt_ras[1]; rightPoint[2]     = rpt_ras[2];
  anteriorPoint[0]  = -apt_ras[0]; anteriorPoint[1]  = -apt_ras[1]; anteriorPoint[2]  = apt_ras[2];
  posteriorPoint[0] = -ppt_ras[0]; posteriorPoint[1] = -ppt_ras[1]; posteriorPoint[2] = ppt_ras[2];

  // Recover original direction matrix

  vnl_matrix<double> oriDirection(3,3);
  oriDirection = image -> GetDirection().GetVnlMatrix();

  // Reorient mean gradient to RAI

  typedef itk::OrientImageFilter< ImageType, ImageType >  FilterType;

  // Reorient sequence to RAI

  FilterType::Pointer filter = FilterType::New();
  filter -> SetInput( image  );
  filter -> UseImageDirectionOn();
  filter -> SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
  filter -> Update();

  ImageType::Pointer image_reo = filter -> GetOutput();

  vnl_matrix<double> reoDirection(3,3);
  reoDirection = image_reo -> GetDirection().GetVnlMatrix();

  vnl_matrix<double> reoDirectionTrans(3,3);
  reoDirectionTrans = reoDirection.transpose();

  // Create the image in standard orientation
  // p = I*S*idx + O_center

  ImageType::Pointer stdImage = image_reo;

  ImageType::DirectionType stdImageDirection  = stdImage -> GetDirection();
  stdImageDirection.SetIdentity();
  stdImage -> SetDirection (stdImageDirection);

  ImageType::SpacingType stdSpacing  = stdImage -> GetSpacing();
  ImageType::SizeType 	 stdSize     = stdImage -> GetLargestPossibleRegion().GetSize();
  ImageType::IndexType   stdIndex    = stdImage -> GetLargestPossibleRegion().GetIndex();

  ImageType::PointType stdOrigin;

  for (unsigned int i=0; i<3; i++)
    stdOrigin[i] = (1.0-stdSize[i])*0.5*stdSpacing[i];

  stdImage -> SetOrigin( stdOrigin );

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

  typedef itk::Euler3DTransform<double> TransformType;
  TransformType::Pointer transform = TransformType::New();

  typedef itk::ContinuousIndex< double, 3 > ContinuousIndexType;
  ContinuousIndexType centerIndex;

  centerIndex[0] = (stdSize[0] - 1.0) / 2.0;
  centerIndex[1] = (stdSize[1] - 1.0) / 2.0;
  centerIndex[2] = (stdSize[2] - 1.0) / 2.0;

  ImageType::PointType centerPoint;

  stdImage -> TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);

  transform -> SetCenter( centerPoint );
  transform -> SetRotationMatrix( NQ );

//  TODO Cannot constrain the rotation because sometimes it is necessary to rotate in xy
//  Something I could do is to round the rotations in xy to n*pi*/2 ...
//  transform -> SetRotation( 0.0, 0.0, transform -> GetAngleZ() );

  // Extend FOV according to the transformation

  ImageType::SizeType  resamplingSize;
  ImageType::IndexType lowerIndex; lowerIndex.Fill(0);
  ImageType::IndexType upperIndex; upperIndex.Fill(0);

  if ( fovSwitch.isSet() )
  {
    typedef itk::ImageRegionIteratorWithIndex< ImageType > ImageIteratorType;
    ImageIteratorType stdImageIt( stdImage, stdImage -> GetLargestPossibleRegion() );

    ImageType::IndexType index;
    ImageType::IndexType transformedIndex;

    ImageType::PointType point;
    ImageType::PointType transformedPoint;

    // Even when only the corners are necessary, the code is shorter in this way
    // (and it's not very expensive to compute it)

    for(stdImageIt.GoToBegin(); !stdImageIt.IsAtEnd(); ++stdImageIt)
    {
      index = stdImageIt.GetIndex();
      stdImage -> TransformIndexToPhysicalPoint(index, point);
      transformedPoint = transform -> TransformPoint(point);
      stdImage -> TransformPhysicalPointToIndex(transformedPoint, transformedIndex);

      for (unsigned int i = 0; i<3; i++)
        if ( transformedIndex[i] < lowerIndex[i] ) lowerIndex[i] = transformedIndex[i];
        else if
           ( transformedIndex[i] > upperIndex[i] ) upperIndex[i] = transformedIndex[i];

    }

    resamplingSize[0] = upperIndex[0] - lowerIndex[0] + 1;
    resamplingSize[1] = upperIndex[1] - lowerIndex[1] + 1;
    resamplingSize[2] = upperIndex[2] - lowerIndex[2] + 1;

  } else
  {
    resamplingSize = stdImage -> GetLargestPossibleRegion().GetSize();
  }

  ImageType::PointType resamplingOrigin;
  stdImage -> TransformIndexToPhysicalPoint(lowerIndex, resamplingOrigin );

  transform -> SetRotationMatrix( NQ.transpose() );

  // Resampling

  typedef itk::ResampleImageFilter< ImageType, ImageType >  ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();

  typedef itk::NearestNeighborInterpolateImageFunction< ImageType >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  if (nn) resampler -> SetInterpolator(interpolator);

  resampler -> SetInput( stdImage );
  resampler -> SetTransform(transform);
  resampler -> SetDefaultPixelValue(0);
  resampler -> SetSize( resamplingSize );
  resampler -> SetOutputOrigin( resamplingOrigin );
  resampler -> SetOutputDirection( stdImage -> GetDirection()  );
  resampler -> SetOutputSpacing( stdImage -> GetSpacing() );

  resampler -> Update();

  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer  imageWriter  = ImageWriterType::New();
  imageWriter -> SetInput( resampler -> GetOutput() );
  imageWriter -> SetFileName(  outputImageFile );
  imageWriter -> Update();

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

