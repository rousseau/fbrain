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

#include "CmdLine.h"
#include "vnl/vnl_cross.h"

#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "vnl/vnl_inverse.h"

// Local includes
#include "btkLandmarksFileReader.h"


int main( int argc, char *argv[] )
{

  try {

  // Parse arguments
  const char *inputImageFile = NULL, *outputImageFile = NULL, *lmksFile = NULL;

  TCLAP::CmdLine cmd("Change the image orientation", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Reoriented image",true,"none","string",cmd);
  TCLAP::ValueArg<std::string> lmksArg("l","landmarks","Landmarks file",true,"none","string",cmd);
  TCLAP::SwitchArg nnSwitch("","nn","Nearest Neighbor interpolation", cmd, false);

  cmd.parse( argc, argv );

  inputImageFile  = inputArg.getValue().c_str();
  outputImageFile = outputArg.getValue().c_str();
  lmksFile        = lmksArg.getValue().c_str();
  bool nn = nnSwitch.getValue();

  const    unsigned int    Dimension = 3;

  // Read image

  typedef itk::Image< short, Dimension >  ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
  imageReader -> SetFileName(  inputImageFile );
  imageReader -> Update();
  ImageType::Pointer image = imageReader->GetOutput();


  // Read landmarks
/*
  FILE* fr;
  fr = fopen( lmksFile, "r" );

  double lpt_ras[3], rpt_ras[3], apt_ras[3], ppt_ras[3];

  fscanf( fr, "%lf %lf %lf\n", &lpt_ras[0], &lpt_ras[1], &lpt_ras[2]);
  fscanf( fr, "%lf %lf %lf\n", &rpt_ras[0], &rpt_ras[1], &rpt_ras[2]);
  fscanf( fr, "%lf %lf %lf\n", &ppt_ras[0], &ppt_ras[1], &ppt_ras[2]);
  fscanf( fr, "%lf %lf %lf\n", &apt_ras[0], &apt_ras[1], &apt_ras[2]);

  fclose (fr);
//*/

    btk::btkLandmarksFileReader::Pointer landRead = btk::btkLandmarksFileReader::New();
    landRead->SetInputFileName(lmksFile);
    landRead->Update();

    double *lpt_ras = landRead->GetOutputLPT();
    double *rpt_ras = landRead->GetOutputRPT();
    double *apt_ras = landRead->GetOutputPPT();
    double *ppt_ras = landRead->GetOutputAPT();

//std::cout << "lpt_ras = " << lpt_ras[0] << "," << lpt_ras[1] << "," << lpt_ras[2] << " | rpt_ras = " << rpt_ras[0] << "," << rpt_ras[1] << "," << rpt_ras[2] << " | apt_ras = " << apt_ras[0] << "," << apt_ras[1] << "," << apt_ras[2] << " | ppt_ras = " << ppt_ras[0] << "," << ppt_ras[1] << "," << ppt_ras[2] << std::endl;


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
//      std::cout << "Polar decomposition used "
//                << ni << " iterations " << std::endl;
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

  typedef itk::ContinuousIndex< double, Dimension > ContinuousIndexType;
  ContinuousIndexType centerIndex;

  centerIndex[0] = (size[0] - 1.0)/2.0;
  centerIndex[1] = (size[1] - 1.0)/2.0;
  centerIndex[2] = (size[2] - 1.0)/2.0;

  ImageType::PointType centerPoint;
  image -> TransformContinuousIndexToPhysicalPoint(centerIndex, centerPoint);
//  std::cout << "Center point = " << centerPoint << std::endl;

  transform -> SetCenter( centerPoint );
  transform -> SetRotationMatrix( NQ );
  transform -> SetRotationMatrix( NQ.transpose() );

//  std::cout << "Rotation matrix corrected = " << std::endl << NQ.transpose() << std::endl;
//  transform -> SetRotation( 0.0, 0.0, transform -> GetAngleZ() );

  // Resampling

  typedef itk::ResampleImageFilter< ImageType, ImageType >  ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();

  typedef itk::NearestNeighborInterpolateImageFunction< ImageType >  InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  resampler -> SetInput( image );
  resampler -> SetReferenceImage( image );
  resampler -> SetUseReferenceImage( true );
  resampler -> SetTransform(transform);
  resampler -> SetDefaultPixelValue(0);

  if (nn)
  {
    resampler -> SetInterpolator(interpolator);
  }


  resampler -> Update();

  ImageType::Pointer reorientedImage = resampler -> GetOutput();


  typedef itk::ImageFileWriter< ImageType > ImageWriterType;
  ImageWriterType::Pointer  imageWriter  = ImageWriterType::New();
  imageWriter->SetInput(  reorientedImage );
  imageWriter->SetFileName(  outputImageFile );
  imageWriter->Update();

  // Transform landmarks

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

