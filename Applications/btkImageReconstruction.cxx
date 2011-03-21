/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 24/01/2011
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

#include "btkLowToHighImageResolutionMethod.h"
#include "btkSliceBySliceRigidRegistration.h"
#include "btkResampleImageByInjectionFilter.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkDivideByConstantImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "itkEuler3DTransform.h"
#include "btkSliceBySliceTransform.h"
#include "itkTransformFileWriter.h"
#include "itkImage.h"
#include "CmdLine.h"

#include "itkImageMaskSpatialObject.h"
#include "btkImageIntersectionCalculator.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <stdio.h>

int main( int argc, char *argv[] )
{

  try {

  std::vector< std::string > input;
  std::vector< std::string > mask;
  std::vector< std::string > transform;
  std::vector< std::string > roi;
  std::vector< std::string > resampled;
  unsigned int itMax;
  double epsilon;
  double margin;

  const char *outImage = NULL;

  std::vector< int > x1, y1, z1, x2, y2, z2;

  // Parse arguments

  TCLAP::CmdLine cmd("Creates a high resolution image from a set of low resolution images", ' ', "Unversioned");

  TCLAP::MultiArg<std::string> inputArg("i","input","Image file",true,"string",cmd);
  TCLAP::MultiArg<std::string> maskArg("m","","Mask file",false,"string",cmd);
  TCLAP::MultiArg<std::string> transformArg("t","transform","Transform output file",false,"string",cmd);
  TCLAP::MultiArg<std::string> roiArg("","roi","roi file (written as mask)",false,"string",cmd);
  TCLAP::MultiArg<std::string> resampledArg("","ir","Resampled image with initial transform",false,"string",cmd);

  TCLAP::ValueArg<std::string> outArg("o","output","High resolution image",true,"none","string",cmd);
  TCLAP::ValueArg<unsigned int> iterArg("n","iter","Maximum number of iterations",false, 30,"unsigned int",cmd);
  TCLAP::ValueArg<double> epsilonArg("e","epsilon","Maximum number of iterations",false, 1e-4,"double",cmd);

  TCLAP::ValueArg<double> marginArg("","margin","Adds a margin to the reconstructed images"
  "to compensate for a small FOV in the reference. The value must be provided in millimeters.",false, 0.0,"double",cmd);

  TCLAP::SwitchArg  boxSwitchArg("","box","Use intersections for roi calculation",false);
  TCLAP::SwitchArg  maskSwitchArg("","mask","Use masks for roi calculation",false);
  TCLAP::SwitchArg  allSwitchArg("","all","Use the whole image FOV",false);

  // xor arguments for roi assessment

  std::vector<TCLAP::Arg*>  xorlist;
  xorlist.push_back(&boxSwitchArg);
  xorlist.push_back(&maskSwitchArg);
  xorlist.push_back(&allSwitchArg);

  cmd.xorAdd( xorlist );

  // Parse the argv array.
  cmd.parse( argc, argv );

  input = inputArg.getValue();
  mask = maskArg.getValue();
  outImage = outArg.getValue().c_str();
  transform = transformArg.getValue();
  roi = roiArg.getValue();
  resampled = resampledArg.getValue();
  itMax = iterArg.getValue();
  epsilon = epsilonArg.getValue();
  margin = marginArg.getValue();


  // typedefs

  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;

  typedef ImageType::RegionType               RegionType;
  typedef std::vector< RegionType >           RegionArrayType;

  typedef itk::ImageFileReader< ImageType  >  ImageReaderType;

  typedef itk::Image< unsigned char, Dimension >    ImageMaskType;
  typedef ImageMaskType::Pointer                    ImageMaskPointer;

  typedef itk::ImageFileReader< ImageMaskType >     MaskReaderType;
  typedef itk::ImageMaskSpatialObject< Dimension >  MaskType;
  typedef MaskType::Pointer  MaskPointer;

  typedef btk::LowToHighImageResolutionMethod<ImageType> LowToHighResFilterType;
  LowToHighResFilterType::Pointer lowToHighResFilter = LowToHighResFilterType::New();

  typedef btk::SliceBySliceRigidRegistration<ImageType> SliceBySliceRegistrationType;
  typedef SliceBySliceRegistrationType::Pointer RegistrationPointer;

  typedef btk::ResampleImageByInjectionFilter< ImageType, ImageType >  ResamplerType;

  typedef btk::SliceBySliceTransform< double, Dimension > TransformType;
  typedef TransformType::Pointer                          TransformPointer;

  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType,
                                                        ImageType > NCMetricType;

  typedef btk::ImageIntersectionCalculator<ImageType> IntersectionCalculatorType;
  IntersectionCalculatorType::Pointer intersectionCalculator = IntersectionCalculatorType::New();

  // Filter setup
  unsigned int numberOfImages = input.size();
  std::vector< ImagePointer >         images(numberOfImages);
  std::vector< ImageMaskPointer >     imageMasks(numberOfImages);
  std::vector< TransformPointer >     transforms(numberOfImages);
  std::vector< RegistrationPointer >  registration(numberOfImages);
  std::vector< MaskPointer >          masks(numberOfImages);
  std::vector< RegionType >           rois(numberOfImages);

  ImagePointer hrImage;
  ImagePointer hrImageOld;
  ImagePointer hrImageIni;

  lowToHighResFilter -> SetNumberOfImages(numberOfImages);
  lowToHighResFilter -> SetTargetImage( 0 );
  lowToHighResFilter -> SetMargin( margin );

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    // Read and set image

    ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input[i].c_str() );
    imageReader -> Update();
    images[i] = imageReader -> GetOutput();

    lowToHighResFilter -> SetImageArray(i, images[i] );

    if ( boxSwitchArg.isSet() )
    {
      intersectionCalculator -> AddImage( images[i] );
    }

  }

  if ( boxSwitchArg.isSet() )
  {
    intersectionCalculator -> Compute();
  }

  // Set roi according to the provided arguments

  for (unsigned int i=0; i<numberOfImages; i++)
  {

    if ( maskSwitchArg.isSet() )
    {
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader -> SetFileName( mask[i].c_str() );
      maskReader -> Update();
      imageMasks[i] = maskReader -> GetOutput();
      lowToHighResFilter -> SetImageMaskArray( i, imageMasks[i] );

      masks[i] = MaskType::New();
      masks[i] -> SetImage( imageMasks[i] );
      rois[i] = masks[i] -> GetAxisAlignedBoundingBoxRegion();

      lowToHighResFilter -> SetRegionArray( i, rois[i] );

    } else if ( boxSwitchArg.isSet() )
      {
        imageMasks[i] = intersectionCalculator -> GetImageMask(i);
        lowToHighResFilter -> SetImageMaskArray( i, imageMasks[i] );

        masks[i] = MaskType::New();
        masks[i] -> SetImage( imageMasks[i] );
        rois[i] = masks[i] -> GetAxisAlignedBoundingBoxRegion();
        lowToHighResFilter -> SetRegionArray( i, rois[i] );

      } else if ( allSwitchArg.isSet() )
        {
          rois[i] = images[i] -> GetLargestPossibleRegion();

          lowToHighResFilter -> SetRegionArray( i, rois[i] );
        }
  }

  // Start registration
  try
    {
    lowToHighResFilter->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Write resampled images
  if ( resampled.size() > 0 )
  {
    for (unsigned int i=0; i<numberOfImages; i++)
    {
      lowToHighResFilter -> WriteResampledImages( i, resampled[i].c_str() );
    }
  }

  // Register slice by slice

  hrImageIni = lowToHighResFilter->GetHighResolutionImage();

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    transforms[i] = TransformType::New();
    transforms[i] -> SetImage( images[i] );
    transforms[i] -> Initialize( lowToHighResFilter -> GetTransformArray(i) );
  }

  unsigned int im = numberOfImages;
  float previousMetric = 0.0;
  float currentMetric = 0.0;

  for(unsigned int it=1; it <= itMax; it++)
  {
    std::cout << "Iteration " << it << std::endl; std::cout.flush();

    // Start registration

    #pragma omp parallel for private(im) schedule(dynamic)

    for (im=0; im<numberOfImages; im++)
    {
      std::cout << "Registering image " << im << " ... "; std::cout.flush();

      registration[im] = SliceBySliceRegistrationType::New();
      registration[im] -> SetFixedImage( images[im] );
      registration[im] -> SetMovingImage( hrImageIni );
      registration[im] -> SetImageMask( imageMasks[im] );
      registration[im] -> SetTransform( transforms[im] );

      try
        {
        registration[im] -> StartRegistration();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
  //      return EXIT_FAILURE;
        }

      transforms[im] = registration[im] -> GetTransform();

      std::cout << "done. "; std::cout.flush();

    }

    std::cout << std::endl; std::cout.flush();

    // Inject images

    std::cout << "Injecting images ... "; std::cout.flush();

    ResamplerType::Pointer resampler = ResamplerType::New();

    for (unsigned int i=0; i<numberOfImages; i++)
    {
      resampler -> AddInput( images[i] );
      resampler -> AddRegion( rois[i] );

      unsigned int nslices = transforms[i] -> GetNumberOfSlices();

      for(unsigned int j=0; j< nslices; j++)
      {
        resampler -> SetTransform(i, j, transforms[i] -> GetSliceTransform(j) ) ;
      }
    }

    resampler -> UseReferenceImageOn();
    resampler -> SetReferenceImage( hrImageIni );
    resampler -> Update();

    if (it == 1)
      hrImageOld = hrImageIni;
    else
      hrImageOld = hrImage;

    hrImage = resampler -> GetOutput() ;

    std::cout << "done. " << std::endl; std::cout.flush();

    // compute error

    typedef itk::MinimumMaximumImageCalculator<ImageType> CalculatorType;

    CalculatorType::Pointer calculatorA = CalculatorType::New();
    calculatorA -> SetImage( hrImageOld );
    calculatorA -> ComputeMaximum();

    typedef itk::DivideByConstantImageFilter<ImageType, float, ImageType> NormalizerType;

    NormalizerType::Pointer normalizerA = NormalizerType::New();
    normalizerA -> SetInput( hrImageOld );
    normalizerA -> SetConstant( calculatorA -> GetMaximum() );
    normalizerA -> Update();

    NormalizerType::Pointer normalizerB = NormalizerType::New();
    normalizerB -> SetInput( hrImage );
    normalizerB -> SetConstant( calculatorA -> GetMaximum() );
    normalizerB -> Update();

    ImagePointer imageA = normalizerA -> GetOutput();
    ImagePointer imageB = normalizerB -> GetOutput();

    typedef itk::Euler3DTransform< double > EulerTransformType;
    EulerTransformType::Pointer identity = EulerTransformType::New();
    identity -> SetIdentity();

    typedef itk::LinearInterpolateImageFunction<
                                      ImageType,
                                      double>     InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();

    NCMetricType::Pointer nc = NCMetricType::New();
    nc -> SetFixedImage(  imageA );
    nc -> SetMovingImage( imageB );
    if ( maskSwitchArg.isSet() )
    {
      nc -> SetFixedImageRegion( rois[0] );
      nc -> SetFixedImageMask( masks[0] );
    }
    else
      nc -> SetFixedImageRegion( imageA -> GetLargestPossibleRegion() );

    nc -> SetTransform( identity );
    nc -> SetInterpolator( interpolator );
    nc -> Initialize();

    previousMetric = currentMetric;
    currentMetric = - nc -> GetValue( identity -> GetParameters() );

    double delta = 0.0;

    if (it >= 2)
      delta = (currentMetric - previousMetric) / previousMetric;
    else
      delta = 1;

    if (delta < epsilon) break;

  }


  // Write HR image

  typedef itk::ImageFileWriter< ImageType >  WriterType;

  WriterType::Pointer writer =  WriterType::New();
  writer-> SetFileName( outImage );
  writer-> SetInput( hrImage );
  writer-> Update();

  // Write transforms

  typedef itk::TransformFileWriter TransformWriterType;

  if ( transform.size() > 0 )
  {
    for (unsigned int i=0; i<numberOfImages; i++)
    {
      TransformWriterType::Pointer transformWriter = TransformWriterType::New();
      transformWriter -> SetInput( transforms[i] );
      transformWriter -> SetFileName ( transform[i].c_str() );

      try
      {
        std::cout << "Writing " << transform[i].c_str() << " ... " ; std::cout.flush();
        transformWriter -> Update();
        std::cout << " done! " << std::endl;
      }
      catch ( itk::ExceptionObject & excp )
      {
        std::cerr << "Error while saving transform" << std::endl;
        std::cerr << excp << std::endl;
        std::cout << "[FAILED]" << std::endl;
        throw excp;
      }

    }
  }

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

