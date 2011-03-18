/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/03/2011
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

#include "tclap/CmdLine.h"

#include "itkImage.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"

#include "itkExtractImageFilter.h"

#include "itkResampleImageFilter.h"

#include "itkImageMaskSpatialObject.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkJoinSeriesImageFilter.h"

#include "btkDiffusionGradientTable.h"


//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer -> GetCurrentStepLength() << "   ";
      std::cout << optimizer->GetValue() << "   ";

      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};


int main( int argc, char *argv[] )
{

  const char *inputName = NULL, *referenceName = NULL, *outputName = NULL;
  const char *toutName = NULL, *invToutName = NULL, *maskName = NULL;
  const char *gTableName = NULL, *cTableName = NULL;

  TCLAP::CmdLine cmd("Registers diffusion to anatomical data.", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Diffusion sequence",true,"","string",cmd);
  TCLAP::ValueArg<std::string> gTableArg("g","gtable","Gradient table",true,"","string",cmd);
  TCLAP::ValueArg<std::string> referenceArg("r","reference","Anatomical image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Registered diffusion sequence",true,"","string",cmd);
  TCLAP::ValueArg<std::string> cTableArg("c","ctable","Corrected table",true,"","string",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","mask","Mask for the reference image",false,"","string",cmd);

  TCLAP::ValueArg<unsigned int> x1Arg("","x1","x min of ROI in diffusion",false,0,"int",cmd);
  TCLAP::ValueArg<unsigned int> x2Arg("","x2","x max of ROI in diffusion",false,0,"int",cmd);

  TCLAP::ValueArg<unsigned int> y1Arg("","y1","y min of ROI in diffusion",false,0,"int",cmd);
  TCLAP::ValueArg<unsigned int> y2Arg("","y2","y max of ROI in diffusion",false,0,"int",cmd);

  TCLAP::ValueArg<unsigned int> z1Arg("","z1","z min of ROI in diffusion",false,0,"int",cmd);
  TCLAP::ValueArg<unsigned int> z2Arg("","z2","z max of ROI in diffusion",false,0,"int",cmd);

  TCLAP::ValueArg<std::string> toutArg("","tout","Estimated transformation",false,"","string",cmd);
  TCLAP::ValueArg<std::string> invToutArg("","itout","Inverse of estimated transformation",false,"","string",cmd);


  // Parse the argv array.
  cmd.parse( argc, argv );

  inputName     = inputArg.getValue().c_str();
  referenceName = referenceArg.getValue().c_str();
  outputName    = outputArg.getValue().c_str();
  maskName      = maskArg.getValue().c_str();

  unsigned int x1 = x1Arg.getValue();
  unsigned int x2 = x2Arg.getValue();

  unsigned int y1 = y1Arg.getValue();
  unsigned int y2 = y2Arg.getValue();

  unsigned int z1 = z1Arg.getValue();
  unsigned int z2 = z2Arg.getValue();

  bool roiIsSet = false;

  if ( (x1>0) || (x2>0) || (y1>0)  || (y2>0)  || (z1>0)  || (z2>0) )
    roiIsSet = true;

  toutName = toutArg.getValue().c_str();
  invToutName = invToutArg.getValue().c_str();

  gTableName = gTableArg.getValue().c_str();
  cTableName = cTableArg.getValue().c_str();

  // Typedefs

  typedef  short  PixelType;

  typedef itk::Image< PixelType, 3 >  ImageType;
  typedef itk::Image< PixelType, 4 >  SequenceType;

  typedef itk::AffineTransform< double, 3 > TransformType;
  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef itk::LinearInterpolateImageFunction<
                                    ImageType,
                                    double             > InterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction<
                                    ImageType,
                                    double             > NNInterpolatorType;
  typedef itk::ImageRegistrationMethod<
                                    ImageType,
                                    ImageType    > RegistrationType;

  typedef itk::MattesMutualInformationImageToImageMetric<
                                          ImageType,
                                          ImageType >    MetricType;


  typedef itk::Image< unsigned char, 3 >    ImageMaskType;
  typedef itk::ImageFileReader< ImageMaskType >     MaskReaderType;
  typedef itk::ImageMaskSpatialObject< 3 >  MaskType;


  // Registration object

  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  NNInterpolatorType::Pointer nn_interpolator  = NNInterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );

  MetricType::Pointer metric = MetricType::New();
  registration->SetMetric( metric  );

  unsigned int numberOfBins = 24;

  metric->SetNumberOfHistogramBins( numberOfBins );
  metric->UseAllPixelsOn();

  // Read input images (sequence and reference)

  typedef itk::ImageFileReader< ImageType  >   ImageReaderType;
  typedef itk::ImageFileReader< SequenceType > SequenceReaderType;

  ImageReaderType::Pointer    imageReader    = ImageReaderType::New();
  SequenceReaderType::Pointer sequenceReader = SequenceReaderType::New();

  imageReader    -> SetFileName( referenceName );
  sequenceReader -> SetFileName( inputName );

  imageReader    -> Update();
  sequenceReader -> Update();

  ImageType::Pointer    reference = imageReader -> GetOutput();
  SequenceType::Pointer input     = sequenceReader -> GetOutput();

  // Extract B0 from sequence

  typedef itk::ExtractImageFilter< SequenceType, ImageType > ExtractorType;
  ExtractorType::Pointer extractor = ExtractorType::New();

  extractor -> SetInput( input );

  SequenceType::SizeType  size  = input -> GetLargestPossibleRegion().GetSize();
  SequenceType::IndexType index = input -> GetLargestPossibleRegion().GetIndex();

  unsigned int inputLength = size[3];

  size[3]  = 0;

  SequenceType::RegionType region;
  region.SetIndex(index);
  region.SetSize(size);

  extractor -> SetExtractionRegion( region );
  extractor -> Update();

  ImageType::Pointer b0 = extractor -> GetOutput();

  registration -> SetFixedImage( b0 );
  registration -> SetMovingImage( reference );

  // Fixed image region

  ImageType::RegionType fixedImageRegion;
  ImageType::IndexType  fixedImageRegionIndex;
  ImageType::SizeType   fixedImageRegionSize;

  if ( strcmp(maskName,"") )
  {
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader -> SetFileName( maskName );
    maskReader -> Update();

    MaskType::Pointer mask = MaskType::New();
    mask -> SetImage( maskReader -> GetOutput() );

    fixedImageRegion = mask -> GetAxisAlignedBoundingBoxRegion();
    metric -> SetFixedImageMask( mask );
  } else if ( roiIsSet )
    {
      fixedImageRegionIndex[0] = x1; fixedImageRegionIndex[1]= y1; fixedImageRegionIndex[2]= z1;
      fixedImageRegionSize[0]  = x2 - x1 + 1; fixedImageRegionSize[1] = y2 - y1 + 1; fixedImageRegionSize[2] = z2 - z1 + 1;

      fixedImageRegion.SetIndex(fixedImageRegionIndex);
      fixedImageRegion.SetSize(fixedImageRegionSize);
    } else
      {
        fixedImageRegion = reference -> GetLargestPossibleRegion();
      }

  fixedImageRegionIndex = fixedImageRegion.GetIndex();
  fixedImageRegionSize  = fixedImageRegion.GetSize();

  registration -> SetFixedImageRegion ( fixedImageRegion );

  // Center of rotation

  ImageType::IndexType  centerIndex;
  ImageType::PointType  centerPoint;

  centerIndex[0] = fixedImageRegionIndex[0] + fixedImageRegionSize[0] / 2.0;
  centerIndex[1] = fixedImageRegionIndex[1] + fixedImageRegionSize[1] / 2.0;
  centerIndex[2] = fixedImageRegionIndex[2] + fixedImageRegionSize[2] / 2.0;

  b0 -> TransformIndexToPhysicalPoint(centerIndex, centerPoint);

  transform -> SetIdentity();
  transform -> SetCenter( centerPoint );


  typedef RegistrationType::ParametersType ParametersType;

  ParametersType initialParameters( transform->GetNumberOfParameters() );
  initialParameters = transform->GetParameters();

  registration->SetInitialTransformParameters( initialParameters );

  std::cout << "Initial parameters = " << initialParameters << std::endl;

  optimizer->MinimizeOn();
  optimizer->SetMaximumStepLength( 0.2 );
  optimizer->SetMinimumStepLength( 0.001 );
  optimizer->SetNumberOfIterations( 300 );

  typedef OptimizerType::ScalesType       OptimizerScalesType;
  OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  1.0;
  optimizerScales[5] =  1.0;
  optimizerScales[6] =  1.0;
  optimizerScales[7] =  1.0;
  optimizerScales[8] =  1.0;
  optimizerScales[9] =  1/1000.0;
  optimizerScales[10] =  1/1000.0;
  optimizerScales[11] =  1/1000.0;

  optimizer->SetScales( optimizerScales );

  // Create the Command observer and register it with the optimizer.
  //
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );


  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  ParametersType finalParameters = registration->GetLastTransformParameters();

  TransformType::Pointer finalTransform = TransformType::New();

  finalTransform->SetParameters( finalParameters );
  finalTransform->SetCenter( centerPoint );

  // Calculate inverse transform (anat -> diffusion)

  TransformType::Pointer inverseTransform = TransformType::New();
  inverseTransform -> SetIdentity();
  inverseTransform -> SetCenter( finalTransform -> GetCenter() );
  inverseTransform -> SetParameters( finalTransform -> GetParameters() );
  inverseTransform -> GetInverse( inverseTransform );

  // Resample sequence

  SequenceType::SpacingType  spacing  = input -> GetSpacing();

  typedef itk::JoinSeriesImageFilter< ImageType, SequenceType > JoinerType;
  JoinerType::Pointer joiner = JoinerType::New();
  joiner -> SetOrigin( 0.0 );
  joiner -> SetSpacing( spacing[3] );

  typedef itk::ResampleImageFilter<
                            ImageType,
                            ImageType >    ResampleFilterType;

  for (unsigned int i = 0; i < inputLength; i++)
  {
     index[3] = i;

     SequenceType::RegionType desiredRegion;
     desiredRegion.SetSize(  size  );
     desiredRegion.SetIndex( index );

     extractor -> SetExtractionRegion( desiredRegion );

     ResampleFilterType::Pointer resample = ResampleFilterType::New();

     resample -> SetTransform( inverseTransform );
     resample -> SetInput( extractor -> GetOutput() );
     resample -> SetUseReferenceImage( true );
     resample -> SetReferenceImage( reference );
     resample -> SetDefaultPixelValue( 0 );
     resample -> Update();

     joiner -> SetInput( i, resample -> GetOutput() );
  }

  joiner -> Update();

  // Write resampled sequence

  typedef itk::ImageFileWriter< SequenceType >  WriterType;

  WriterType::Pointer writer =  WriterType::New();
  writer->SetFileName( outputName );
  writer->SetInput( joiner -> GetOutput() );

  if (strcmp(outputName,"")) writer->Update();

  // Change and write gradient table

  typedef btk::DiffusionGradientTable< SequenceType > GradientTableType;
  GradientTableType::Pointer gradientTable = GradientTableType::New();

  vnl_matrix<double> R(3,3);
  R.set_identity();

  R = inverseTransform -> GetMatrix().GetVnlMatrix();

  vnl_matrix<double> PQ = R;
  vnl_matrix<double> NQ = R;
  vnl_matrix<double> PQNQDiff;

  for(unsigned int ni = 0; ni < 100; ni++ )
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

  typedef itk::Euler3DTransform<double> EulerTransformType;
  EulerTransformType::Pointer eulerTransform = EulerTransformType::New();

  eulerTransform -> SetRotationMatrix( NQ );

  gradientTable -> SetNumberOfGradients( inputLength );
  gradientTable -> SetImage( input );
  gradientTable -> SetTransform( eulerTransform );
  gradientTable -> LoadFromFile( gTableName );
  gradientTable -> RotateGradientsInWorldCoordinates();
  gradientTable -> SaveToFile( cTableName);

  // Write transformation if required

  typedef itk::TransformFileWriter TransformWriterType;

  if (strcmp(toutName,""))
  {
    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter -> SetFileName ( toutName );
    transformWriter -> AddTransform( finalTransform );

    try
    {
      transformWriter -> Update();
    }
    catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while saving the transforms" << std::endl;
      std::cerr << excp << std::endl;
      std::cout << "[FAILED]" << std::endl;
    }
  }

  // Write inverse transformation if required

  if (strcmp(invToutName,""))
  {

    TransformWriterType::Pointer transformWriter = TransformWriterType::New();
    transformWriter -> SetFileName ( invToutName );
    transformWriter -> AddTransform( inverseTransform );

    try
    {
      transformWriter -> Update();
    }
    catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while saving the transforms" << std::endl;
      std::cerr << excp << std::endl;
      std::cout << "[FAILED]" << std::endl;
    }
  }

  return EXIT_SUCCESS;
}

