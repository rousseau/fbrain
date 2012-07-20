/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/03/2011
  Author(s): Estanislao Oubel (oubel@unistra.fr)
             Julien Pontabry (pontabry@unistra.fr)

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


// STL includes
#include "cstdlib"
#include "iostream"
#include "string"
#include "sstream"

// TCLAP includes
#include <tclap/CmdLine.h>

// ITK includes
#include "itkImage.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"

#include "itkExtractImageFilter.h"

#include "itkResampleImageFilter.h"

#include "itkImageMaskSpatialObject.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkJoinSeriesImageFilter.h"

#include "itkAffineTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"

// Local includes
#include "btkFileNameTools.h"
#include "btkFileHelper.h"
#include "btkImageHelper.h"
#include "btkCommandIterationUpdate.h"
#include "btkMacro.h"
#include "btkDiffusionGradientTable.h"


// Image and sequence definitions
typedef short PixelType;
typedef itk::Image< PixelType,3 >               Image;
typedef itk::Image< PixelType,4 >               Sequence;
typedef itk::Image< unsigned char,3 >           ImageMask;
typedef itk::ImageMaskSpatialObject< 3 >        Mask;
typedef btk::DiffusionGradientTable< Sequence > GradientTable;

// Registration objects definitions
typedef double PrecisionType;
typedef itk::MatrixOffsetTransformBase< PrecisionType,3,3 >           Transform;
typedef itk::AffineTransform< PrecisionType,3 >                       AffineTransform;
typedef itk::TransformFileReader                                      AffineTransformFileReader;
typedef itk::TransformFileWriter                                      AffineTransformFileWriter;
typedef itk::Euler3DTransform< PrecisionType >                        EulerTransform;
typedef itk::RegularStepGradientDescentOptimizer                      RegularStepGradientDescentOptimizer;
typedef itk::ImageRegistrationMethod< Image,Image >                   RegistrationMethod;
typedef RegistrationMethod::ParametersType                            RegistrationParameters;
typedef itk::MattesMutualInformationImageToImageMetric< Image,Image > MattesMutualInformationMetric;

// Interpolators definitions
typedef itk::LinearInterpolateImageFunction< Image,PrecisionType >                LinearInterpolator;
typedef itk::NearestNeighborInterpolateImageFunction< Image,PrecisionType >       NearestNeighborInterpolator;
typedef itk::BSplineInterpolateImageFunction< Image,PrecisionType,PrecisionType > BSplineInterpolator;

// Filters definitions
typedef itk::ExtractImageFilter< Sequence,Image >    SequenceExtractor;
typedef itk::JoinSeriesImageFilter< Image,Sequence > ImageSequenceJoiner;
typedef itk::ResampleImageFilter< Image,Image >      ImageResampler;


int main(int argc, char *argv[])
{
    try
    {
        //
        // Parse program's arguments
        //

        // Define command line object for program
        TCLAP::CmdLine cmd("Registers diffusion sequence to anatomical data.", ' ', "1.0");

        // Define arguments
        TCLAP::ValueArg<std::string>  inputSequenceFileNameArg("i", "input", "Input diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputSequenceFileNameArg("o", "output", "Registred output diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> referenceImageFileNameArg("r", "reference", "Anatomical image, reference for registration", true, "", "string", cmd);

        TCLAP::ValueArg<std::string> maskImageFileNameArg("m", "mask", "Mask image for baseline image of diffusion sequence", false, "", "string", cmd);
        TCLAP::ValueArg<std::string> transformFileNameArg("t", "transformation", "Transformation from diffusion sequence to anatomical reference", false, "", "string", cmd);
        TCLAP::ValueArg<std::string> initialTransformFileNameArg("", "initial_transform", "Initial affine transformation to apply", false, "", "string", cmd);

        TCLAP::SwitchArg verboseModeArg("v", "verbose", "Verbose mode", cmd, false);

        // Parse command line
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string  inputSequenceFileName = inputSequenceFileNameArg.getValue();
        std::string outputSequenceFileName = outputSequenceFileNameArg.getValue();
        std::string referenceImageFileName = referenceImageFileNameArg.getValue();

        std::string        maskImageFileName = maskImageFileNameArg.getValue();
        std::string        transformFileName = transformFileNameArg.getValue();
        std::string initialTransformFileName = initialTransformFileNameArg.getValue();

        bool verboseMode = verboseModeArg.getValue();


        //
        // Preprocessing
        //

        // Register the transformation
        itk::TransformFactory< Transform >::RegisterTransform();

        // Set B-values and gradients' filenames according to diffusion sequences' filenames
        std::string     inputSequenceFileNameRadix = btk::GetRadixOf(inputSequenceFileName);
        std::string inputSequenceGradientsFileName = inputSequenceFileNameRadix + ".bvec";
        std::string inputSequenceBValuesFileName = inputSequenceFileNameRadix + ".bval";

        std::string     outputSequenceFileNameRadix = btk::GetRadixOf(outputSequenceFileName);
        std::string outputSequenceGradientsFileName = outputSequenceFileNameRadix + ".bvec";
        std::string   outputSequenceBValuesFileName = outputSequenceFileNameRadix + ".bval";

        // Read transformation if file exist
        bool transformProvided = false;
        AffineTransform::Pointer finalTransform;

        if(!transformFileName.empty())
        {
            if(btk::FileHelper::FileExist(transformFileName))
            {
                transformProvided = true;

                // Read transformation file (ITK style)
                AffineTransformFileReader::Pointer reader = AffineTransformFileReader::New();
                reader->SetFileName(transformFileName);
                reader->Update();
                finalTransform = static_cast< AffineTransform* >(reader->GetTransformList()->front().GetPointer());
            }
        }

        // Read initial transform if provided (and file exist)
        // TODO test if file exists
        AffineTransform::Pointer initialTransform = AffineTransform::New();
        initialTransform->SetIdentity();

        if(!initialTransformFileName.empty())
        {
            // Read transformation file (ITK style)
            AffineTransformFileReader::Pointer reader = AffineTransformFileReader::New();
            reader->SetFileName(initialTransformFileName);
            reader->Update();
            initialTransform = static_cast< AffineTransform* >(reader->GetTransformList()->front().GetPointer());
        }

        // TODO Test if files do exist


        // Read input images
        Image::Pointer   referenceImage = btk::ImageHelper< Image >::ReadImage(referenceImageFileName);
        Sequence::Pointer inputSequence = btk::ImageHelper< Sequence >::ReadImage(inputSequenceFileName);

        // Read mask if any
        Mask::Pointer mask;
        if(!maskImageFileName.empty())
        {
            ImageMask::Pointer imageMask = btk::ImageHelper< ImageMask >::ReadImage(maskImageFileName);

            mask = Mask::New();
            mask->SetImage(imageMask);
        }

        std::cout << "Preparing registration..." << std::endl;

        // Extract baseline from diffusion sequence
        Sequence::RegionType extractionRegion = inputSequence->GetLargestPossibleRegion();
        extractionRegion.SetSize(3, 0);

        SequenceExtractor::Pointer extractor = SequenceExtractor::New();
        extractor->SetInput(inputSequence);
        extractor->SetExtractionRegion(extractionRegion);
        extractor->SetDirectionCollapseToSubmatrix();
        extractor->Update();

        Image::Pointer baselineImage = extractor->GetOutput();

        std::cout << "done." << std::endl;


        //
        // Registration
        //

        if(!transformProvided)
        {
            std::cout << "Registering..." << std::endl;

            // Define registration method
            RegistrationMethod::Pointer registration = RegistrationMethod::New();
            registration->SetFixedImage(baselineImage);
            registration->SetMovingImage(referenceImage);

            // Define initial transformation
            initialTransform->GetInverse(initialTransform);
            registration->SetTransform(initialTransform);

            // Define metric used for registration
            MattesMutualInformationMetric::Pointer metric = MattesMutualInformationMetric::New();
            registration->SetMetric(metric);

            // Define image region for registration
            Image::RegionType fixedImageRegion;

            if(!maskImageFileName.empty())
            {
                // Read mask image
                ImageMask::Pointer imageMask = btk::ImageHelper< ImageMask >::ReadImage(maskImageFileName);

                // Create mask object
                Mask::Pointer mask = Mask::New();
                mask->SetImage(imageMask);

                // Set mask
                fixedImageRegion = mask->GetAxisAlignedBoundingBoxRegion();
                metric->SetFixedImageMask(mask);
            }
            else // no mask provided
            {
                fixedImageRegion = baselineImage->GetLargestPossibleRegion();
            }

            // Define the center of rotation
            Image::PointType centerPoint = initialTransform->GetCenter();

            if(initialTransformFileName.empty())
            {
                Image::IndexType centerIndex;
                centerIndex[0] = fixedImageRegion.GetIndex(0) + fixedImageRegion.GetSize(0) / 2.0;
                centerIndex[1] = fixedImageRegion.GetIndex(1) + fixedImageRegion.GetSize(1) / 2.0;
                centerIndex[2] = fixedImageRegion.GetIndex(2) + fixedImageRegion.GetSize(2) / 2.0;

                baselineImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);
                initialTransform->SetCenter(centerPoint);
            }

            // Define initial transformation parameters
            RegistrationParameters initialParameters(initialTransform->GetNumberOfParameters());
            initialParameters = initialTransform->GetParameters();
            registration->SetInitialTransformParameters(initialParameters);

            // Define optimizer used for registration
            RegularStepGradientDescentOptimizer::Pointer optimizer = RegularStepGradientDescentOptimizer::New();
            registration->SetOptimizer(optimizer);

            optimizer->MaximizeOn();
//            optimizer->MinimizeOn();
            optimizer->SetMaximumStepLength(0.2);
            optimizer->SetMinimumStepLength(0.0001);
            optimizer->SetNumberOfIterations(10000);
            optimizer->SetGradientMagnitudeTolerance(0.00001);

            RegularStepGradientDescentOptimizer::ScalesType optimizerScales(initialTransform->GetNumberOfParameters());
            optimizerScales[0] =  1.0;
            optimizerScales[1] =  1.0;
            optimizerScales[2] =  1.0;
            optimizerScales[3] =  1.0;
            optimizerScales[4] =  1.0;
            optimizerScales[5] =  1.0;
            optimizerScales[6] =  1.0;
            optimizerScales[7] =  1.0;
            optimizerScales[8] =  1.0;
            optimizerScales[9] =  1.0;
            optimizerScales[10] =  1.0;
            optimizerScales[11] =  1.0;
            optimizer->SetScales(optimizerScales);

            // Define interpolation used for registration
            LinearInterpolator::Pointer interpolator = LinearInterpolator::New();
            interpolator->SetInputImage(referenceImage);
            registration->SetInterpolator(interpolator);

            // Define the command observer in verbose mode
            if(verboseMode)
            {
                btk::CommandIterationUpdate::Pointer registrationObserver = btk::CommandIterationUpdate::New();
                optimizer->AddObserver(itk::IterationEvent(), registrationObserver);
            }

            // Start registration process
            registration->Update();

            // Get the final transformation parameters
            RegistrationParameters finalParameters = registration->GetLastTransformParameters();

            // Get the final transform
            finalTransform = AffineTransform::New();
            finalTransform->SetCenter(centerPoint);
            finalTransform->SetParameters(finalParameters);
            finalTransform->GetInverse(finalTransform);

            std::cout << "done." << std::endl;
        }
        else // no transformation provided
        {
            finalTransform->Compose(initialTransform, true);
        }


        //
        // Warp sequence to anatomical reference
        //

        std::cout << "Warping..." << std::endl;

        // Define image sequence joiner
        ImageSequenceJoiner::Pointer joiner = ImageSequenceJoiner::New();
        joiner->SetOrigin(0);
        joiner->SetSpacing(1);

        // Each image of diffusion sequence is warped to anatomical reference and resampled
        for(unsigned int i = 0; i < inputSequence->GetLargestPossibleRegion().GetSize(3) && i < 1; i++)
        {
            // Define the current region of interest (1 image is sequence) for extractor
            Sequence::RegionType currentRegion = inputSequence->GetLargestPossibleRegion();
            currentRegion.SetSize(3,0);
            currentRegion.SetIndex(3,i);

            extractor->SetExtractionRegion(currentRegion);
            extractor->Update();

            // Define resampler
            ImageResampler::Pointer resampler = ImageResampler::New();
            resampler->SetTransform(finalTransform);
            resampler->SetInput(extractor->GetOutput());
            resampler->SetUseReferenceImage(true);
            resampler->SetReferenceImage(referenceImage);
            resampler->SetDefaultPixelValue(0);

            // Define interpolation
            BSplineInterpolator::Pointer interpolator = BSplineInterpolator::New();
            interpolator->SetSplineOrder(3);
            resampler->SetInterpolator(interpolator);

            // Resample image and add to output sequence
            resampler->Update();
            joiner->SetInput(i, resampler->GetOutput());
        }

        joiner->Update();

        std::cout << "done." << std::endl;


        //
        // Update gradient table
        //

        std::cout << "Updating gradient table..." << std::endl;

        // Decompose affine transform into rigid and affine part
        // (only rigid part is necessary for gradient table reorientation)
        vnl_matrix< PrecisionType > R(3,3);
        R.set_identity();
        R = finalTransform->GetMatrix().GetVnlMatrix();

        vnl_matrix< PrecisionType > PQ = R;
        vnl_matrix< PrecisionType > NQ = R;
        vnl_matrix< PrecisionType > PQNQDiff;

        for(unsigned int ni = 0; ni < 100; ni++ )
        {
            // Average current Qi with its inverse transpose
            NQ = ( PQ + vnl_inverse_transpose( PQ ) ) / 2.0;
            PQNQDiff = NQ - PQ;

            if( PQNQDiff.frobenius_norm() < 1e-7 )
            {
                break;
            }
            else
            {
                PQ = NQ;
            }
        }

        // Define rigid transformation for gradient table reorientation
        EulerTransform::Pointer gradientTableTransform = EulerTransform::New();
        gradientTableTransform->SetRotationMatrix(NQ);

        // Load input gradient table, rotate and save to new gradient table file
        GradientTable::Pointer inputGradientTable = GradientTable::New();
        inputGradientTable->SetNumberOfGradients(inputSequence->GetLargestPossibleRegion().GetSize(3));
        inputGradientTable->SetImage(inputSequence);
        inputGradientTable->SetTransform(gradientTableTransform);
        inputGradientTable->LoadFromFile(inputSequenceGradientsFileName.c_str());
        inputGradientTable->RotateGradientsInWorldCoordinates();

        std::cout << "done." << std::endl;


        //
        // Write output sequence
        //

        // Save diffusion sequence
        btk::ImageHelper< Sequence >::WriteImage(joiner->GetOutput(), outputSequenceFileName);

        // Save gradient table
        inputGradientTable->SaveToFile(outputSequenceGradientsFileName.c_str());

        // Save b-values
        std::stringstream stream;
        stream << "cp " << inputSequenceBValuesFileName << " " << outputSequenceBValuesFileName;
        std::system(stream.str().c_str());

        // Save transformation (if asked)
        if(!transformFileName.empty())
        {
            AffineTransformFileWriter::Pointer transformWriter = AffineTransformFileWriter::New();
            transformWriter->SetFileName(transformFileName);
            transformWriter->AddTransform(finalTransform);
            transformWriter->Update();
        }
    }
    catch(itk::ExceptionObject &exception)// TODO: manage exceptions
    {
        std::cerr << "ITK error:" << std::endl;
        std::cerr << exception << std::endl;
        exit(EXIT_FAILURE);
    }

  return EXIT_SUCCESS;
}

