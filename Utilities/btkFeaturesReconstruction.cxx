/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 16/12/2013
  Author(s): Julien Pontabry (pontabry@unistra.fr)

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


/**
 * @file btkFeaturesReconstruction.cxx
 * @author Julien Pontabry
 * @date 18/03/2013
 * @ingroup FeatureSelection
 * @brief Reconstruct features dataset given the selected features.
 * For each selected feature, the feature value is made variying in its range
 * ([min,max] of the dataset) and the corresponding dataset is reconstructed
 * consequently. Use this program to show the modes of the dataset.
 */

// TCLAP includes
#include <tclap/CmdLine.h>

// STL includes
#include "fstream"

// ITK includes
#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegionIterator.h"
#include "itkDisplacementFieldTransform.h"

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"
#include "btkFeatureSelectionCostFunction.h"
#include "btkNadarayaWatsonReconstructionErrorFunction.h"

// Definitions
typedef itk::DisplacementFieldTransform< double,3 >::DisplacementFieldType DisplacementField;
typedef itk::ImageRegionIterator< DisplacementField >                      DisplacementFieldIterator;

typedef itk::Image< unsigned char,3 >         ImageMask;
typedef itk::ImageMaskSpatialObject< 3 >      Mask;
typedef itk::ImageRegionIterator< ImageMask > ImageMaskIterator;


/**
 * @brief Main function of the program.
 */
int main(int argc, char *argv[])
{
    try
    {
        //
        // Program's command line
        //

        // Defines the command line parser
        TCLAP::CmdLine cmd("BTK features dataset reconstruction", ' ', "1.0");

        // Defines arguments
        TCLAP::UnlabeledMultiArg< std::string > imagesFileNamesArg("image", "Input displacement field images' filenames", true, "string", cmd);
        TCLAP::MultiArg< double >                 imagesWeightsArg("w", "weight", "Input weights corresponding to images", false, "weights", cmd);
        TCLAP::ValueArg< std::string >            maskFileNamesArg("m", "mask", "Mask filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > reducedMaskImageFileNamesArg("r", "reduced_mask", "Selected features filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >            outputFileNameArg("o", "output", "Output filename prefix (Default: \"reconstruction\")", false, "reconstruction", "string", cmd);

        // Optional arguments
        TCLAP::ValueArg< unsigned int > numberOfSamplesPerFeatureArg("s","number_of_samples", "Number of samples per feature (for reconstruction, default; 3)",false, 3, "positive integer", cmd);
        TCLAP::ValueArg< float >                  kernelBandwidthArg("b", "kernel_bandwidth", "Bandwidth of the kernel use for cost function", false, 1.0, "positive real", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::vector< std::string > imagesFileNames = imagesFileNamesArg.getValue();
        std::vector< double >        imagesWeights = imagesWeightsArg.getValue();
        std::string                   maskFileName = maskFileNamesArg.getValue();
        std::string       reducedMaskImageFileName = reducedMaskImageFileNamesArg.getValue();
        std::string                 outputFileName = outputFileNameArg.getValue();

        unsigned int numberOfSamplesPerFeature = numberOfSamplesPerFeatureArg.getValue();
        float                  kernelBandwidth = kernelBandwidthArg.getValue();


        //
        // Read input file
        //

        // Test if there are as many weights as images
        if(imagesWeights.empty())
        {
            imagesWeights = std::vector< double >(imagesFileNames.size(), 1.0);
        }
        else if(imagesFileNames.size() != imagesWeights.size())
        {
            btkException("The number of weights is not the same as the number of images !");
        }

        // Read images
        std::vector< DisplacementField::Pointer > &images = btk::ImageHelper< DisplacementField >::ReadImage(imagesFileNames);

        // Test if images are in the same space
        if(!btk::ImageHelper< DisplacementField >::IsInSamePhysicalSpace(images))
        {
            btkException("Input images are not in the same physical space !");
        }

        // Read mask image
        ImageMask::Pointer maskImage = btk::ImageHelper< ImageMask >::ReadImage(maskFileName);

        // Test if images and mask are in the same space
        if(!btk::ImageHelper< DisplacementField,ImageMask >::IsInSamePhysicalSpace(images[0],maskImage, 1e-4))
        {
            btkException("Input images and mask image are not in the same physical space !");
        }

        // Read selected features mask
        ImageMask::Pointer reducedMaskImage = btk::ImageHelper< ImageMask >::ReadImage(reducedMaskImageFileName);

        // Test if images and mask are in the same space
        if(!btk::ImageHelper< DisplacementField,ImageMask >::IsInSamePhysicalSpace(images[0],reducedMaskImage, 1e-4))
        {
            btkException("Input images and reduced mask image are not in the same physical space !");
        }


        //
        // Preprocessing
        //

        std::cout << "Preprocessing..." << std::endl;

        // Compute the object mask
        Mask::Pointer mask = Mask::New();
        mask->SetImage(maskImage);
        mask->Update();

        // Get the masked region
        Mask::RegionType maskedRegion = mask->GetAxisAlignedBoundingBoxRegion();

        // Compute the number of parameters
        unsigned int numberOfParameters = 0, numberOfReducedParameters = 0;
        ImageMaskIterator maskIt(maskImage, maskedRegion);
        ImageMaskIterator reducedMaskIt(reducedMaskImage, maskedRegion);

        for(maskIt.GoToBegin(), reducedMaskIt.GoToBegin(); !maskIt.IsAtEnd() && !reducedMaskIt.IsAtEnd(); ++maskIt, ++reducedMaskIt)
        {
            if(maskIt.Get() > 0)
            {
                numberOfParameters++;

                if(reducedMaskIt.Get() > 0)
                {
                    numberOfReducedParameters++;
                }
            }
        } // for each voxels

        btkCoutMacro("\tThere are " << numberOfParameters << " initial parameters and " << numberOfReducedParameters << " reduced parameters.");


        // Build input and reduced parameters
        vnl_matrix< double > Y(numberOfParameters*3, images.size());
        vnl_matrix< double > X(numberOfReducedParameters*3, images.size());
        vnl_vector< double > w(images.size());
        unsigned int j = 0;

        for(std::vector< DisplacementField::Pointer >::iterator image = images.begin(); image != images.end(); image++, j++)
        {
            // Fill image weights
            w(j) = imagesWeights[j];

            // Fill parameters
            unsigned int i = 0, i2 = 0;

            DisplacementFieldIterator imageIt(*image, maskedRegion);

            for(maskIt.GoToBegin(), reducedMaskIt.GoToBegin(), imageIt.GoToBegin(); !maskIt.IsAtEnd() && !reducedMaskIt.IsAtEnd() && !imageIt.IsAtEnd(); ++maskIt, ++reducedMaskIt, ++imageIt)
            {
                if(maskIt.Get() > 0)
                {
                    Y(i,j) = imageIt.Get()[0];
                    Y(i+numberOfParameters,j) = imageIt.Get()[1];
                    Y(i+2*numberOfParameters,j) = imageIt.Get()[2];
                    i++;

                    if(reducedMaskIt.Get() > 0)
                    {
                        X(i2,j) = imageIt.Get()[0];
                        X(i2+numberOfReducedParameters,j) = imageIt.Get()[1];
                        X(i2+2*numberOfReducedParameters,j) = imageIt.Get()[2];
                        i2++;
                    }
                }
            } // for each voxels

            // Clean memory
            *image = NULL;
        } // for each input subject

        // Clean memory
        mask = NULL;
        reducedMaskImage = NULL;

        std::cout << "done." << std::endl;


        //
        // Processing
        //

        btkCoutMacro("Computing selected features range...");

        // Compute range of selected features
        vnl_matrix< double > range(X.rows(),numberOfSamplesPerFeature);
        range.fill(0.0);

        for(unsigned int i = 0; i < range.rows(); i++)
        {
            // min range
            range(i,0) = X.get_row(i).min_value();

            // max range
            range(i,range.columns()-1) = X.get_row(i).max_value();

            // increment
            double increment = (range(i,range.columns()-1)-range(i,0)) / static_cast< double >(numberOfSamplesPerFeature-1);

            // Fill the other columns of the matrix
            for(unsigned int j = 1; j < range.cols()-1; j++)
            {
                range(i,j) += range(i,j-1) + increment;
            } // for each column
        } // for each row

        btkCoutMacro("done.");

        btkCoutMacro("Computing mean vector...");

        // Compute the mean vector of reduced parameters
        vnl_vector< double > meanParameters(X.rows());

        for(int i = 0; i < meanParameters.size(); i++)
        {
            meanParameters(i) = X.get_row(i).mean();
        } // for each row

        std::cerr << std::endl;

        btkCoutMacro("done.");

        btkCoutMacro("Computing the varying parameters matrix...");

        // Compute the varying parameters matrix
        vnl_matrix< double > parameters(range.rows(),range.rows()*numberOfSamplesPerFeature);

        for(unsigned int j = 0; j < parameters.columns(); j++)
        {
            parameters.set_column(j, meanParameters);
        }

        for(unsigned int i = 0; i < range.rows(); i++)
        {
            for(unsigned int j = 0; j < range.columns(); j++)
            {
                parameters(i,numberOfSamplesPerFeature*i + j) = range(i,j);
            }
        }

        btkCoutMacro("done.");

        // Compute bandwidth matrix diagonal matrix
        vnl_diag_matrix< double > H(numberOfReducedParameters*3);
        H.fill_diagonal(kernelBandwidth);

        // Define cost function
        btk::FeatureSelectionCostFunction::Pointer costFunction = NULL;

        // We use a kernel based cost function
        btk::NadarayaWatsonReconstructionErrorFunction::Pointer kernelCostFunction = btk::NadarayaWatsonReconstructionErrorFunction::New();
        kernelCostFunction->SetBandwidthMatrix(&H);
        kernelCostFunction->SetImagesWeightVector(&w);
        costFunction = kernelCostFunction;

        // Initialize cost function
        costFunction->SetInputParameters(&Y);
        costFunction->SetCurrentParameters(&X);
        costFunction->Initialize();

        btkCoutMacro("Computing reconstruction...");

        // Compute reconstruction
        vnl_matrix< double > reconstruction = costFunction->GetReconstruction(parameters);

        btkCoutMacro("done.");

        // Save the reconstructed images
        for(unsigned int j = 0; j < reconstruction.columns(); j++)
        {
            DisplacementField::Pointer reconstructedImage = btk::ImageHelper< ImageMask,DisplacementField >::CreateNewImageFromPhysicalSpaceOf(maskImage);

            unsigned int i = 0;

            DisplacementFieldIterator imageIt(reconstructedImage, maskedRegion);

            for(maskIt.GoToBegin(), imageIt.GoToBegin(); !maskIt.IsAtEnd() && !imageIt.IsAtEnd(); ++maskIt, ++imageIt)
            {
                DisplacementField::PixelType pixel;

                if(maskIt.Get() > 0)
                {
                    pixel[0] = reconstruction(i,j);
                    pixel[1] = reconstruction(i+numberOfParameters,j);
                    pixel[2] = reconstruction(i+2*numberOfParameters,j);
                    imageIt.Set(pixel);
                    i++;
                }
                else // maskIt.Get <= 0
                {
                    pixel[0] = pixel[1] = pixel[2] = 0.0;
                    imageIt.Set(pixel);
                }
            }

            // Write file
            std::stringstream s;
            s << outputFileName << "_feature-" << j / numberOfSamplesPerFeature << "_increment-" << j % numberOfSamplesPerFeature << ".nii.gz";

            btk::ImageHelper< DisplacementField >::WriteImage(reconstructedImage, s.str());
        }
    }
    catch(TCLAP::ArgException &exception)
    {
        btkCoutMacro("Command line error: " << exception.error() << " for argument " << exception.argId());
        exit(EXIT_FAILURE);
    }
    catch(itk::ExceptionObject &e)
    {
        btkCoutMacro("Error (" << e.GetLocation() << "): " << e.GetDescription() << std::endl);
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}

