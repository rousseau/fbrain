/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 18/03/2013
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
 * @file btkComputeFeatureSelectionResiduals.cxx
 * @author Julien Pontabry
 * @date 18/03/2013
 * @ingroup FeatureSelection
 * @brief Compute the residuals of a feature selection.
 */

// TCLAP includes
#include <tclap/CmdLine.h>

// STL includes
#include "cstdlib"
#include "iostream"
#include "string"
#include "vector"
#include "cfloat"
#include "sstream"

// ITK includes
#include "itkImage.h"
#include "itkArray.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_matrix.h"
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

typedef itk::Image< double,3 >                  EnergyImage;
typedef itk::ImageRegionIterator< EnergyImage > EnergyImageIterator;

typedef itk::Image< double,3 >                    ResidualImage;
typedef itk::ImageRegionIterator< ResidualImage > ResidualImageIterator;

/**
 * @brief Sort the selected features.
 * Features are sorted using their explained variance.
 * @param explainedVariance Explained variance of features.
 * @return Vector of indices corresponding to sorted features.
 */
vnl_vector< int > SortFeatures(vnl_vector< double > &explainedVariance)
{
    // Set an indices vector for selected features
    vnl_vector< int > indices(explainedVariance.size());

    for(int i = 0; i < indices.size(); i++)
    {
        indices(i) = i;
    }

    // Bubble sort
    bool convergence;

    do
    {
        convergence = true;

        for(int i = 0; i < explainedVariance.size()-1; i++)
        {
            if(explainedVariance(i) < explainedVariance(i+1))
            {
                int     tmpI = indices(i);
                indices(i)   = indices(i+1);
                indices(i+1) = tmpI;

                double            tmpE = explainedVariance(i);
                explainedVariance(i)   = explainedVariance(i+1);
                explainedVariance(i+1) = tmpE;

                convergence = false;
            }
        }
    } while(!convergence);

    // Return indices corresponding to sorted features
    return indices;
}


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
        TCLAP::CmdLine cmd("BTK residuals computation program of space reduction", ' ', "1.0");

        // Defines arguments
        TCLAP::UnlabeledMultiArg< std::string > imagesFileNamesArg("image", "Input displacement field images' filenames", true, "string", cmd);
        TCLAP::ValueArg< std::string >           energyFileNameArg("e", "energy", "Energy of selected features", true, "", "string", cmd);
        TCLAP::MultiArg< double >                 imagesWeightsArg("w", "weight", "Input weights corresponding to images", false, "weights", cmd);
        TCLAP::ValueArg< std::string >            maskFileNamesArg("m", "mask", "Mask filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >     reducedmaskFileNamesArg("r", "reduced_mask", "Reduced mask filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >           outputFileNameArg("o", "output", "Output residuals filename", false, "reduction_residuals.nii.gz", "string", cmd);

        // Options
        TCLAP::ValueArg< float > kernelBandwidthArg("b", "kernel_bandwidth", "Bandwidth of the kernel use for cost function", false, 1.0, "positive real", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::vector< std::string > imagesFileNames = imagesFileNamesArg.getValue();
        std::string                 energyFileName = energyFileNameArg.getValue();
        std::vector< double >        imagesWeights = imagesWeightsArg.getValue();
        std::string                   maskFileName = maskFileNamesArg.getValue();
        std::string           reducedmaskFileNames = reducedmaskFileNamesArg.getValue();
        std::string                 outputFileName = outputFileNameArg.getValue();

        float kernelBandwidth = kernelBandwidthArg.getValue();


        //
        // Read input files
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

        // Read reduced mask image
        ImageMask::Pointer reducedMaskImage = btk::ImageHelper< ImageMask >::ReadImage(reducedmaskFileNames);

        // Test if images and mask are in the same space
        if(!btk::ImageHelper< DisplacementField,ImageMask >::IsInSamePhysicalSpace(images[0],reducedMaskImage, 1e-4))
        {
            btkException("Input images and reduced mask image are not in the same physical space !");
        }

        // Read energy image
        EnergyImage::Pointer energyImage = btk::ImageHelper< EnergyImage >::ReadImage(energyFileName);

        // Test if images and energy image are in the same space
        if(!btk::ImageHelper< DisplacementField,EnergyImage >::IsInSamePhysicalSpace(images[0], energyImage, 1e-4))
        {
            btkException("Input images and energy image are not in the same physical space !");
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


        // Build input parameters
        vnl_matrix< double > Y(numberOfParameters*3, images.size());
        vnl_matrix< double > X(numberOfReducedParameters*3, images.size());
        vnl_vector< double > w(images.size());
        unsigned int j = 0;

        EnergyImageIterator energyIt(energyImage, maskedRegion);
        vnl_vector< double > energy(numberOfReducedParameters);

        for(std::vector< DisplacementField::Pointer >::iterator image = images.begin(); image != images.end(); image++, j++)
        {
            // Fill image weights
            w(j) = imagesWeights[j];

            // Fill parameters
            unsigned int i = 0, i2 = 0;
            unsigned int k = 0;

            DisplacementFieldIterator imageIt(*image, maskedRegion);

            for(maskIt.GoToBegin(), reducedMaskIt.GoToBegin(), energyIt.GoToBegin(), imageIt.GoToBegin(); !maskIt.IsAtEnd() && !reducedMaskIt.IsAtEnd() && !energyIt.IsAtEnd() && !imageIt.IsAtEnd(); ++maskIt, ++reducedMaskIt, ++energyIt, ++imageIt)
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

                        energy(k++) = energyIt.Get();
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

        std::cout << "Processing..." << std::endl;

        // Sort selected features by their energy (or explained variance)
        vnl_vector< int > sortedFeatures = SortFeatures(energy);

        // Compute bandwidth matrix diagonal matrix
        vnl_diag_matrix< double > H(numberOfReducedParameters*3);
        H.fill_diagonal(kernelBandwidth);

        // Setup activation diagonal matrix
        vnl_diag_matrix< double > W(numberOfReducedParameters*3);
        W.fill(0.0);

        int actualNumberOfReducedParameters = 0;

        // No elements in selected features matrix at first
        vnl_matrix< double > tmpX = W * X;

        // compute and save residuals and add elements by ascending order of energy
        do
        {
            // Define cost function
            btk::FeatureSelectionCostFunction::Pointer costFunction = NULL;

            // We use a kernel based cost function
            btk::NadarayaWatsonReconstructionErrorFunction::Pointer kernelCostFunction = btk::NadarayaWatsonReconstructionErrorFunction::New();
            kernelCostFunction->SetBandwidthMatrix(&H);
            kernelCostFunction->SetImagesWeightVector(&w);
            costFunction = kernelCostFunction;

            // Initialize cost function
            costFunction->SetInputParameters(&Y);
            costFunction->SetCurrentParameters(&tmpX);
            costFunction->Initialize();
            std::cerr << costFunction->OptimizeParameters() << std::endl;

            // Compute residuals
            vnl_matrix< double > residuals = costFunction->GetResiduals();

            // Save the residuals images
            for(unsigned int j = 0; j < residuals.columns(); j++)
            {
                ResidualImage::Pointer residualsImage = btk::ImageHelper< ImageMask,ResidualImage >::CreateNewImageFromPhysicalSpaceOf(maskImage);

                unsigned int i = 0;

                ResidualImageIterator imageIt(residualsImage, maskedRegion);

                for(maskIt.GoToBegin(), imageIt.GoToBegin(); !maskIt.IsAtEnd() && !imageIt.IsAtEnd(); ++maskIt, ++imageIt)
                {
                    if(maskIt.Get() > 0)
                    {
                        imageIt.Set(residuals(i,j)*residuals(i,j) + residuals(i+numberOfParameters,j)*residuals(i+numberOfParameters,j) + residuals(i+2*numberOfParameters,j)*residuals(i+2*numberOfParameters,j));
                        i++;
                    }
                    else // maskIt.Get <= 0
                    {
                        imageIt.Set(0.0);
                    }
                }

                // Write file
                std::stringstream s;
                s << outputFileName << imagesFileNames[j] << "_" << actualNumberOfReducedParameters << ".nii.gz";

                btk::ImageHelper< ResidualImage >::WriteImage(residualsImage, s.str());
            }

            if(actualNumberOfReducedParameters < numberOfReducedParameters)
            {
                // Add next parameter to the set
                W(sortedFeatures(actualNumberOfReducedParameters))                               = 1.0;
                W(sortedFeatures(actualNumberOfReducedParameters) + numberOfReducedParameters)   = 1.0;
                W(sortedFeatures(actualNumberOfReducedParameters) + 2*numberOfReducedParameters) = 1.0;

                // Update reduced set of parameters
                tmpX = W * X;
            }

            actualNumberOfReducedParameters++;
        }
        while(actualNumberOfReducedParameters <= numberOfReducedParameters);

        std::cout << "done." << std::endl;
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
