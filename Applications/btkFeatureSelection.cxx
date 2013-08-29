/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 19/12/2012
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
 * @file btkFeatureSelection.cxx
 * @author Julien Pontabry
 * @date 19/12/2012
 * @ingroup FeatureSelection
 * @brief Space reduction program by feature selection.
 *
 * Reduce the space of a set of displacement fields by a greedy feature selection algorithm. The features are the vectors.
 * As input, there should be at least the displacement fields (and eventually weights for these fields) and a mask.
 * The program produces two outputs: a mask pointing the location of the selected features and an explained variance map,
 * expressed as percent of total variance.
 *
 * The complexity of the complete algorithm (i.e. the sequential floating search algorithm and the cost function) is in O(n p² K²),
 * where n, p and K are respectively the number of selected features, the number of initial features and the number of samples.
 *
 * This program is intended to shape change analysis on brain fetuses by feature selection on displacement fields.
 * You can find more information in the following publication:
 *
 * Pontabry et al. (2013) Sélection de caractéristiques pour l'étude de la maturation cérébrale, in ORASIS: journées francophones
 * des jeunes chercheurs en vision par ordinateurs, Cluny (FRANCE), Oral presentation.
 *
 * @todo Use multi-resolution to improve the speedness of the program.
 * @todo Extend the program to an arbitrary number of components.
 */

// TCLAP includes
#include <tclap/CmdLine.h>

// STL includes
#include "cstdlib"
#include "iostream"
#include "string"
#include "vector"
#include "cfloat"

// ITK includes
#include "itkImage.h"
#include "itkArray.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_matrix.h"
#include "itkDisplacementFieldTransform.h"

// Local includes
#include "btkImageHelper.h"
#include "btkFeatureSelectionAlgorithm.h"
#include "btkGreedyFeatureSelectionAlgorithm.h"
#include "btkFeatureSelectionCostFunction.h"
#include "btkNadarayaWatsonReconstructionErrorFunction.h"
#include "btkFeatureSelectionCommandIterationUpdate.h"


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
        TCLAP::CmdLine cmd("BTK space reduction program for a set of displacement fields (space reduction by feature selection)", ' ', "1.0");

        // Defines arguments
        TCLAP::MultiArg< std::string > imagesFileNamesArg("i","image", "Input displacement field images' filenames", true, "string", cmd);
        TCLAP::MultiArg< double >        imagesWeightsArg("w", "weight", "Input weights corresponding to images", false, "weights", cmd);
        TCLAP::ValueArg< std::string >   maskFileNamesArg("m", "mask", "Mask filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >    outputPrefixArg("o", "output", "Output reduced filenames prefix", false, "reduced", "string", cmd);

        // Options
        TCLAP::ValueArg< unsigned int > maxNumberOfParametersArg("n", "max_number_of_parameters", "Maximal number of parameters", false, 10, "unsigned int", cmd);
        TCLAP::ValueArg< float >              kernelBandwidthArg("b", "kernel_bandwidth", "Bandwidth of the kernel use for cost function", false, 1.0, "positive real", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::vector< std::string > imagesFileNames = imagesFileNamesArg.getValue();
        std::vector< double >        imagesWeights = imagesWeightsArg.getValue();
        std::string                   maskFileName = maskFileNamesArg.getValue();
        std::string                   outputPrefix = outputPrefixArg.getValue();

        unsigned int maxNumberOfParameters = maxNumberOfParametersArg.getValue();
        float              kernelBandwidth = kernelBandwidthArg.getValue();


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
        if(!btk::ImageHelper< DisplacementField,ImageMask >::IsInSamePhysicalSpace(images[0],maskImage, 10e-5))
        {
            btkException("Input images and the mask image are not in the same physical space !");
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
        unsigned int numberOfParameters = 0;
        ImageMaskIterator maskIt(maskImage, maskedRegion);

        for(maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt)
        {
            if(maskIt.Get() > 0)
            {
                numberOfParameters++;
            }
        } // for each voxels

        btkCoutMacro("\tThere are " << numberOfParameters << " parameters and " << images.size() << " learning samples.");


        // Build input parameters
        vnl_matrix< double > Y(numberOfParameters*3, images.size());
        vnl_vector< double > w(images.size());
        unsigned int j = 0;

        for(std::vector< DisplacementField::Pointer >::iterator image = images.begin(); image != images.end(); image++, j++)
        {
            // Fill image weights
            w(j) = imagesWeights[j];

            // Fill parameters
            unsigned int i = 0;

            DisplacementFieldIterator imageIt(*image, maskedRegion);

            for(maskIt.GoToBegin(), imageIt.GoToBegin(); !maskIt.IsAtEnd() && !imageIt.IsAtEnd(); ++maskIt, ++imageIt)
            {
                if(maskIt.Get() > 0)
                {
                    Y(i,j) = imageIt.Get()[0];
                    Y(i+numberOfParameters,j) = imageIt.Get()[1];
                    Y(i+2*numberOfParameters,j) = imageIt.Get()[2];
                    i++;
                }
            } // for each voxels

            // Clean memory
            *image = NULL;
        } // for each input subject

        // Clean memory
        mask = NULL;

        std::cout << "done." << std::endl;


        //
        // Processing
        //

        std::cout << "Processing..." << std::endl;

        vnl_diag_matrix< double > H(numberOfParameters*3);
        H.fill_diagonal(kernelBandwidth);

        // Define cost function
        btk::FeatureSelectionCostFunction::Pointer costFunction = NULL;

        // We use a kernel based cost function
        btk::NadarayaWatsonReconstructionErrorFunction::Pointer kernelCostFunction = btk::NadarayaWatsonReconstructionErrorFunction::New();
        kernelCostFunction->SetBandwidthMatrix(&H);
        kernelCostFunction->SetImagesWeightVector(&w);
        costFunction = kernelCostFunction;

        // Define optimizer for feature selection
        btk::FeatureSelectionAlgorithm::Pointer spaceReduction = NULL;

        // Use greedy algorithm
        btk::GreedyFeatureSelectionAlgorithm::Pointer greedySpaceReduction = btk::GreedyFeatureSelectionAlgorithm::New();
        spaceReduction = greedySpaceReduction;

        btkCoutMacro("\tWill use a greedy algorithm (sequential floating forward-backward search).");

        // Set general paramters for optimizer
        spaceReduction->SetInputParameters(&Y);
        spaceReduction->SetCostFunction(costFunction);
        spaceReduction->SetMaxNumberOfParameters(maxNumberOfParameters);

        // Set a command-line observer
        btk::FeatureSelectionCommandIterationUpdate::Pointer observer = btk::FeatureSelectionCommandIterationUpdate::New();
        spaceReduction->AddObserver(itk::IterationEvent(), observer);

        // Start reduction algorithm
        spaceReduction->Update();

        // Get weigths vectors
        vnl_vector< short > &weights = *spaceReduction->GetWeightsVector();

        // Get energy vector
        vnl_vector< double > &energy = *spaceReduction->GetEnergyVector();

        std::cout << "done." << std::endl;


        //
        // Export output paramters
        //

        std::cout << "Creating and saving outputs from reduced set of parameters..." << std::endl;

        // Create output images
        ImageMask::Pointer                maskOutput = btk::ImageHelper< ImageMask >::CreateNewImageFromPhysicalSpaceOf(maskImage);
        itk::Image< double,3 >::Pointer energyOutput = btk::ImageHelper< ImageMask,itk::Image< double,3 > >::CreateNewImageFromPhysicalSpaceOf(maskImage);

        ImageMaskIterator maskOutputIt(maskOutput, maskedRegion);
        itk::ImageRegionIterator< itk::Image< double,3 > > energyOutputIt(energyOutput, maskedRegion);

        // Fill it with the reduced parameters
        unsigned int i = 0;

        for(maskIt.GoToBegin(), maskOutputIt.GoToBegin(), energyOutputIt.GoToBegin(); !maskIt.IsAtEnd() && !maskOutputIt.IsAtEnd() && !energyOutputIt.IsAtEnd(); ++maskIt, ++maskOutputIt, ++energyOutputIt)
        {
            if(maskIt.Get() > 0)
            {
                if(weights(i) > 0)
                {
                    energyOutputIt.Set(energy(i));
                    maskOutputIt.Set(static_cast< unsigned char >(1));
                }

                i++;
            }
        } // for each voxels


        // Write it
        btk::ImageHelper< ImageMask >::WriteImage(maskOutput, outputPrefix+"_selected-features.nii.gz");
        btk::ImageHelper< itk::Image< double,3 > >::WriteImage(energyOutput, outputPrefix+"_explained-variance.nii.gz");

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
