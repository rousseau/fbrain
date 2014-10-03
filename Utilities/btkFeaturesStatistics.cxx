/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 10/02/2014
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
 * @file btkFeaturesStatistics.cxx
 * @author Julien Pontabry
 * @date 10/02/2014
 * @ingroup FeatureSelection
 * @brief Compute statistics on features
 * For each feature, statistics are computed (mean, std for instance).
 * The dynamic of reduced feature is computed too (if any).
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

typedef itk::Image< double,3 >                  ScalarImage;
typedef itk::ImageRegionIterator< ScalarImage > ScalarImageIterator;

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
        TCLAP::ValueArg< std::string >            maskFileNamesArg("m", "mask", "Mask filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >            outputFileNameArg("o", "output", "Output filename prefix (Default: \"statistics\")", false, "statistics", "string", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::vector< std::string > imagesFileNames = imagesFileNamesArg.getValue();
        std::string                   maskFileName = maskFileNamesArg.getValue();
        std::string                 outputFileName = outputFileNameArg.getValue();


        //
        // Read input file
        //

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

        btkCoutMacro("\tThere are " << numberOfParameters << " parameters.");


        // Build input and reduced parameters
        vnl_matrix< double > Y(numberOfParameters*3, images.size());
        unsigned int j = 0;

        for(std::vector< DisplacementField::Pointer >::iterator image = images.begin(); image != images.end(); image++, j++)
        {
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

        btkCoutMacro("Computing statistics...");

        // Compute the mean vector of parameters
        vnl_vector< double > meanParameters(Y.rows());

        for(int i = 0; i < meanParameters.size(); i++)
        {
            meanParameters(i) = Y.get_row(i).mean();
        } // for each row


        // Compute the std deviation of parameters
        vnl_vector< double > stdParameters(Y.rows());

        for(int i = 0; i < stdParameters.size(); i++)
        {
            double sum = 0.0;

            for(int j = 0; j < Y.columns(); j++)
            {
                sum += (Y(i,j) - meanParameters(i)) * (Y(i,j) - meanParameters(i));
            }

            sum /= Y.columns()-1.0;

            stdParameters(i) = std::sqrt(sum);
        }

        // Compute the mean of magnitudes for each sample
        vnl_vector< double > sampleMeansParameters(Y.cols());

        for(int j = 0; j < sampleMeansParameters.size(); j++)
        {
            double sum = 0.0;

            for(int p = 0; p < numberOfParameters; p++)
            {
                // Make the sum of magnitudes
                sum += std::sqrt(
                    Y(p,j)*Y(p,j) +
                    Y(p+numberOfParameters,j)*Y(p+numberOfParameters,j) +
                    Y(p+2*numberOfParameters,j)*Y(p+2*numberOfParameters,j)
                );
            }

            // Make the average
            sampleMeansParameters(j) = sum / static_cast< double >(numberOfParameters);
        } // for each column

        btkCoutMacro("done.");

        // Save the statistics images
        DisplacementField::Pointer meanImage    = btk::ImageHelper< ImageMask,DisplacementField >::CreateNewImageFromPhysicalSpaceOf(maskImage);
        DisplacementField::Pointer stdImage     = btk::ImageHelper< ImageMask,DisplacementField >::CreateNewImageFromPhysicalSpaceOf(maskImage);
        ScalarImage::Pointer meanMagnitudeImage = btk::ImageHelper< ImageMask,ScalarImage >::CreateNewImageFromPhysicalSpaceOf(maskImage);

        unsigned int i = 0;

        DisplacementFieldIterator meanImageIt(meanImage, maskedRegion);
        DisplacementFieldIterator stdImageIt(stdImage, maskedRegion);
        ScalarImageIterator meanMagnitudeImageIt(meanMagnitudeImage, maskedRegion);

        for(maskIt.GoToBegin(), meanImageIt.GoToBegin(), stdImageIt.GoToBegin(), meanMagnitudeImageIt.GoToBegin(); !maskIt.IsAtEnd() && !meanImageIt.IsAtEnd() && !stdImageIt.IsAtEnd() && !meanMagnitudeImageIt.IsAtEnd(); ++maskIt, ++meanImageIt, ++stdImageIt, ++meanMagnitudeImageIt)
        {
            DisplacementField::PixelType pixelMean, pixelStd;
            ScalarImage::PixelType pixelMeanMagnitude;

            if(maskIt.Get() > 0)
            {
                pixelMean[0] = meanParameters(i);
                pixelMean[1] = meanParameters(i+numberOfParameters);
                pixelMean[2] = meanParameters(i+2*numberOfParameters);
                meanImageIt.Set(pixelMean);

                pixelStd[0] = stdParameters(i);
                pixelStd[1] = stdParameters(i+numberOfParameters);
                pixelStd[2] = stdParameters(i+2*numberOfParameters);
                stdImageIt.Set(pixelStd);

                pixelMeanMagnitude = std::sqrt(
                            meanParameters(i)*meanParameters(i) +
                            meanParameters(i+numberOfParameters)*meanParameters(i+numberOfParameters) +
                            meanParameters(i+2*numberOfParameters)*meanParameters(i+2*numberOfParameters)
                );
                meanMagnitudeImageIt.Set(pixelMeanMagnitude);
                i++;
            }
            else // maskIt.Get <= 0
            {
                pixelMean[0] = pixelMean[1] = pixelMean[2] = 0.0;
                meanImageIt.Set(pixelMean);

                pixelStd[0] = pixelStd[1] = pixelStd[2] = 0.0;
                stdImageIt.Set(pixelStd);

                pixelMeanMagnitude = 0.0;
                meanMagnitudeImageIt.Set(pixelMeanMagnitude);
            }
        }

        // Write files
        std::stringstream sMean;
        sMean << outputFileName << "_mean.nii.gz";

        btk::ImageHelper< DisplacementField >::WriteImage(meanImage, sMean.str());

        std::stringstream sStd;
        sStd << outputFileName << "_std.nii.gz";

        btk::ImageHelper< DisplacementField >::WriteImage(stdImage, sStd.str());

        std::stringstream sMeanMagnitude;
        sMeanMagnitude << outputFileName << "_mean-magnitude.nii.gz";

        btk::ImageHelper< ScalarImage >::WriteImage(meanMagnitudeImage, sMeanMagnitude.str());


        // Save sample statistics
        std::stringstream sSampleMeanMagnitude;
        sSampleMeanMagnitude << outputFileName << "_sample-mean-magnitude.csv";

        std::ofstream file(sSampleMeanMagnitude.str().c_str());

        if(file.is_open())
        {
            for(int j = 0; j < Y.cols(); j++)
            {
                file << sampleMeansParameters(j) << std::endl;
            } // for each column

            // Close file
            file.close();
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


