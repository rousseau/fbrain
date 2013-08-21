/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 27/05/2013
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
 * @file btkComputeCoefficientOfDeterminationForAtlas.cxx
 * @author Julien Pontabry
 * @date 27/05/2013
 * @brief Compute the coefficient of determination for atlas built with btk scripts.
 */

// STL includes
#include "string"
#include "vector"
#include "algorithm"

// ITK includes
#include "itkImage.h"
#include "itkDisplacementFieldTransform.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

// TCLAP includes
#include <tclap/CmdLine.h>

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"


// Type definitions
typedef double ScalarType;
const unsigned int Dimension = 3;

typedef itk::Image< ScalarType,Dimension >                                             IntensitiesImage;
typedef itk::DisplacementFieldTransform< ScalarType,Dimension >::DisplacementFieldType DisplacementField;
typedef itk::Image< short,Dimension >                                                  MaskImage;


/**
 * @brief Templated function for processing.
 * This function is templated over the input image type.
 * @param inputFileNames Vector of input image filenames.
 * @param modeFileNames Vector of corresponding model image filenames.
 * @param outputFileName Output image filename.
 * @param maskFileName Mask image filename.
 */
template< class T >
void ComputeDeterminationCoefficient(std::vector< std::string > &inputFileNames, std::vector< std::string > &modeFileNames, std::string &outputFileName, std::string &maskFileName);


/**
 * @brief Main function of the program.
 */
int main(int argc, char *argv[])
{
    try
    {
        //
        // Command line processing
        //

        // Define command line
        TCLAP::CmdLine cmd("Compute the coefficient of determination for atlas built with btk scripts.", ' ', "1.0", true);

        // Define arguments
        TCLAP::MultiArg< std::string > inputFileNamesArg("i", "input_image", "Input images' filenames (sample)", true, "string", cmd);
        TCLAP::MultiArg< std::string > modelFileNamesArg("m", "model_image", "Model images' filenames (regression corresponding to sample)", true, "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output_image", "Output image filename", false, "determination.nii.gz", "string", cmd);

        // Option
        TCLAP::ValueArg< std::string > maskFileNameArg("", "mask", "Mask image filename", false, "", "string", cmd);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each argument
        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();
        std::vector< std::string > modelFileNames = modelFileNamesArg.getValue();
        std::string                outputFileName = outputFileNameArg.getValue();

        std::string maskFileName = maskFileNameArg.getValue();


        //
        // Processing depending on image type
        //

        // Checking input data
        if(inputFileNames.size() != modelFileNames.size())
        {
            btkException("Error: There is not the same number of input and model images ! There should be the same number exactly.");
        }

        // Opening header
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFileNames[0].c_str(), itk::ImageIOFactory::ReadMode);
        imageIO->SetFileName(inputFileNames[0]);
        imageIO->ReadImageInformation();

        // The dimension of the images should be 3
        if(imageIO->GetNumberOfDimensions() != 3)
        {
            btkException("Error: Unsupported image dimension ! Only dimension 3 is supported.");
        }

        // Selecting appropriate structures
        switch(imageIO->GetPixelType())
        {
            case itk::ImageIOBase::VECTOR:
                ComputeDeterminationCoefficient< DisplacementField >(inputFileNames, modelFileNames, outputFileName, maskFileName);
                break;

            case itk::ImageIOBase::SCALAR:
                ComputeDeterminationCoefficient< IntensitiesImage >(inputFileNames, modelFileNames, outputFileName, maskFileName);
                break;

            default:
                btkException("Error: Unsupported pixel type ! Only vector and scalar type are supported.");
        }
    }
    catch (TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch(std::string &err)
    {
        std::cerr << err << std::endl;
    }

    return EXIT_SUCCESS;
}


template< class T >
void ComputeDeterminationCoefficient(std::vector< std::string > &inputFileNames, std::vector< std::string > &modeFileNames, std::string &outputFileName, std::string &maskFileName)
{
    // Define types
    typedef itk::ImageRegionIterator< T >                Iterator;
    typedef itk::ImageRegionIterator< IntensitiesImage > IntensitiesIterator;
    typedef itk::ImageRegionIterator< MaskImage >        MaskIterator;


    // Read images
    std::vector< typename T::Pointer > inputImages = btk::ImageHelper< T >::ReadImage(inputFileNames);
    std::vector< typename T::Pointer > modelImages = btk::ImageHelper< T >::ReadImage(modeFileNames);
    MaskImage::Pointer                        mask = btk::ImageHelper< MaskImage >::ReadImage(maskFileName);


    // Variables for progress indications
    double percent     = 0.0;
    double percentStep = 100.0 / static_cast< double >(inputImages.size());


    // Variable for output filename
    std::string outputPrefix = btk::FileHelper::GetRadixOf(outputFileName);


    // Check image physical space
    std::cout << "Checking images..." << std::flush;

    btk::ImageHelper< T >::IsInSamePhysicalSpace(inputImages);
    btk::ImageHelper< T >::IsInSamePhysicalSpace(modelImages);
    btk::ImageHelper< T >::IsInSamePhysicalSpace(inputImages[0], modelImages[0]);
    btk::ImageHelper< T,MaskImage >::IsInSamePhysicalSpace(inputImages[0], mask);

    std::cout << "done." << std::endl;

    std::cout << "Processing..." << std::endl;


    // Compute the mean image
    std::cout << std::setprecision(2) << "\t-> Compute mean image: " << percent << " %\r" << std::flush;

    typename T::Pointer mean = btk::ImageHelper< T >::CreateNewImageFromPhysicalSpaceOf(inputImages[0]);

    Iterator meanIt(mean, mean->GetLargestPossibleRegion());
    MaskIterator maskIt(mask, mask->GetLargestPossibleRegion());

    double normalizationCoefficient = 1.0 / static_cast< double >(inputImages.size());

    for(unsigned int i = 0; i < inputImages.size(); i++)
    {
        Iterator sampleIt(inputImages[i], inputImages[i]->GetLargestPossibleRegion());

        for(meanIt.GoToBegin(), sampleIt.GoToBegin(), maskIt.GoToBegin(); !meanIt.IsAtEnd() && !sampleIt.IsAtEnd() && !maskIt.IsAtEnd(); ++meanIt, ++sampleIt, ++maskIt)
        {
            if(maskIt.Get() > 0)
            {
                meanIt.Set(meanIt.Get() + sampleIt.Get() * normalizationCoefficient);
            }
        }

        percent += percentStep;
        std::cout << "\t-> Compute mean image: " << percent << " %\r" << std::flush;
    }

    std::cout << "\t-> Compute mean image: 100 %" << std::endl;


    // Compute the total sum of squares
    percent = 0.0;
    std::cout << "\t-> Compute total sum of squares: " << percent << " %\r" << std::flush;

    IntensitiesImage::Pointer totalSumSquares = btk::ImageHelper< T,IntensitiesImage >::CreateNewImageFromPhysicalSpaceOf(inputImages[0]);

    IntensitiesIterator totalSumSquaresIt(totalSumSquares, totalSumSquares->GetLargestPossibleRegion());

    for(unsigned int i = 0; i < inputImages.size(); i++)
    {
        Iterator sampleIt(inputImages[i], inputImages[i]->GetLargestPossibleRegion());

        for(totalSumSquaresIt.GoToBegin(), meanIt.GoToBegin(), sampleIt.GoToBegin(), maskIt.GoToBegin(); !totalSumSquaresIt.IsAtEnd() && !meanIt.IsAtEnd() && !sampleIt.IsAtEnd() && !maskIt.IsAtEnd(); ++totalSumSquaresIt, ++meanIt, ++sampleIt, ++maskIt)
        {
            if(maskIt.Get() > 0)
            {
                typename T::PixelType difference = sampleIt.Get() - meanIt.Get();
                totalSumSquaresIt.Set(totalSumSquaresIt.Get() + difference * difference);
            }
        }

        percent += percentStep;
        std::cout << "\t-> Compute total sum of squares: " << percent << " %\r" << std::flush;
    }

    std::cout << "\t-> Compute total sum of squares: 100 %" << std::endl;

    btk::ImageHelper< IntensitiesImage >::WriteImage(totalSumSquares, outputPrefix+"_total-error.nii.gz");


    // Compute the residuals sum of squares
    percent = 0.0;
    std::cout << "\t-> Compute residuals sum of squares: " << percent << " %\r" << std::flush;

    IntensitiesImage::Pointer residualsSumSquares = btk::ImageHelper< T,IntensitiesImage >::CreateNewImageFromPhysicalSpaceOf(inputImages[0]);

    IntensitiesIterator residualsSumSquaresIt(residualsSumSquares, residualsSumSquares->GetLargestPossibleRegion());

    for(unsigned int i = 0; i < inputImages.size(); i++)
    {
        Iterator sampleIt(inputImages[i], inputImages[i]->GetLargestPossibleRegion());
        Iterator modelIt(modelImages[i], modelImages[i]->GetLargestPossibleRegion());

        for(residualsSumSquaresIt.GoToBegin(), modelIt.GoToBegin(), sampleIt.GoToBegin(), maskIt.GoToBegin(); !residualsSumSquaresIt.IsAtEnd() && !modelIt.IsAtEnd() && !sampleIt.IsAtEnd() && !maskIt.IsAtEnd(); ++residualsSumSquaresIt, ++modelIt, ++sampleIt, ++maskIt)
        {
            if(maskIt.Get() > 0)
            {
                typename T::PixelType difference = sampleIt.Get() - modelIt.Get();
                residualsSumSquaresIt.Set(residualsSumSquaresIt.Get() + difference * difference);
            }
        }

        percent += percentStep;
        std::cout << "\t-> Compute residuals sum of squares: " << percent << " %\r" << std::flush;
    }

    std::cout << "\t-> Compute residuals sum of squares: 100 %" << std::endl;

    btk::ImageHelper< IntensitiesImage >::WriteImage(residualsSumSquares, outputPrefix+"_residual-error.nii.gz");


    // Compute coefficient of determination
    percent = 0.0;
    percentStep = 100.0 / static_cast< double >(inputImages[0]->GetLargestPossibleRegion().GetNumberOfPixels());
    std::cout << "\t-> Compute coefficient of determination: " << percent << " %\r" << std::flush;

    IntensitiesImage::Pointer coefficientOfDetermination = btk::ImageHelper< T,IntensitiesImage >::CreateNewImageFromPhysicalSpaceOf(inputImages[0]);

    IntensitiesIterator coefficientOfDeterminationIt(coefficientOfDetermination, coefficientOfDetermination->GetLargestPossibleRegion());

    for(coefficientOfDeterminationIt.GoToBegin(), totalSumSquaresIt.GoToBegin(), residualsSumSquaresIt.GoToBegin(), maskIt.GoToBegin(); !coefficientOfDeterminationIt.IsAtEnd() && !totalSumSquaresIt.IsAtEnd() && !residualsSumSquaresIt.IsAtEnd() && !maskIt.IsAtEnd(); ++coefficientOfDeterminationIt, ++totalSumSquaresIt, ++residualsSumSquaresIt, ++maskIt)
    {
        if(maskIt.Get() > 0)
        {
            IntensitiesImage::PixelType ratio = residualsSumSquaresIt.Get() / totalSumSquaresIt.Get();
            coefficientOfDeterminationIt.Set(1 - ratio);
        }
    }

    std::cout << "\t-> Compute coefficient of determination: 100 %" << std::endl;

    std::cout << "done." << std::endl;


    // Write output image
    btk::ImageHelper< IntensitiesImage >::WriteImage(coefficientOfDetermination, outputFileName);
}

