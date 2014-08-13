/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 05/10/2011
 Author(s): François Rousseau (rousseau@unistra.fr)
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
#include "string"
#include "vector"
#include "numeric"

// TCLAP includes
#include <tclap/CmdLine.h>

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkDisplacementFieldTransform.h"
#include "itkResampleImageFilter.h"

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"
#include "btkResampleImagesToBiggestImageFilter.h"
#include "btkWeightedSumOfImagesFilter.h"


// Weight sum over template images
template< class TImage >
void Process(std::vector< std::string > &inputFileNames, std::vector< float > &weights, std::string &outputFileName)
{
    // Read images
    std::vector< typename TImage::Pointer > inputImages = btk::ImageHelper< TImage >::ReadImage(inputFileNames);

    // Verify sizes of images and resample on largest lattice
    if(!btk::ImageHelper< TImage >::IsInSamePhysicalSpace(inputImages))
    {
        btkCoutMacro("Resampling images in same physical space... ");

        typename btk::ResampleImagesToBiggestImageFilter< TImage >::Pointer resampleFilter = btk::ResampleImagesToBiggestImageFilter< TImage >::New();
        resampleFilter->SetInputs(inputImages);
        resampleFilter->Update();
        inputImages = resampleFilter->GetOutputs();

        btkCoutMacro("done.");
    }

    // Compute weighted sum
    btkCoutMacro("Summing images... ");

    typename btk::WeightedSumOfImagesFilter< TImage >::Pointer sumFilter = btk::WeightedSumOfImagesFilter< TImage >::New();
    sumFilter->SetInputs(inputImages);
    sumFilter->SetWeights(weights);
    sumFilter->Update();

    btkCoutMacro("done.");

    // Write output image
    btk::ImageHelper< TImage >::WriteImage(sumFilter->GetOutput(), outputFileName);
}


// Main function
int main(int argc, char *argv[])
{ 
    try
    {

        //
        // Command line parser
        //

        // Command line
        TCLAP::CmdLine cmd("Compute a weighted mean of 3D images", ' ', "2.0", true);
    
        // Arguments
        TCLAP::MultiArg< std::string > inputFileNamesArg("i", "input", "Input image filenames", true, "string", cmd);
        TCLAP::MultiArg< float >              weightsArg("w", "weight", "Image weight (default: 1/N)", false, "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output image filename (default: \"out.nii.gz\")", false, "out.nii.gz", "string", cmd);

        TCLAP::SwitchArg  noNormalizationArg("", "no_normalization", "Desactivate the normalization of the weights (default: false -> meaning that there is a normalization of weights)", cmd, false);
    
        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();
        std::vector< float >              weights = weightsArg.getValue();
        std::string                outputFileName = outputFileNameArg.getValue();

        bool  noNormalization = noNormalizationArg.getValue();
 

        //
        // Processing
        //

        // Verify number of images and number of weights
        if(weights.size() > 0 && inputFileNames.size() != weights.size())
        {
            btkException("Main: The number of input images is different than the number of weights !");
        }
        else if(weights.empty())
        {
            weights = std::vector< float >(inputFileNames.size(), 1.0f/static_cast< float >(inputFileNames.size()));
        }

        // Normalize weights if needed
        if(!noNormalization)
        {
            float sum = std::accumulate(weights.begin(), weights.end(), 0.0f);

            for(int i = 0; i < weights.size(); i++)
            {
                weights[i] /= sum;
            }
        }

        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFileNames[0].c_str(), itk::ImageIOFactory::ReadMode);
        imageIO->SetFileName(inputFileNames[0]);
        imageIO->ReadImageInformation();

        if(imageIO->GetNumberOfDimensions() != 3)
        {
            btkCerrMacro("Error: Unsupported image dimension !");
            exit(EXIT_FAILURE);
        }

        // Switch on pixel type
        switch(imageIO->GetPixelType())
        {
            case itk::ImageIOBase::VECTOR:

                Process< itk::DisplacementFieldTransform< float,3 >::DisplacementFieldType >(inputFileNames, weights, outputFileName);
                break;

            case itk::ImageIOBase::SCALAR:

                switch(imageIO->GetComponentType())
                {
                    case itk::ImageIOBase::SHORT:
                        Process< itk::Image< short,3 > >(inputFileNames, weights, outputFileName);
                        break;
                        
                    case itk::ImageIOBase::USHORT:
                        Process< itk::Image< unsigned short,3 > >(inputFileNames, weights, outputFileName);
                        break;

                    case itk::ImageIOBase::FLOAT:
                        Process< itk::Image< float,3 > >(inputFileNames, weights, outputFileName);
                        break;

                    case itk::ImageIOBase::DOUBLE:
                        Process< itk::Image< double,3 > >(inputFileNames, weights, outputFileName);
                        break;

                    default:
                        btkCerrMacro("Unsupported component type !");
                        exit(EXIT_FAILURE);
                }

                break;

            default:
                btkCerrMacro("Error: Unsupported pixel type !");
                exit(EXIT_FAILURE);
        }
    }
    catch(TCLAP::ArgException &e)
    {
        btkCoutMacro("Exception: " << e.error() << " for arg " << e.argId());
    }
    catch(itk::ExceptionObject &object)
    {
        btkCoutMacro(object);
    }
  
  return EXIT_SUCCESS;
}
