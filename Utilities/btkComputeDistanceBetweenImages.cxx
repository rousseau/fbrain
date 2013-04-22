/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 17/04/2013
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

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"


template< class TImage >
void ProcessScalar(std::vector< std::string > &inputFileNames)
{
    // Read images
    std::vector< typename TImage::Pointer > inputImages = btk::ImageHelper< TImage >::ReadImage(inputFileNames);

    // Verify sizes of images and resample on largest lattice
    if(!btk::ImageHelper< TImage >::IsInSamePhysicalSpace(inputImages))
    {
        btkException("Process: Images are not in the same physical space !");
    }

    // Compute weighted sum
    btkCoutMacro("Computing the euclidian distance... ");

    typedef itk::ImageRegionConstIterator< TImage > ImageIterator;

    ImageIterator it1(inputImages[0], inputImages[0]->GetLargestPossibleRegion());
    ImageIterator it2(inputImages[1], inputImages[1]->GetLargestPossibleRegion());

    double distance = 0.0;

    for(it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd() && !it2.IsAtEnd(); ++it1, ++it2)
    {
        double difference = it1.Get() - it2.Get();
        distance += difference*difference;
    }

    btkCoutMacro("done.");

    // Write output
    std::cerr << distance / inputImages[0]->GetLargestPossibleRegion().GetNumberOfPixels() << std::endl;
}

void ProcessVector(std::vector< std::string > &inputFileNames)
{
    typedef itk::DisplacementFieldTransform< float,3 >::DisplacementFieldType TImage;

    // Read images
    std::vector< TImage::Pointer > inputImages = btk::ImageHelper< TImage >::ReadImage(inputFileNames);

    // Verify sizes of images and resample on largest lattice
    if(!btk::ImageHelper< TImage >::IsInSamePhysicalSpace(inputImages))
    {
        btkException("Process: Images are not in the same physical space !");
    }

    // Compute weighted sum
    btkCoutMacro("Computing the euclidian distance... ");

    typedef itk::ImageRegionConstIterator< TImage > ImageIterator;

    ImageIterator it1(inputImages[0], inputImages[0]->GetLargestPossibleRegion());
    ImageIterator it2(inputImages[1], inputImages[1]->GetLargestPossibleRegion());

    TImage::PixelType distance = itk::NumericTraits< TImage::PixelType >::ZeroValue();

    for(it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd() && !it2.IsAtEnd(); ++it1, ++it2)
    {
        TImage::PixelType difference = it1.Get() - it2.Get();
        distance += difference*difference;
    }

    btkCoutMacro("done.");

    double distanceValue = distance.GetNorm();

    // Write output
    std::cerr << distanceValue / inputImages[0]->GetLargestPossibleRegion().GetNumberOfPixels() << std::endl;
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
        TCLAP::CmdLine cmd("Compute the euclidian distance between a pair of 3D images", ' ', "1.0", true);

        // Arguments
        TCLAP::UnlabeledMultiArg< std::string > inputFileNamesArg("input", "Input image filenames (2 files)", true, "string", cmd);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();


        //
        // Processing
        //

        // Verify number of images and number of weights
        if(inputFileNames.size() != 2)
        {
            btkException("Main: There should be 2 images !");
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

                ProcessVector(inputFileNames);
                break;

            case itk::ImageIOBase::SCALAR:

                switch(imageIO->GetComponentType())
                {
                    case itk::ImageIOBase::SHORT:
                        ProcessScalar< itk::Image< short,3 > >(inputFileNames);
                        break;

                    case itk::ImageIOBase::FLOAT:
                        ProcessScalar< itk::Image< float,3 > >(inputFileNames);
                        break;

                    case itk::ImageIOBase::DOUBLE:
                        ProcessScalar< itk::Image< double,3 > >(inputFileNames);
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
