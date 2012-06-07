/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 07/06/2012
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


// TCLAP includes
#include "CmdLine.h"

// STL includes
#include "iostream"
#include "string"
#include "cstdlib"

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"

// Local includes
#include "btkImageHelper.h"


// Types definitions
typedef float ProbabilityMapPixelType;
typedef unsigned short TissueSegmentationPixelType;
const unsigned int ImageDimension = 3;

typedef itk::Image< ProbabilityMapPixelType,ImageDimension > ProbabilityMap;
typedef itk::ImageRegionIterator< ProbabilityMap >           ProbabilityMapIterator;

typedef itk::Image< TissueSegmentationPixelType,ImageDimension > TissueSegmentation;
typedef itk::ImageRegionIterator<TissueSegmentation>             TissueSegmentationIterator;


int main(int argc, char *argv[])
{
    try
    {
        //
        // Parse program's arguments
        //

        // Define command line parser
        TCLAP::CmdLine cmd("Binarize tissue probability maps", ' ', "1.0");

        // Define command line arguments
        TCLAP::ValueArg< std::string > inputGMFileNameArg("g", "grey_matter", "Input grey matter filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > inputWMFileNameArg("w", "white_matter", "Input white matter filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > inputCerveletFileNameArg("c", "cervelet", "Input cervelet filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > inputBrainstemFileNameArg("b", "brainstem", "Input brainstem filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > inputCSFFileNameArg("r", "csf", "Input CSF filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputTissueFileNameArg("o", "output", "Output filename", false, "outputTissue.nii.gz", "string", cmd);

        // Parse arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string inputGMFileName        = inputGMFileNameArg.getValue();
        std::string inputWMFileName        = inputWMFileNameArg.getValue();
        std::string inputCerveletFileName  = inputCerveletFileNameArg.getValue();
        std::string inputBrainstemFileName = inputBrainstemFileNameArg.getValue();
        std::string inputCSFFileName       = inputCSFFileNameArg.getValue();
        std::string outputTissueFileName   = outputTissueFileNameArg.getValue();


        //
        // Load images
        //

        ProbabilityMap::Pointer        inputGM = btk::ImageHelper< ProbabilityMap >::ReadImage(inputGMFileName);
        ProbabilityMap::Pointer        inputWM = btk::ImageHelper< ProbabilityMap >::ReadImage(inputWMFileName);
        ProbabilityMap::Pointer  inputCervelet = btk::ImageHelper< ProbabilityMap >::ReadImage(inputCerveletFileName);
        ProbabilityMap::Pointer inputBrainstem = btk::ImageHelper< ProbabilityMap >::ReadImage(inputBrainstemFileName);
        ProbabilityMap::Pointer       inputCSF = btk::ImageHelper< ProbabilityMap >::ReadImage(inputCSFFileName);


        //
        // Processing image
        //

        std::cout << "Processing images..." << std::endl;

        // Create a new tissue map
        TissueSegmentation::Pointer outputTissues = btk::ImageHelper< ProbabilityMap,TissueSegmentation >::CreateNewFromPhysicalSpaceOf(inputGM);

        // Iterators
        ProbabilityMapIterator     inputGMIterator(inputGM, inputGM->GetLargestPossibleRegion());
        ProbabilityMapIterator     inputWMIterator(inputWM, inputWM->GetLargestPossibleRegion());
        ProbabilityMapIterator     inputCerveletIterator(inputCervelet, inputCervelet->GetLargestPossibleRegion());
        ProbabilityMapIterator     inputBrainstemIterator(inputBrainstem, inputBrainstem->GetLargestPossibleRegion());
        ProbabilityMapIterator     inputCSFIterator(inputCSF, inputCSF->GetLargestPossibleRegion());
        TissueSegmentationIterator outputTissueIterator(outputTissues, outputTissues->GetLargestPossibleRegion());


        // Rebuild the tissue map by picking the maximum probability of each voxel
        for(inputGMIterator.GoToBegin(), inputWMIterator.GoToBegin(), inputCerveletIterator.GoToBegin(), inputBrainstemIterator.GoToBegin(), inputCSFIterator.GoToBegin(), outputTissueIterator.GoToBegin();
            !inputGMIterator.IsAtEnd() && !inputWMIterator.IsAtEnd() && !inputCerveletIterator.IsAtEnd() && !inputBrainstemIterator.IsAtEnd() && !inputCSFIterator.IsAtEnd() && !outputTissueIterator.IsAtEnd();
            ++inputGMIterator, ++inputWMIterator, ++inputCerveletIterator, ++inputBrainstemIterator, ++inputCSFIterator, ++outputTissueIterator)
        {
            // Normalization
            ProbabilityMapPixelType       sum = inputGMIterator.Get() + inputWMIterator.Get() + inputCerveletIterator.Get() + inputBrainstemIterator.Get() + inputCSFIterator.Get();
            ProbabilityMapPixelType        GM = inputGMIterator.Get() / sum;
            ProbabilityMapPixelType        WM = inputWMIterator.Get() / sum;
            ProbabilityMapPixelType  Cervelet = inputCerveletIterator.Get() / sum;
            ProbabilityMapPixelType Brainstem = inputBrainstemIterator.Get() / sum;
            ProbabilityMapPixelType       CSF = inputCSFIterator.Get() / sum;

            // 1: GM, 2: WM, 3: Cervelet, 4: Brainstem, 5: CSF
            ProbabilityMapPixelType maxValue = inputGMIterator.Get();
            outputTissueIterator.Set(1);

            if(inputWMIterator.Get() > maxValue)
            {
                maxValue = inputWMIterator.Get();
                outputTissueIterator.Set(2);
            }

            if(inputCerveletIterator.Get() > maxValue)
            {
                maxValue = inputCerveletIterator.Get();
                outputTissueIterator.Set(3);
            }

            if(inputBrainstemIterator.Get() > maxValue)
            {
                maxValue = inputBrainstemIterator.Get();
                outputTissueIterator.Set(4);
            }

            if(inputCSFIterator.Get() > maxValue)
            {
                outputTissueIterator.Set(5);
            }
        }

        std::cout << "done." << std::endl;


        //
        // Write output image
        //

        btk::ImageHelper< TissueSegmentation >::WriteImage(outputTissues, outputTissueFileName);
    }
    catch(TCLAP::ArgException &e)
    {
        std::cerr << "Error (arguments): " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch(itk::ImageFileReaderException &e)
    {
        std::cerr << "Error (input images): an image filename is wrong or does not exist ; program exit !" << std::endl;
        std::cerr << "ITK error trace:" << e;
        exit(EXIT_FAILURE);
    }
    catch(std::string &e)
    {
        std::cerr << std::endl << e << " ; program exit !" << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
