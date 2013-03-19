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
#include <tclap/CmdLine.h>

// STL includes
#include "iostream"
#include "string"
#include "cstdlib"

// ITK includes
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"

// Local includes
#include "btkImageHelper.h"


// Types definitions
typedef float ProbabilityMapPixelType;
typedef short TissueSegmentationPixelType;
const unsigned int ImageDimension = 3;

typedef itk::Image< ProbabilityMapPixelType,ImageDimension >                   ProbabilityMap;
typedef itk::ResampleImageFilter< ProbabilityMap,ProbabilityMap >              ProbabilityMapResampleFilter;
typedef itk::AddImageFilter< ProbabilityMap,ProbabilityMap,ProbabilityMap >    ProbabilityMapAddFilter;
typedef itk::DivideImageFilter< ProbabilityMap,ProbabilityMap,ProbabilityMap > ProbabilityMapDivideFilter;

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
        TCLAP::MultiArg< std::string > inputLabelFileNamesArg("i", "input", "Input label filenames", true, "string", cmd);
        TCLAP::ValueArg< std::string > outputTissueFileNameArg("o", "output", "Output filename", false, "outputTissue.nii.gz", "string", cmd);

        // Parse arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::vector< std::string > inputLabelFileNames = inputLabelFileNamesArg.getValue();
        std::string               outputTissueFileName = outputTissueFileNameArg.getValue();


        //
        // Load images
        //

        std::vector< ProbabilityMap::Pointer > inputs = btk::ImageHelper< ProbabilityMap >::ReadImage(inputLabelFileNames);

        //
        // Checking physical space
        //

        if(!btk::ImageHelper< ProbabilityMap >::IsInSamePhysicalSpace(inputs))
        {
            throw std::string("The input images do not have the same physical space !");
        }


        //
        // Processing image
        //

        std::cout << "Processing images..." << std::endl;
/*
        // Normalize images
        ProbabilityMap::Pointer normalizationImage = btk::ImageHelper< ProbabilityMap >::CreateNewImageFromPhysicalSpaceOf(inputs[0].GetPointer());

        for(unsigned int i = 0; i < inputs.size(); i++)
        {
            ProbabilityMapAddFilter::Pointer addFilter = ProbabilityMapAddFilter::New();
            addFilter->SetInput1(normalizationImage);
            addFilter->SetInput2(inputs[i]);
            addFilter->Update();

            normalizationImage = addFilter->GetOutput();
        } // for each image

        for(unsigned int i = 0; i < inputs.size(); i++)
        {
            ProbabilityMapDivideFilter::Pointer divideFilter = ProbabilityMapDivideFilter::New();
            divideFilter->SetInput1(inputs[i]);
            divideFilter->SetInput2(normalizationImage);
            divideFilter->Update();

            inputs[i] = divideFilter->GetOutput();
        } // for each image
*/

        // Create a new tissue map
        TissueSegmentation::Pointer outputTissues = btk::ImageHelper< ProbabilityMap,TissueSegmentation >::CreateNewImageFromPhysicalSpaceOf(inputs[0].GetPointer());


        // Compute new tissue map by majority voting
        ProbabilityMap::SizeType size = inputs[0]->GetLargestPossibleRegion().GetSize();

        for(unsigned int z = 0; z < size[2]; z++)
        {
            for(unsigned int x = 0; x < size[0]; x++)
            {
                for(unsigned int y = 0; y < size[1]; y++)
                {
                    ProbabilityMap::IndexType index;
                    index[0] = x; index[1] = y; index[2] = z;

                    ProbabilityMap::PixelType maxValue = inputs[0]->GetPixel(index);
                    TissueSegmentation::PixelType label = static_cast< TissueSegmentation::PixelType >(1);

                    for(unsigned int i = 1; i < inputs.size(); i++)
                    {
                        ProbabilityMap::PixelType value = inputs[i]->GetPixel(index);

                        if(value > maxValue)
                        {
                            maxValue = value;
                            label    = static_cast< TissueSegmentation::PixelType >(i+1);
                        }
                    } // for each image
										if(fabs(maxValue)<0.0001)
											label = 0;
                    outputTissues->SetPixel(index, label);
                } // for each voxel
            } // for each line
        } // for each slice


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
