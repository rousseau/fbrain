/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 06/03/2012
 Author(s): François Rousseau (rousseau@unistra.fr)
 
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
#include "iomanip"

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// TCLAP includes
#include <tclap/CmdLine.h>

// Local includes
#include "btkImageHelper.h"


// Type definitions
typedef float PixelType;
const unsigned int Dimension = 3;

typedef itk::Image< PixelType,Dimension >          Image;
typedef itk::ImageRegionIteratorWithIndex< Image > ImageIterator;


int main(int argc, char *argv[])
{
    try
    {
        //
        // Command line processing
        //

        TCLAP::CmdLine cmd("Normalize in a voxelwise manner a set of (positive) probability maps", ' ', "1.0", true);
    
        TCLAP::MultiArg<std::string> inputImageArg("i", "input_image", "Input images' filenames", true, "string", cmd);
        TCLAP::MultiArg<std::string> outputImageArg("o", "output_image"," Output images' filenames", true, "string", cmd);
        TCLAP::ValueArg<float> thresholdArg  ("t", "threshold", "Threshold on the sum (for each voxel) to avoid calculation error) (default = 0.00001)", false, 0.00001, "float", cmd);
    
        // Parse the args.
        cmd.parse( argc, argv );

        // Get the value parsed by each arg.
        std::vector<std::string> inputFileNames  = inputImageArg.getValue();
        std::vector<std::string> outputFileNames = outputImageArg.getValue();
        float threshold                          = thresholdArg.getValue();

        // Verify data
        if(inputFileNames.size() != outputFileNames.size())
            throw(std::string("Error: the number of input filenames is not the same as the number of output filenames !"));
     

        //
        // Reading images
        //

        std::vector< Image::Pointer > inputImages  = btk::ImageHelper< Image >::ReadImage(inputFileNames);
        std::vector< Image::Pointer > outputImages = btk::ImageHelper< Image >::CreateNewImageFromPhysicalSpaceOf(inputImages);
    

        //
        // Processing
        //

        ImageIterator it(inputImages[0], inputImages[0]->GetLargestPossibleRegion());

        for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            double sum = 0;
            Image::IndexType index = it.GetIndex();

            for(int i=0; i<inputFileNames.size(); i++)
            {
                if(inputImages[i]->GetPixel(index) > 0)
                {
                    sum += inputImages[i]->GetPixel(index);
                }
            }

            if(sum > threshold)
            {
                for(int i=0; i<outputFileNames.size(); i++)
                {
                    if(inputImages[i]->GetPixel(index) > 0)
                    {
                        outputImages[i]->SetPixel(index,inputImages[i]->GetPixel(index) / sum);
                    }
                }
            }
        }


        //
        // Write images
        //

        btk::ImageHelper< Image >::WriteImage(outputImages, outputFileNames);
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




