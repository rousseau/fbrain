/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 15/11/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr)

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
#include "string"
#include "algorithm"

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"


// Define label image type
const unsigned int Dimension = 3;

typedef short                                       PixelType;
typedef itk::Image< PixelType,Dimension >           LabelImage;
typedef itk::ImageRegionIterator< LabelImage >      LabelImageIterator;
typedef itk::ImageRegionConstIterator< LabelImage > LabelImageConstIterator;


int main(int argc, char *argv[])
{
    try
    {
        //
        // Read command line
        //

        TCLAP::CmdLine cmd("Splits a label image into binary components", ' ', "2.0");

        TCLAP::ValueArg< std::string >  inputFileNameArg("i", "input", "Input image", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output image", true, "", "string", cmd);
        //TCLAP::ValueArg< short >                labelArg("l", "label", "Label value", true, 0, "short", cmd);
	TCLAP::MultiArg< short >                labelArg("l", "label", "Label value", true, "short", cmd);

        // Parse the command line
        cmd.parse( argc, argv );

        std::string inputFileName  = inputFileNameArg.getValue();
        std::string outputFileName = outputFileNameArg.getValue();
        //short                label = labelArg.getValue();
	std::vector<short> labels = labelArg.getValue();


        //
        // Read input image
        //
        LabelImage::ConstPointer inputImage = btk::ImageHelper< LabelImage >::ReadConstImage(inputFileName);


        //
        // Create new image from input
        //
        LabelImage::Pointer outputImage = btk::ImageHelper< LabelImage >::CreateNewImageFromPhysicalSpaceOfConst(inputImage);


        //
        // Binarize label
        //
        LabelImageIterator outputIt(outputImage, outputImage->GetLargestPossibleRegion());
        LabelImageConstIterator inputIt(inputImage, inputImage->GetLargestPossibleRegion());

        for(inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd() && !outputIt.IsAtEnd(); ++inputIt, ++outputIt)
        {
	  if( std::find(labels.begin(), labels.end(), inputIt.Get()) != labels.end() ) 
            {
                outputIt.Set(1);
            }
            else // inputIt.Get() != label
            {
                outputIt.Set(0);
            }
        } // for each voxel


        //
        // Write output image
        //

        btk::ImageHelper< LabelImage >::WriteImage(outputImage, outputFileName);


        return EXIT_SUCCESS;

    }
    catch(TCLAP::ArgException &e)
    {
        btkCoutMacro("Exception: " << e.error() << " for arg " << e.argId());
    }
    catch(std::string &message)
    {
        btkCoutMacro("Exception: " << message);
    }
}

