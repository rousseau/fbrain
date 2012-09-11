/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 03/07/2012
  Author(s): Youssef Taleb

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


/* ITK */
#include "itkImage.h"

/* BTK */
#include "btkImageHelper.h"
#include "btkImageIntersectionCalculator.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>

int main(int argc, char * argv[])
{
    /* Typedefs */

    const unsigned int Dimension = 3;
    typedef short PixelType;
    typedef itk::Image<PixelType, Dimension>            ImageType;
    typedef btk::ImageIntersectionCalculator<ImageType> IntersectionCalculatorType;
    typedef itk::Image<unsigned char, Dimension>        ImageMaskType;


    IntersectionCalculatorType::Pointer intersector=IntersectionCalculatorType::New();

    //TCLAP Commands for arguments
    TCLAP::CmdLine cmd("Extract masks using intersections of input images (Bounding Box)", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","input image files",true,"string",cmd);
    TCLAP::MultiArg<std::string> outArg  ("o","output","output mask files",true,"string",cmd);


    std::vector< std::string > inputFileImage;
    std::vector< std::string > outputFileImage;
    std::vector< ImageType::Pointer > inputImage;

    // Parse the argv array.
    cmd.parse( argc, argv );
    inputFileImage = inputArg.getValue();
    outputFileImage = outArg.getValue();


    std::cout<<"Initialization... "<<std::endl<<std::endl;


    // Error output if wrong argument number
    if(inputFileImage.size() == 0)
    {
        std::cout<<"You should have at least one image in input  !"<<std::endl;
        return EXIT_FAILURE;
    }

    // Read input arguments into input image
    inputImage = btk::ImageHelper<ImageType>::ReadImage(inputFileImage);

    // Intersection Process
    for (unsigned int i=0;i<inputFileImage.size();i++)
    {
        intersector->AddImage(inputImage[i]) ;
    }

    intersector->Compute();
    for (int i=0;i<outputFileImage.size();i++)
    {
        ImageMaskType* mask = intersector-> GetImageMask(i);
        btk::ImageHelper<ImageMaskType>::WriteImage(mask,outputFileImage[i]);
    }


    std::cout<<"End of execution..."<<std::endl;

    return EXIT_SUCCESS;

}
