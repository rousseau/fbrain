/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 25/06/2012
  Author(s): Marc Schweitzer (marc.schweitzer@unistra.fr)

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
#include "itkSubtractImageFilter.h"

/* BTK */

#include "btkImageHelper.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>


int main(int argc, char * argv[])
{
    /* Typedefs */

    const unsigned int Dimension = 3;
    typedef float PixelType;

    typedef itk::Image<PixelType, Dimension> itkImage;

    typedef itk::SubtractImageFilter<itkImage> SubstractFilter;

    //TCLAP

    TCLAP::CmdLine cmd("Pixel-wise subtraction of two image, or one image and a constant (-i im1 -i im2 -o result or -i im1 -c value -o result", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",true,"","string",cmd);
    TCLAP::ValueArg<float> cstArg  ("c","constant","Super resolution output image",false,-9999,"string",cmd);

    std::vector< std::string > inputFileImages;
    std::string outputFileImage;
    std::vector< itkImage::Pointer > inputsImages;
    itkImage::Pointer outputImage;

    // Parse the argv array.
    cmd.parse( argc, argv );
    inputFileImages = inputArg.getValue();
    outputFileImage = outArg.getValue().c_str();
    float cst = cstArg.getValue();

    std::cout<<"Begin the subtraction : "<<std::endl<<std::endl;

    if(inputFileImages.size() != 2 && cst == -9999)
    {
        std::cout<<"You should have 2 image in input or 1 image and a constant value(-c) !"<<std::endl;
        return EXIT_FAILURE;
    }
    else if(cst != -9999)
    {
        std::cout<<"Constant value = "<<cst<<std::endl;
    }


    inputsImages = btk::ImageHelper<itkImage>::ReadImage(inputFileImages);

    SubstractFilter::Pointer subtractFilter = SubstractFilter::New ();
    subtractFilter->SetInput1(inputsImages[0]);

    if(inputFileImages.size() == 1 && cst != -9999)
    {
        subtractFilter->SetConstant2(cst);
    }
    else
    {
        subtractFilter->SetInput2(inputsImages[1]);
    }

    subtractFilter->Update();

    outputImage = subtractFilter->GetOutput();

    btk::ImageHelper<itkImage>::WriteImage(outputImage,outputFileImage);

    std::cout<<"End Subtraction..."<<std::endl;

    return EXIT_SUCCESS;











}
