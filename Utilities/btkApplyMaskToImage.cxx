/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 26/06/2012
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
#include "itkMaskImageFilter.h"

/* BTK */
#include "btkImageHelper.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>


int main(int argc, char * argv[])
{
    const unsigned int Dimension = 3;
    typedef short PixelType;
    typedef itk::Image<PixelType, Dimension> itkImage;
    typedef itk::Image<unsigned char, Dimension> itkImageMask;
    typedef itk::MaskImageFilter<itkImage,itkImageMask, itkImage> MaskImageFilter;

    TCLAP::CmdLine cmd("btkApplyMaskToImage: Apply a mask(non-zero values) to a image", ' ', "1.0", true);
    TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string", cmd);
    TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string", cmd);
    TCLAP::ValueArg<std::string> inputMaskArg("m","mask_file","filename of the mask image (dimension = 3)",false,"","string", cmd);

    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::string input_file       = inputImageArg.getValue();
    std::string output_file      = outputImageArg.getValue();
    std::string mask_file        = inputMaskArg.getValue();

    itkImage::Pointer image = itkImage::New();
    itkImageMask::Pointer mask = itkImageMask::New();

    //read input image and mask
    image = btk::ImageHelper<itkImage>::ReadImage(input_file);
    mask = btk::ImageHelper<itkImageMask>::ReadImage(mask_file);

    MaskImageFilter::Pointer maskImageFilter = MaskImageFilter::New();

    maskImageFilter->SetInput(image);
    maskImageFilter->SetMaskImage(mask);
    maskImageFilter->Update();

    //write the result
    btk::ImageHelper<itkImage>::WriteImage(maskImageFilter->GetOutput(),output_file);

    return EXIT_SUCCESS;
}
