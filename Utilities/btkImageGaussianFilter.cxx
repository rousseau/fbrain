/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 01/06/2012
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
#include "cstdlib"
#include "string"

// ITK includes
#include "itkImage.h"
#include "itkDiscreteGaussianImageFilter.h"

// BTK includes
#include "btkImageHelper.h"


// Define label image type
const unsigned int Dimension = 3;

typedef unsigned short LabelPixelType;
typedef float          ProbabilityPixelType;

typedef itk::Image< LabelPixelType,Dimension >       LabelImage;
typedef itk::Image< ProbabilityPixelType,Dimension > ProbabilityImage;

typedef itk::DiscreteGaussianImageFilter< LabelImage,ProbabilityImage > ImageGaussianFilter;


int main(int argc, char *argv[])
{
    try
    {
        //
        // Command line arguments
        //

        // Set up command line parser
        TCLAP::CmdLine cmd("Apply a gaussian filter on an image (convolutive filter)", ' ', "1.0", true);

        // Set up arguments
        TCLAP::ValueArg< std::string > inputFileNameArg("i", "input_image", "Input image filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output_image", "Output image filename", false, "Blurred-image.nii.gz", "string", cmd);

        // Parse arguments
        cmd.parse(argc, argv);

        // Get the parsed values
        std::string inputFileName  = inputFileNameArg.getValue();
        std::string outputFileName = outputFileNameArg.getValue();


        //
        // Read image
        //

        LabelImage::Pointer inputImage = btk::ImageHelper< LabelImage >::ReadImage(inputFileName);


        //
        // Process image
        //

        ImageGaussianFilter::Pointer filter = ImageGaussianFilter::New();
        filter->SetInput(inputImage);

        ImageGaussianFilter::ArrayType variance;
        variance[0] = 1.0;
        variance[1] = 1.0;
        variance[2] = 1.0;
        filter->SetVariance(variance);

        filter->SetUseImageSpacingOff();
        filter->Update();

        ProbabilityImage::Pointer outputImage = filter->GetOutput();


        //
        // Write image
        //

        btk::ImageHelper< ProbabilityImage >::WriteImage(outputImage, outputFileName);
    }
    catch(std::string &error)
    {

    }

    return EXIT_SUCCESS;
}
