/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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
#include "tclap/CmdLine.h"

#include "itkImageIOBase.h"
#include "itkImageIOFactory.h"
#include "itkDisplacementFieldTransform.h"

#include "btkImageHelper.h"
#include "btkResolutionSamplerFilter.h"


typedef short PixelType;
typedef itk::Image<PixelType,3> TImage;

int main(int argc, char *argv[])
{
    try
    {

        //////////////////////////////////////////////////////////////////////////
        //
        // Parse program's arguments
        //

        // Define command line object for program
        TCLAP::CmdLine cmd("Resolution sampler.", ' ', "1.0");

        // Define arguments
        TCLAP::ValueArg< std::string >  inputNameArg("i", "input_image", "Input diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >  outputFileNameArg("o", "output_image", "Output diffusion sequence", true, "", "string", cmd);

        TCLAP::ValueArg<unsigned int>  resolutionFactorArg("","resolution","mm side length of voxel",true,1,"int",cmd);

        // Parse command line
        cmd.parse(argc, argv);



        // Get back arguments' values
        std::string inputFileName       = inputNameArg.getValue();
        std::string outputFileName      = outputFileNameArg.getValue();

        unsigned int resolutionFactor     = resolutionFactorArg.getValue();

        //
        // Processing
        //

        // Check image informations
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFileName.c_str(), itk::ImageIOFactory::ReadMode);
        imageIO->SetFileName(inputFileName);
        imageIO->ReadImageInformation();

        int Dimension = imageIO->GetNumberOfDimensions();
        btkCoutMacro("Dimension of input image :" << Dimension);

        TImage::Pointer inputImage = btk::ImageHelper< TImage >::ReadImage(inputFileName);
        typedef btk::ResolutionSamplerFilter<TImage> ResolutionFilterType;
        ResolutionFilterType::Pointer ResolutionSampler = ResolutionFilterType::New();
        ResolutionSampler -> SetInputImage(inputImage.GetPointer());
        ResolutionSampler -> SetResolution(resolutionFactor);
        ResolutionSampler -> Update();


        btk::ImageHelper<  TImage >::WriteImage(ResolutionSampler->GetOutput(),outputFileName);



    }
    catch(itk::ExceptionObject &exception)
    {
        std::cerr << "ITK error:" << std::endl;
        std::cerr << exception << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
