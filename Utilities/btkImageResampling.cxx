/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 27/01/2012
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


/* Standard includes */
#include "string"
#include "iomanip"
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"

/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPandoraBoxImageFilters.h"

int main(int argc, char *argv[])
{
  try{

    TCLAP::CmdLine cmd("Resample an image using: 1) a reference image (better option) or 2) specific size or spacing.", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputImageArg    ("i","input_image","input image file (will be casted in float)",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outputImageArg   ("o","output_image","output image file",true,"","string",cmd);
    TCLAP::ValueArg<std::string> referenceImageArg("r","reference_image","reference image file",false,"","string",cmd);
    TCLAP::MultiArg<float>       spacingArg       ("","spacing","resolution (spacing)) in mm (1 by default)",false,"float",cmd);
    TCLAP::MultiArg<float>       sizeArg          ("","size","size in voxels",false,"float");
    TCLAP::ValueArg<int>         orderArg         ("","order","order of the bspline interpolator (default:1)",false,1,"int",cmd);
             
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg ----------------------------------------------------
    std::string input_file       = inputImageArg.getValue();
    std::string output_file      = outputImageArg.getValue();
    std::string reference_file   = referenceImageArg.getValue();
            
    std::vector<float> size      = sizeArg.getValue();
    std::vector<float> spacing   = spacingArg.getValue();
    int order                    = orderArg.getValue();
    
    //ITK declaration ----------------------------------------------------------------------
    typedef float PixelType;
    const   unsigned int        Dimension = 3;
    
    typedef itk::Image< PixelType, Dimension >    itkImage;
  
    std::cout<<"Reading the input image:"<<input_file<<"\n";
    itkImage::Pointer inputImage = btk::ImageHelper< itkImage > ::ReadImage(input_file);

    itkImage::Pointer outputImage;

    if(reference_file != ""){
      std::cout<<"Using the following reference image:"<<reference_file<<"\n";
      itkImage::Pointer referenceImage = btk::ImageHelper< itkImage > ::ReadImage(reference_file);
      btk::PandoraBoxImageFilters::ResampleImageUsingReference(inputImage,outputImage, referenceImage, order);
    }
    else{
      if( (size.size() < Dimension ) && (spacing.size() < Dimension ) ){
        std::cout<<"Please provide at least size or scaling factors for each dimension.\n";
        return 0;
      }
      itkImage::SizeType itkSize;
      itkImage::SpacingType itkSpacing;

      if( size.size() == Dimension ){
        for(unsigned int i=0; i<Dimension; i++)
          itkSize[i] = size[i];
        btk::PandoraBoxImageFilters::ResampleImageUsingSize(inputImage,outputImage, itkSize, order);
      }
      if( spacing.size() == Dimension ){
        for(uint i=0; i<Dimension; i++)
          itkSpacing[i] = spacing[i];
        btk::PandoraBoxImageFilters::ResampleImageUsingSpacing(inputImage,outputImage, itkSpacing, order);
      }
    }

    btk::ImageHelper<itkImage>::WriteImage(outputImage, output_file);
    
    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;

}
