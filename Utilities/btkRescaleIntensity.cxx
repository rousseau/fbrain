/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 22/11/2011
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
#include "itkImageRegionIterator.h"


#include "btkImageHelper.h"


int main(int argc, char *argv[])
{
  
  try{

    TCLAP::CmdLine cmd("Rescale the intensity values of an image using short values", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputImageArg("i","input_file","input image file (short)",true,"","string", cmd);
    TCLAP::ValueArg<std::string> maskImageArg("m","mask_file","mask image file (short)",true,"","string", cmd);
    TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string", cmd);
    TCLAP::ValueArg< short > minArg("","min","minimum value for the rescaling (default is 0)",false,0,"short", cmd);
    TCLAP::ValueArg< short > maxArg("","max","maximum value for the rescaling (default is 255)",false,255,"short", cmd);
    
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::string input_file       = inputImageArg.getValue();
    std::string mask_file        = maskImageArg.getValue();
    std::string output_file      = outputImageArg.getValue();    
    
    short minValue              = minArg.getValue();
    short maxValue              = maxArg.getValue();
    
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::ImageRegionIterator< ImageType > itkIterator;
      
  //Reading the input file
  ImageType::Pointer inputImage = btk::ImageHelper<ImageType>::ReadImage(input_file);
  ImageType::Pointer maskImage  = btk::ImageHelper<ImageType>::ReadOrCreateImage(mask_file, inputImage, 1);
  ImageType::Pointer outputImage= btk::ImageHelper<ImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage);
  outputImage->FillBuffer(0.0);
  
  
  float currentMin = 32767;
  float currentMax = -32768;

  itkIterator maskImageIt( maskImage, maskImage->GetLargestPossibleRegion());
  itkIterator inputImageIt( inputImage, inputImage->GetLargestPossibleRegion());    
  itkIterator outputImageIt( outputImage, outputImage->GetLargestPossibleRegion());    
    
  for(inputImageIt.GoToBegin(), maskImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt)
  {
    if(maskImageIt.Get() > 0)
    {
      if(currentMin > inputImageIt.Get() ) currentMin = inputImageIt.Get();
      if(currentMax < inputImageIt.Get() ) currentMax = inputImageIt.Get();        
    }  
  }
  std::cout<<"Image values (inside the mask) max: "<<currentMax<<", min: "<<currentMin<<std::endl;    
    
    

  //rescale coefficients xnew = a*xold + b
  float a = (maxValue-minValue) / (currentMax-currentMin);
  float b = minValue - a*currentMin;
  std::cout<<"linear coefficient for rescaling (xnew = a * xold + b). a="<<a<<", b="<<b<<std::endl;
    
  for(inputImageIt.GoToBegin(), maskImageIt.GoToBegin(), outputImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt, ++outputImageIt)
    if(maskImageIt.Get() > 0)
      outputImageIt.Set( a*inputImageIt.Get() + b );

  btk::ImageHelper<ImageType>::WriteImage(outputImage, output_file);

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}




