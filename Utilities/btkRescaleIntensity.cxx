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



#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"

#include <string>
#include <iomanip>

#include <tclap/CmdLine.h>


int main(int argc, char *argv[])
{
  
  try{

    TCLAP::CmdLine cmd("Rescale the intensity values of an image (basic ITK filter) using short values", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputImageArg("i","input_file","input image file (short)",true,"","string");
    cmd.add( inputImageArg );
    TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string");
    cmd.add( outputImageArg );
    TCLAP::ValueArg< short > minArg("","min","minimum value for the rescaling (default is 0)",false,0,"short");
    cmd.add( minArg );
    TCLAP::ValueArg< short > maxArg("","max","maximum value for the rescaling (default is 255)",false,255,"short");
    cmd.add( maxArg );  
    
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::string input_file       = inputImageArg.getValue();
    std::string output_file      = outputImageArg.getValue();    
    
    short minimum              = minArg.getValue();
    short maximum              = maxArg.getValue();
    
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
      
  //Reading the input file
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( input_file  );
  reader->Update();

  ImageType::Pointer input_image = reader->GetOutput();
  
  itk::MinimumMaximumImageCalculator< ImageType >::Pointer minmax = itk::MinimumMaximumImageCalculator< ImageType >::New();
  minmax->SetImage(input_image);
  minmax->Compute();
  std::cout<<"Dynamique - input image :"<<minmax->GetMaximum()<<" "<<minmax->GetMinimum()<<"\n";

  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetOutputMinimum(minimum);    
  filter->SetOutputMaximum(maximum);
  filter->SetInput(input_image);  
  filter->UpdateLargestPossibleRegion();
  
  ImageType::Pointer output_image = filter->GetOutput();  
    
  //Writing the output file
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( output_file  );
  writer -> SetInput( output_image );
  writer -> Update();

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}




