/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

6 december 2010
rousseau@unistra.fr

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
*/



/* Standard includes */
#include <tclap/CmdLine.h>
#include "iostream"
#include "fstream"
#include "string"
#include "iomanip"
#include "map"

/* Itk includes */
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"


int main(int argc, char** argv)
{

try {  

	TCLAP::CmdLine cmd("Command description message", ' ', "1.0", true);

	TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
	cmd.add( inputImageArg );
	TCLAP::ValueArg<std::string> inputTableArg("t","table_file","input look up table file (ascii)",true,"","string");
	cmd.add( inputTableArg );
	TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string");
	cmd.add( outputImageArg );

	// Parse the args.
	cmd.parse( argc, argv );

	// Get the value parsed by each arg. 
  std::string image_file = inputImageArg.getValue();
  std::string table_file = inputTableArg.getValue();
  std::string output_file = outputImageArg.getValue();

  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  std::cout<<"Reading the input file\n";
  ReaderType::Pointer image_reader = ReaderType::New();
  image_reader->SetFileName( image_file  );
  image_reader->Update();

  ImageType::Pointer input_image = image_reader->GetOutput();
  ImageType::Pointer output_image = image_reader->GetOutput();  //by default, the output image is set to the input image

  std::ifstream lookupTableFile(table_file.c_str(), std::ios::in);  // on ouvre le fichier en lecture
  std::map<PixelType, PixelType> lookupTable;

  std::cout<<"Loading the look-up table\n";
  PixelType index = 0;
  PixelType value = 0;
  while( !lookupTableFile.eof() ){
    lookupTableFile >> index;
    lookupTableFile >> value;
    lookupTable[index] = value;
  }

  typedef itk::ImageRegionIterator< ImageType > Iterator;
  Iterator itData( input_image, input_image->GetRequestedRegion() );
  Iterator itOutput( output_image, output_image->GetRequestedRegion() );

  for(itData.GoToBegin(), itOutput.GoToBegin(); !itData.IsAtEnd(); ++itData, ++itOutput)
    itOutput.Set( lookupTable.find( itData.Get() )->second  );
    
  lookupTableFile.close();
  std::cout<<"Writing the output image\n";
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( output_file ); 
  writer->SetInput( output_image ); 
  writer->Update();

  return 1;

	} catch (TCLAP::ArgException &e)  // catch any exceptions
	{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}




