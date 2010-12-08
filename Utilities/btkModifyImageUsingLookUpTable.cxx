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
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImage.h"

#include <iostream>
#include <fstream>

#include <string>
#include <iomanip>

#include <map>



int main(int ac,char *av[])
{
  std::string image_file = "";
  std::string table_file = "";
  std::string output_file = "";

  //Parsing the command line using Boost
  po::options_description desc("Usage:");
  desc.add_options()
    ("help", "produce help message")
    ("image-file,i", po::value< std::string >(&image_file), "input image file (short)")
    ("table-file,t", po::value< std::string >(&table_file), "input look up table file (ascii)")
    ("output-file,o", po::value< std::string >(&output_file), "output file")
    ;
  
  po::variables_map vm;        
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
    std::cout << "This program modifies the values of the input image using a look up table.\n";
    std::cout << desc << "\n";
    return 1;
  }
  if (vm.count("image-file")) {
    std::cout<<"   input image file "<<image_file<<"\n";
  }
  else{
    std::cout<<"No input image file given. Exit.\n";
    return 1;
  }  
  if (vm.count("table-file")) {
    std::cout<<"   input look-up table file "<<table_file<<"\n";
  }
  else{
    std::cout<<"No table file given. Exit.\n";
    return 1;
  }  
  //end of parsing command line

  

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
}




