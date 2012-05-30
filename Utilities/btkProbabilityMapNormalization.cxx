/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 06/03/2012
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
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <string>
#include <iomanip>

#include <tclap/CmdLine.h>


int main(int argc, char *argv[])
{
  
  try{

    TCLAP::CmdLine cmd("Normalize in a voxelwise manner a set of (positive) probability maps", ' ', "1.0", true);
    
    TCLAP::MultiArg<std::string> inputImageArg("i","input_image","input images",true,"string",cmd);
    TCLAP::MultiArg<std::string> outputImageArg("o","output_image","output images",true,"string",cmd);
    TCLAP::ValueArg<float> thresholdArg  ("t","threshold","Threshold on the sum (for each voxel) to avoid calculation error) (default = 0.00001)",false, 0.00001,"float",cmd);
    
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::vector<std::string> input_file  = inputImageArg.getValue();
    std::vector<std::string> output_file = outputImageArg.getValue();
    float threshold = thresholdArg.getValue();
    
  
    //ITK declaration
    typedef float PixelType;
    const   unsigned int        Dimension = 3;
    
    typedef itk::Image< PixelType, Dimension >    ImageType;
    typedef itk::ImageFileReader< ImageType >     ReaderType;
    typedef itk::ImageFileWriter< ImageType >     WriterType;
    typedef ImageType::Pointer itkPointer;
    typedef itk::ImageRegionIteratorWithIndex< ImageType > itkIteratorWithIndex;
        
    //Set image pointers
    std::vector<itkPointer>     inputImages;
    std::vector<itkPointer>     outputImages;
    
    inputImages.resize(input_file.size());
    outputImages.resize(output_file.size());
    
    std::cout<<"Reading the input files\n";
    for(int i=0; i<input_file.size(); i++){
      std::cout<<input_file[i]<<"\n";
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( input_file[i]  );
      reader->Update();
      inputImages[i] = reader->GetOutput();     
    }
    
    std::cout<<"Initialize the output images\n";
    for(int i=0; i<output_file.size(); i++)
    {
      outputImages[i] = ImageType::New();
      outputImages[i]->SetRegions( inputImages[i]->GetLargestPossibleRegion() );
      outputImages[i]->SetSpacing( inputImages[i]->GetSpacing() );
      outputImages[i]->SetOrigin( inputImages[i]->GetOrigin() );
      outputImages[i]->SetDirection( inputImages[i]->GetDirection() );
      outputImages[i]->Allocate();
      outputImages[i]->FillBuffer(0);       
    }
    
    
    itkIteratorWithIndex it( inputImages[0], inputImages[0]->GetLargestPossibleRegion());
    
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it){
      double sum = 0;
      ImageType::IndexType index = it.GetIndex();
      for(int i=0; i<input_file.size(); i++){
        if(inputImages[i]->GetPixel(index) > 0)
          sum += inputImages[i]->GetPixel(index);  
      }  
      
      if(sum > threshold)
        for(int i=0; i<output_file.size(); i++)
        {
          if(inputImages[i]->GetPixel(index) > 0)
            outputImages[i]->SetPixel(index,inputImages[i]->GetPixel(index) / sum); 
        }
    }

    std::cout<<"Writing the output images\n";
    for(int i=0; i<output_file.size(); i++)
    {
      std::cout<<output_file[i]<<"\n";
      WriterType::Pointer writer = WriterType::New();  
      writer->SetFileName( output_file[i] );
      writer->SetInput( outputImages[i] );
      writer->Update();      
    }
    
 
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}




