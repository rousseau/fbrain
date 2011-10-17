/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 05/10/2011
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
#include "itkImageRegionIterator.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <string>
#include <iomanip>

#include <tclap/CmdLine.h>


int main(int argc, char *argv[])
{
  
  try{

    TCLAP::CmdLine cmd("Compute a weighted mean using 2 input images and one weight image (or scalar): output = w * I1 + (1-w) * I2.", ' ', "1.0", true);
    
    TCLAP::MultiArg<std::string> inputImageArg("i","input_image","input image (short) (possible multiple inputs) ",true,"string");
    cmd.add( inputImageArg );
    TCLAP::ValueArg<std::string> weightImageArg("w","weight_image","weight image file (float)",false,"","string");
    cmd.add( weightImageArg );
    TCLAP::ValueArg< float > scalarArg("s","scalar_weight","weight as a scalar if no weight image is used (default is 0.5)",false,0.5,"int");
    cmd.add( scalarArg );
    TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string");
    cmd.add( outputImageArg ); 
    
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::vector<std::string> input_file = inputImageArg.getValue();
    std::string weight_file      = weightImageArg.getValue();
    std::string output_file      = outputImageArg.getValue();
        
    float scalar                 = scalarArg.getValue();
    
  
    //ITK declaration
    typedef short PixelType;
    const   unsigned int        Dimension = 3;
    
    typedef itk::Image< PixelType, Dimension >    ImageType;
    typedef itk::ImageFileReader< ImageType >     ReaderType;
    typedef itk::ImageFileWriter< ImageType >     WriterType;
    typedef ImageType::Pointer ImagePointer;
    
    typedef itk::Image< float, Dimension >             FloatImageType;
    typedef itk::ImageFileReader< FloatImageType >     FloatReaderType;
    typedef FloatImageType::Pointer FloatImagePointer;
    
    //Set image pointers
    ImagePointer image1;
    ImagePointer image2;
    FloatImagePointer weightImage;
    
    //Reading the 2 input files
    ReaderType::Pointer input_reader1 = ReaderType::New();
    input_reader1->SetFileName( input_file[0]  );
    input_reader1->Update();
    image1 = input_reader1->GetOutput();
    
    ReaderType::Pointer input_reader2 = ReaderType::New();
    input_reader2->SetFileName( input_file[1]  );
    input_reader2->Update();
    image2 = input_reader2->GetOutput();
    
    //Reading the weight image or create one using the scalar value
    if(weight_file != ""){
      FloatReaderType::Pointer weight_reader = FloatReaderType::New();
      weight_reader->SetFileName( weight_file  );
      weight_reader->Update();
      weightImage = weight_reader->GetOutput();
    }
    else {
      weightImage->SetRegions( image1->GetLargestPossibleRegion() );
      weightImage->Allocate();
      weightImage->FillBuffer(scalar);
    }

    //Set output image pointer
    ImagePointer outputImage = ImageType::New();
    outputImage->SetRegions( image1->GetLargestPossibleRegion() );
    outputImage->SetSpacing( image1->GetSpacing() );
    outputImage->SetOrigin( image1->GetOrigin() );
    outputImage->SetDirection( image1->GetDirection() );
    outputImage->Allocate();
    outputImage->FillBuffer(0); 
    
    //Iterator declarations
    typedef itk::ImageRegionIterator< ImageType > shortIterator;
    typedef itk::ImageRegionIterator< FloatImageType > floatIterator;
    
    //Compute the weighted mean using the weight image
    //output = w * I1 + (1-w) * I2
    
    shortIterator inputIt1(image1, image1->GetLargestPossibleRegion());
    shortIterator inputIt2(image2, image2->GetLargestPossibleRegion());
    floatIterator weightIt(weightImage, weightImage->GetLargestPossibleRegion());
    shortIterator outputIt(outputImage,outputImage->GetLargestPossibleRegion());
    
    for ( inputIt1.GoToBegin(), inputIt2.GoToBegin(), weightIt.GoToBegin(), outputIt.GoToBegin(); !inputIt1.IsAtEnd(); ++inputIt1, ++inputIt2, ++weightIt, ++outputIt)
      outputIt.Set( weightIt.Get() * inputIt1.Get() + (1-weightIt.Get()) * inputIt2.Get());


    //Write the result image
    WriterType::Pointer writer = WriterType::New();  
    writer->SetFileName( output_file );
    writer->SetInput( outputImage );
    writer->Update();  
    
 
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}




