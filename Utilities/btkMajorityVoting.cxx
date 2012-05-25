/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/09/2010
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


/*  Standard includes */
#include "vector"
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"


int main( int argc, char *argv[] )
{
  try {
    
    TCLAP::CmdLine cmd("Compute label map using majority voting rule", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg ("i","input","Label image file.", true,"string",cmd);
    TCLAP::ValueArg<std::string> weightArg("w","weight","Weight image file. ",false,"","string",cmd);
    TCLAP::ValueArg<std::string> outArg   ("o","output","Output file of the estimated label image.", true,"","string",cmd);
    
    // Parse the argv array.
    cmd.parse( argc, argv );
    
    // Get the value parsed by each arg. 
    std::vector<std::string> label_file       = inputArg.getValue();
    std::string output_file      = outArg.getValue();
    std::string weight_file      = weightArg.getValue();

  //ITK declaration
  typedef short InputPixelType;
  typedef float OutputPixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< InputPixelType, Dimension >    InputImageType;
  typedef itk::Image< OutputPixelType, Dimension >    OutputImageType;


  typedef InputImageType::Pointer InputImagePointer;
  typedef OutputImageType::Pointer OutputImagePointer;

  typedef itk::ImageFileReader< InputImageType >  ReaderType;
  typedef itk::ImageFileWriter< InputImageType >  WriterType;
  typedef itk::ImageFileWriter< OutputImageType >  OutputWriterType;


  std::vector<InputImagePointer> labelImages;
  labelImages.resize(label_file.size());

  //Reading input data
  for(uint i=0;i<label_file.size();i++){
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( label_file[i]  );
    reader->Update();
    labelImages[i] = reader->GetOutput();
  }
  InputImageType::RegionType region;
  InputImageType::SizeType size;
  InputImageType::SpacingType spacing;

  //we assume that all images have the same size, spacing etc.
  region = labelImages[0]->GetLargestPossibleRegion();
  size = region.GetSize();
  spacing = labelImages[0]->GetSpacing();


  //Allocate output images
  InputImagePointer outputImage = InputImageType::New();
  outputImage->SetRegions( region );
  outputImage->SetSpacing( spacing );
  outputImage->SetOrigin( labelImages[0]->GetOrigin() );
  outputImage->SetDirection( labelImages[0]->GetDirection() );
  outputImage->Allocate();
  outputImage->FillBuffer(0); 

  OutputImagePointer weightImage = OutputImageType::New();
  weightImage->SetRegions( region );
  weightImage->SetSpacing( spacing );
  weightImage->SetOrigin( labelImages[0]->GetOrigin() );
  weightImage->SetDirection( labelImages[0]->GetDirection() );
  weightImage->Allocate();
  weightImage->FillBuffer(0); 


  int x,y,z;
  uint n = labelImages.size();

  #pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)size[2]; z++)
    for(y=0; y < (int)size[1]; y++)
      for(x=0; x < (int)size[0]; x++){

        InputImageType::IndexType p;
        p[0] = x;
        p[1] = y;
        p[2] = z;

        InputPixelType label = 0;
        InputPixelType outputLabel = 0;
        std::map<InputPixelType, float> map;
        float weight = 1.0;

        for(uint i=0; i < n; i++){
            label = labelImages[i]->GetPixel(p);
            map[label] += weight;
        }
        std::map<InputPixelType, float>::iterator mapIt;
        float wmax = 0;
        for(mapIt = map.begin (); mapIt != map.end (); ++mapIt){
          if( (*mapIt).second > wmax ){
            wmax = (*mapIt).second;
            outputLabel = (*mapIt).first;
          }
        }
 
        outputImage->SetPixel( p, outputLabel );
        weightImage->SetPixel( p, wmax );	  	  
	  	  

  }


  //Write the result 
  WriterType::Pointer writer = WriterType::New();  
  writer->SetFileName( output_file );
  writer->SetInput( outputImage );
  writer->Update();  

  OutputWriterType::Pointer outputWriter = OutputWriterType::New();  
  outputWriter->SetFileName( weight_file );
  outputWriter->SetInput( weightImage );
  outputWriter->Update();  



  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}
