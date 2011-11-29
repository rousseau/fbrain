/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 28/11/2011
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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImage.h"
#include <tclap/CmdLine.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char *argv[] )
{

  try {

  TCLAP::CmdLine cmd("Compare several reconstruction methods: creates a label image explaining which reconstructed image is the best at each voxel", ' ', "Unversioned");

  TCLAP::MultiArg<std::string> inputArg("i","input","image file (possibly multiple inputs)",true,"string",cmd);
  TCLAP::ValueArg<std::string> outArg("o","output","label image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> refArg("r","reference","Reference image",true,"","string",cmd);

  // Parse the args.
  cmd.parse( argc, argv );
   
  // Get the value parsed by each arg. 
  std::vector<std::string> input_file = inputArg.getValue();
  std::string reference_file          = refArg.getValue();
  std::string output_file             = outArg.getValue();
 
  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;
  typedef itk::ImageFileReader< ImageType  >  ImageReaderType;

  unsigned int numberOfImages = input_file.size();

  std::vector< ImagePointer >	inputImages;

  std::cout<<"Reading input (reconstruction) images:"<<std::endl;
  for (unsigned int i = 0; i<numberOfImages; i++)
  {
    std::cout<<"reading : "<<input_file[i]<<std::endl;
    ImageReaderType::Pointer  imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input_file[i] );
    imageReader -> Update();
    inputImages.push_back( imageReader -> GetOutput() );
  }

  std::cout<<"Reading the reference image:"<<reference_file<<std::endl;
  ImageReaderType::Pointer refImageReader = ImageReaderType::New();
  refImageReader -> SetFileName( reference_file );
  refImageReader->Update();
  ImageType::Pointer referenceImage = refImageReader -> GetOutput();

  std::cout<<"Compute the output label image\n";
  typedef itk::ImageRegionIteratorWithIndex< ImageType > itkIterator;
    
  ImageType::Pointer outputImage = ImageType::New();    
  outputImage->SetRegions(referenceImage->GetLargestPossibleRegion());
  outputImage->SetSpacing(referenceImage->GetSpacing() );
  outputImage->SetOrigin(referenceImage->GetOrigin() );
  outputImage->SetDirection(referenceImage->GetDirection() );
  outputImage->Allocate();
  outputImage->FillBuffer(0);    

  std::vector<double> counter(numberOfImages+1);
  for(uint i = 0; i<numberOfImages+1; i++)
    counter[i] = 0;
  double numberOfVoxels = 0;
    
  itkIterator itref( referenceImage, referenceImage->GetRequestedRegion() );
  ImageType::IndexType index;
  for(itref.GoToBegin(); !itref.IsAtEnd(); ++itref){
    index = itref.GetIndex();
    PixelType outputValue = 1; //best image is set to one 
    PixelType minError = fabs(itref.Get() - inputImages[0]->GetPixel(index));  
    
    for (unsigned int i = 1; i<numberOfImages; i++){
      PixelType error = fabs(itref.Get() - inputImages[i]->GetPixel(index));
      if(error == minError){
        outputValue = 0;
      }
      if(error < minError){
        minError = error;
        outputValue = i+1;
      }
    }
    counter[outputValue]+=1;
    numberOfVoxels+=1;
    outputImage->SetPixel(index,outputValue);
  }
  
  for(uint i = 0; i<numberOfImages+1; i++){
    std::cout<<"percentage of voxels with label "<<i<<": "<<counter[i]/numberOfVoxels*100<<"\n";
  }
    
  
  std::cout<<"Writing output image: "<<output_file<<std::endl;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer writer =  WriterType::New();
  writer->SetFileName( output_file );
  writer->SetInput( outputImage );
  writer->Update();


  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

