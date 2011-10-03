/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 27/09/2011
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


#include "itkImageFileReader.h"

#include "itkImage.h"
#include "itkMinimumMaximumImageCalculator.h"

#include <string>
#include <iomanip>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <tclap/CmdLine.h>


int main(int argc, char *argv[])
{
  
  try{

    TCLAP::CmdLine cmd("Compute PSNR between the input image and a reference image.", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
    cmd.add( inputImageArg );
    TCLAP::ValueArg<std::string> referenceImageArg("r","reference_file","reference image file (short)",true,"","string");
    cmd.add( referenceImageArg );
    TCLAP::ValueArg< int > borderArg("b","border","border size not considered during the computation (default is 0)",false,0,"int");
    cmd.add( borderArg );
    TCLAP::ValueArg< float > thresholdArg("t","threshold","threshold acting as a padding value (default is -1)",false,-1,"float");
    cmd.add( thresholdArg );
    TCLAP::ValueArg<std::string> inputMaskArg("m","mask_file","filename of the mask image",false,"","string");
    cmd.add( inputMaskArg );
    TCLAP::ValueArg< int > fuzzyArg("f","fuzzy","fuzzy (1) or hard (0) psnr (using the mask) (default is 0)",false,0,"int");
    cmd.add( fuzzyArg );    
    
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::string input_file       = inputImageArg.getValue();
    std::string reference_file   = referenceImageArg.getValue();    
    std::string mask_file        = inputMaskArg.getValue();
    
    float threshold              = thresholdArg.getValue();
    int border                   = borderArg.getValue();
    int fuzzy                    = fuzzyArg.getValue();
    
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  //Reading the 2 input files
  ReaderType::Pointer reference_reader = ReaderType::New();
  reference_reader->SetFileName( reference_file  );
  reference_reader->Update();

  ReaderType::Pointer input_reader = ReaderType::New();
  input_reader->SetFileName( input_file  );
  input_reader->Update();

  ImageType::Pointer reference_image = reference_reader->GetOutput();
  ImageType::Pointer input_image = input_reader->GetOutput();

  ImageType::RegionType region;
  region.SetSize(reference_image->GetLargestPossibleRegion().GetSize());


  ImageType::Pointer mask_image;
  if(mask_file !="")
  {
    ReaderType::Pointer mask_reader = ReaderType::New();
    mask_reader->SetFileName( mask_file  );
    mask_reader->Update();
    mask_image = mask_reader->GetOutput();
  }
  else
  {
    mask_image = ImageType::New();
    mask_image->SetRegions(region);
    mask_image->Allocate();
    mask_image->FillBuffer(1);
  }  
  
  itk::MinimumMaximumImageCalculator< ImageType >::Pointer minmax = itk::MinimumMaximumImageCalculator< ImageType >::New();
  minmax->SetImage(reference_image);
  minmax->Compute();
  std::cout<<"Dynamique - reference image :"<<minmax->GetMaximum()<<" "<<minmax->GetMinimum()<<"\n";
  if(threshold == -1 ) 
    threshold = minmax->GetMinimum();
  float d = minmax->GetMaximum() - minmax->GetMinimum();

  minmax->SetImage(input_image);
  minmax->Compute();
  std::cout<<"Dynamique - input image :"<<minmax->GetMaximum()<<" "<<minmax->GetMinimum()<<"\n";

  minmax->SetImage(mask_image);
  minmax->Compute();
  std::cout<<"Dynamique - mask image :"<<minmax->GetMaximum()<<" "<<minmax->GetMinimum()<<"\n";
  PixelType max_mask = minmax->GetMaximum();

  double res = 0.0;
  double tmp = 0.0;
  double nbp = 0.0;
  ImageType::SizeType imageSize;
  imageSize = region.GetSize();  
  ImageType::IndexType pixelIndex;	

  int bip = border;
  if(fuzzy==0){
    for(int k=bip;k<(int)(imageSize[2]-bip);k++)
      for(int j=bip;j<(int)(imageSize[1]-bip);j++)
        for(int i=bip;i<(int)(imageSize[0]-bip);i++)
        {
	        pixelIndex[0] = i;
	        pixelIndex[1] = j;
	        pixelIndex[2] = k;
        if( (reference_image->GetPixel(pixelIndex) >= threshold) && (mask_image->GetPixel(pixelIndex) > 0) ){
          tmp = reference_image->GetPixel(pixelIndex) - input_image->GetPixel(pixelIndex);
          res += tmp*tmp;
          nbp++;
        }
        }
  }
  if(fuzzy==1){
    for(int k=bip;k<(int)(imageSize[2]-bip);k++)
      for(int j=bip;j<(int)(imageSize[1]-bip);j++)
        for(int i=bip;i<(int)(imageSize[0]-bip);i++)
        {
          pixelIndex[0] = i;
	        pixelIndex[1] = j;
	        pixelIndex[2] = k;
        if( (reference_image->GetPixel(pixelIndex) >= threshold) && (mask_image->GetPixel(pixelIndex) > 0) ){
          tmp = reference_image->GetPixel(pixelIndex) - input_image->GetPixel(pixelIndex);
          float weight = mask_image->GetPixel(pixelIndex) / max_mask;
          res += tmp*tmp * weight;
          nbp += weight;
        }
        }
  }

  std::cout<<"Number of points used to compute PSNR : "<<nbp<<"\n";
  std::cout<<"MSE : "<<res/nbp<<"\n";
  res = 10 * log10(d*d*nbp/res);
  std::cout<<"Range : "<<d<<"\n";
  std::cout<<"PSNR : "<<res<<"\n";


  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}




