/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/08/2012
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
#include <vector>
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"


/*Btk includes*/
#include "btkImageHelper.h"
//#include "btkHistogram.h" 
#include "btkMidwayImageEqualization.h" 


int main (int argc, char* argv[])
{
  
  //Principle of the algorithm: apply a midway histogram equalization to all input images 
  //See these papers: 
  //Cox et al. Dynamic histogram warping of image pairs for constant image brightness, ICIP 1995
  //Delon, Midway Image Equalization, Journal of Mathematical Imaging and Vision (JMIV), vol.21, no.2, pp.119-134, 2004.
  
  //Algorithm:
  //Rescale all images between 0 and 1
  //Compute the histograms of each image
  //Compute the cumulative distribution functions (cdf) H_i
  //Compute the inverse cdf
  //Average all the inverse cdf 
  //Inverse this mean inverse cdf (called H_new)
  //Inverse H_new
  //Apply to all images H_new^-1 o H_i
  
  
  try {  

  TCLAP::CmdLine cmd("It performs midway histogram equalization for a set of 3D images.", ' ', "1.0", true);
  TCLAP::MultiArg<std::string> inputImageArg("i","input_file","input anatomical image files (possible multiple inputs) ",true,"string", cmd);
  TCLAP::MultiArg<std::string> maskImageArg("m","mask_file","input mask image files (possible multiple inputs) ",false,"string", cmd);
  TCLAP::MultiArg<std::string> outputImageArg("o","output_file","output anatomical image files (possible multiple inputs) ",true,"string", cmd);
  TCLAP::ValueArg< unsigned int > binArg("b","bin","number of bins for image intensity (default is 10000)",false,10000,"unsigned int",cmd);
  TCLAP::ValueArg< unsigned int > samArg("s","sample","number of bins for number of samples (voxels) (default is 10000)",false,10000,"unsigned int",cmd);
  TCLAP::MultiArg<std::string> refImageArg("r","ref_file","reference anatomical image files (possible multiple inputs) ",false,"string", cmd);
  TCLAP::MultiArg<std::string> maskRefImageArg("","mask_ref_file","input mask reference image files (possible multiple inputs) ",false,"string", cmd);


  // Parse the args.
  cmd.parse( argc, argv );
  
  std::vector<std::string> input_file    = inputImageArg.getValue();
  std::vector<std::string> mask_file     = maskImageArg.getValue();
  std::vector<std::string> output_file   = outputImageArg.getValue();
  unsigned int numberOfBins              = binArg.getValue();
  unsigned int sampleQuantification      = samArg.getValue();
  std::vector<std::string> ref_file      = refImageArg.getValue();
  std::vector<std::string> mask_ref_file = maskRefImageArg.getValue();
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef ImageType::Pointer                    ImagePointer;
  typedef itk::ImageDuplicator< ImageType >     itkDuplicator;
  typedef itk::ImageRegionIterator< ImageType > itkIterator;
  typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > rescaleFilterType;

  //Initialize image vectors (input, mask, output)
  std::vector<ImagePointer> inputImages;
  inputImages.resize(input_file.size());

  std::vector<ImagePointer> maskImages;
  maskImages.resize(mask_file.size());
  
  std::vector<ImagePointer> outputImages;
  outputImages.resize(input_file.size());
    
  std::vector<ImagePointer> refImages;
  refImages.resize(ref_file.size());

  std::vector<ImagePointer> maskRefImages;
  maskRefImages.resize(mask_ref_file.size());
  
  std::cout<<"Read input images \n";
  inputImages = btk::ImageHelper<ImageType>::ReadImage(input_file);
  
  if(ref_file.size() > 0)
  {
    std::cout<<"Read reference images \n";
    refImages = btk::ImageHelper<ImageType>::ReadImage(ref_file);  
  }
  
  if(mask_file.size()>0)
  {
    std::cout<<"Read mask images \n";
    maskImages = btk::ImageHelper<ImageType>::ReadImage(mask_file);
  }
  else
  {
    std::cout<<"No mask images provided: by default, mask images correspond to the entire image domain\n";    
    maskImages.resize(input_file.size());
    for(unsigned int i=0; i<input_file.size(); i++)
    {
      itkDuplicator::Pointer duplicator = itkDuplicator::New();
      duplicator->SetInputImage( inputImages[i] );
      duplicator->Update();
      maskImages[i] = duplicator->GetOutput();
      maskImages[i]->FillBuffer(1);  
    }
  }
  
  if(mask_ref_file.size()>0)
  {
    std::cout<<"Read mask reference images \n";
    maskRefImages = btk::ImageHelper<ImageType>::ReadImage(mask_ref_file);
  }
  else if(ref_file.size() > 0)
  {
    std::cout<<"No mask reference images provided: by default, mask images correspond to the entire image domain\n";    
    maskRefImages.resize(ref_file.size());
    
    for(unsigned int i=0; i<ref_file.size(); i++)
    {
        maskRefImages[i] = btk::ImageHelper<ImageType>::CreateNewImageFromPhysicalSpaceOf(refImages[i].GetPointer());
      maskRefImages[i]->FillBuffer(1);  
    }
  }


  //Duplicate input images to initialize output images
  for(unsigned int i=0; i<input_file.size(); i++)
  {
    itkDuplicator::Pointer duplicator = itkDuplicator::New();
    duplicator->SetInputImage( inputImages[i] );
    duplicator->Update();
    outputImages[i] = duplicator->GetOutput();
    outputImages[i]->FillBuffer(0);        
  }

  btk::MidwayImageEqualization<PixelType> myTool;
  myTool.SetNumberOfBins(numberOfBins);
  myTool.SetSampleQuantification(sampleQuantification);
  
  std::cout<<"histogram parameters"<<std::endl;
  std::cout<<"number of bins : "<<numberOfBins<<std::endl;
  std::cout<<"sample quantification : "<<sampleQuantification<<std::endl;
  
  if(ref_file.size() > 0)
    myTool.DoWithReference(inputImages, refImages, maskImages, maskRefImages, outputImages);
  else
    myTool.Do(inputImages, maskImages, outputImages);
  
  btk::ImageHelper<ImageType>::WriteImage(outputImages, output_file);
    
  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return 1;
}
