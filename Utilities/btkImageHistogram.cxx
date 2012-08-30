/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 30/08/2012
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
#include "btkHistogram.h" 


int main(int argc, char *argv[])
{
  
  try{

    TCLAP::CmdLine cmd("Analysis of images through histograms", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputImageArg("i","input_file","input image file (short)",true,"","string", cmd);
    TCLAP::ValueArg<std::string> maskImageArg("m","mask_file","mask image file (short)",true,"","string", cmd);
    TCLAP::ValueArg< unsigned int > binArg("b","bin","number of bins for image intensity (default is 10000)",false,10000,"unsigned int",cmd);
    TCLAP::ValueArg< unsigned int > samArg("s","sample","number of bins for number of samples (voxels) (default is 10000)",false,10000,"unsigned int",cmd);
    TCLAP::ValueArg<std::string> outputHistogramArg("","output_histogram","file name of the image histogram",false,"","string", cmd);
    TCLAP::ValueArg<std::string> outputNormalizedHistogramArg("","output_normalized_histogram","file name of the normalized image histogram",false,"","string", cmd);
    TCLAP::ValueArg<std::string> outputCDFArg("","output_cdf","file name of the cumulative density function",false,"","string", cmd);
    TCLAP::ValueArg<std::string> outputNormalizedCDFArg("","output_normalized_cdf","file name of the normalized cumulative density function",false,"","string", cmd);
    TCLAP::ValueArg<std::string> outputICDFArg("","output_icdf","file name of the inverse cumulative density function",false,"","string", cmd);
    TCLAP::ValueArg<std::string> outputNormalizedICDFArg("","output_normalized_icdf","file name of the normalized inverse cumulative density function",false,"","string", cmd);

    
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg.
    std::string input_file       = inputImageArg.getValue();
    std::string mask_file        = maskImageArg.getValue();
    unsigned int numberOfBins         = binArg.getValue();
    unsigned int sampleQuantification = samArg.getValue();

  std::string output_histogram_file            = outputHistogramArg.getValue();       
  std::string output_normalized_histogram_file = outputNormalizedHistogramArg.getValue();       
  std::string output_cdf_file                  = outputCDFArg.getValue();       
  std::string output_normalized_cdf_file       = outputNormalizedCDFArg.getValue();       
  std::string output_icdf_file                 = outputICDFArg.getValue();       
  std::string output_normalized_icdf_file      = outputNormalizedICDFArg.getValue();       

  
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
    
  for(inputImageIt.GoToBegin(), maskImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt)
  {
    if(maskImageIt.Get() > 0)
    {
      if(currentMin > inputImageIt.Get() ) currentMin = inputImageIt.Get();
      if(currentMax < inputImageIt.Get() ) currentMax = inputImageIt.Get();        
    }  
  }
  std::cout<<"Image values (inside the mask) max: "<<currentMax<<", min: "<<currentMin<<std::endl;    

  btk::Histogram histogram;

  histogram.SetNumberOfBins(numberOfBins);
  histogram.SetSampleQuantification(sampleQuantification);
  histogram.SetLowerBound(currentMin);
  histogram.SetUpperBound(currentMax);    
  histogram.Setup();

  for(inputImageIt.GoToBegin(), maskImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt)
  {
    if(maskImageIt.Get() > 0)
    {
      histogram.AddSample( inputImageIt.Get() );        
    }  
  }

  //COMPUTE THINGS ------------------------------------------------------------------------------------------------------
  histogram.NormalizeData();
  histogram.ComputeCumulativeDistributionFunction();
  histogram.ComputeInverseCumulativeDistributionFunction();
  histogram.ComputeNormalizedCumulativeDistributionFunction();
  histogram.ComputeNormalizedInverseCumulativeDistributionFunction();      

  // SAVE OUTPUT --------------------------------------------------------------------------------------------------------
  if(output_histogram_file != "")
    histogram.SaveHistogram(output_histogram_file);  
  
  if(output_normalized_histogram_file != "")
    histogram.SaveNormalizedHistogram(output_normalized_histogram_file);  
  
  if(output_cdf_file != "")
    histogram.SaveCumulativeDistributionFunction(output_cdf_file);  
  
  if(output_normalized_cdf_file != "")
    histogram.SaveNormalizedCumulativeDistributionFunction(output_normalized_cdf_file);  
  
  if(output_icdf_file != "")
    histogram.SaveInverseCumulativeDistributionFunction(output_icdf_file);  
  
  if(output_normalized_icdf_file != "")
    histogram.SaveNormalizedInverseCumulativeDistributionFunction(output_normalized_icdf_file);  
  


  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}




