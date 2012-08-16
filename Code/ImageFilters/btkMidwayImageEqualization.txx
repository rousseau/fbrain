#ifndef __btkMidwayImageEqualization_TXX__
#define __btkMidwayImageEqualization_TXX__

#include "btkMidwayImageEqualization.h"

namespace btk
{
  
template <typename T>
void MidwayImageEqualization<T>::SetNumberOfBins(unsigned int n)
{
  m_numberOfBins = n;
}

template <typename T>
void MidwayImageEqualization<T>::SetSampleQuantification(unsigned int n)
{
  m_sampleQuantification = n;
}    

template <typename T>
void MidwayImageEqualization<T>::Do(std::vector<itkTPointer> inputImages, std::vector<itkTPointer> maskImages, std::vector<itkTPointer> outputImages)
{
  std::vector<itkTPointer> rescaledImages;
  rescaledImages.resize(inputImages.size());  
  
  //Duplicate input images to initialize rescaled images
  for(unsigned int i=0; i<inputImages.size(); i++)
  {
    typename itkTDuplicator::Pointer duplicator = itkTDuplicator::New();
    duplicator->SetInputImage( inputImages[i] );
    duplicator->Update();
    rescaledImages[i] = duplicator->GetOutput();
    rescaledImages[i]->FillBuffer(0);        
  }

  //Rescaling between 0 and 32767
  //Don't ITK for rescaling because we need to take into account the mask of the image
  T minValue = 0;
  T maxValue = 32767;

  for(unsigned int i=0; i<inputImages.size(); i++)
  {
    short currentMin = 32767;
    short currentMax = -32768;
    itkTIterator maskImageIt( maskImages[i], maskImages[i]->GetLargestPossibleRegion());
    itkTIterator inputImageIt( inputImages[i], inputImages[i]->GetLargestPossibleRegion());    
    
    for(inputImageIt.GoToBegin(), maskImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt)
    {
      if(maskImageIt.Get() > 0)
      {
        if(currentMin > inputImageIt.Get() ) currentMin = inputImageIt.Get();
        if(currentMax < inputImageIt.Get() ) currentMax = inputImageIt.Get();        
      }  
    }
    std::cout<<"Image "<<i<<" (inside the mask) max: "<<currentMax<<", min: "<<currentMin<<std::endl;    
    
    std::cout<<"Rescale current image between 0 and 1 (32 767 for  short).\n";
    itkTIterator rescaledImageIt( rescaledImages[i], rescaledImages[i]->GetLargestPossibleRegion());    

    //rescale coefficients xnew = a*xold + b
    float a = (maxValue-minValue) / (currentMax-currentMin);
    float b = minValue - a*currentMin;

    for(inputImageIt.GoToBegin(), maskImageIt.GoToBegin(), rescaledImageIt.GoToBegin(); !inputImageIt.IsAtEnd(); ++inputImageIt, ++maskImageIt, ++rescaledImageIt)
      if(maskImageIt.Get() > 0)
        rescaledImageIt.Set( a*inputImageIt.Get() + b );
        
  }
      
  T min = 0;
  T max = 1;
  
  //FROM THE ITK WIKI WEBSITE: The regular (non-masked) ImageToHistogramFilter has some problems, so no sense in trying this until those are resolved (14/08/2012). 
  //Meaning that there are some issues with image to histogram filter
  //-> do it on our own

  std::cout<<"Initialize histogram parameters\n";
  
  std::vector< btk::Histogram > histograms;
  histograms.resize( inputImages.size() );

  for(unsigned int i=0; i<inputImages.size(); i++)
  {
    histograms[i].SetNumberOfBins(m_numberOfBins);
    histograms[i].SetSampleQuantification(m_sampleQuantification);
    histograms[i].SetLowerBound(min);
    histograms[i].SetUpperBound(max);    
    histograms[i].Setup();    
  }
  

  for(unsigned int i=0; i<inputImages.size(); i++)
  {
    std::cout<<"Analyzing Image "<<i+1<<std::endl;
    
    std::cout<<" Filling histograms "<<std::endl; 
    itkTIterator maskImageIt( maskImages[i], maskImages[i]->GetLargestPossibleRegion());
    itkTIterator rescaledImageIt( rescaledImages[i], rescaledImages[i]->GetLargestPossibleRegion());
    
    for(rescaledImageIt.GoToBegin(), maskImageIt.GoToBegin(); !rescaledImageIt.IsAtEnd(); ++rescaledImageIt, ++maskImageIt)
    {
      if(maskImageIt.Get() > 0)
      {
        float inputValue = rescaledImageIt.Get()  / (maxValue*1.0);  
        histograms[i].AddSample( inputValue );
      }
    }
    
    std::cout<<" Normalize the histogram "<<std::endl;
    histograms[i].NormalizeData();
    
    std::cout<<" Compute cumulative distribution functions (cdf) and inverse cdf "<<std::endl;
    histograms[i].ComputeNormalizedCumulativeDistributionFunction();
    histograms[i].ComputeNormalizedInverseCumulativeDistributionFunction();            
  }  
  
  
  std::cout<<"Compute mean normalized inverse cdf \n";
  std::vector<float> meanNormalizedInverseCumulativeDistributionFunction(m_sampleQuantification,0);
  for(unsigned int i=0; i<inputImages.size(); i++)
    for(unsigned int j=0; j<m_sampleQuantification; j++)
      meanNormalizedInverseCumulativeDistributionFunction[j] = histograms[i].m_normalizedInverseCumulativeDistributionFunction[j];
  for(unsigned int j=0; j<m_sampleQuantification; j++)
      meanNormalizedInverseCumulativeDistributionFunction[j] /= inputImages.size();    
   
  
  for(unsigned int i=0; i<inputImages.size(); i++)
  {
    std::cout<<"Applying on Image "<<i+1<<std::endl;
    
    itkTIterator maskImageIt( maskImages[i], maskImages[i]->GetLargestPossibleRegion());
    itkTIterator rescaledImageIt( rescaledImages[i], rescaledImages[i]->GetLargestPossibleRegion());
    itkTIterator outputImageIt( outputImages[i], outputImages[i]->GetLargestPossibleRegion());
    
    for(rescaledImageIt.GoToBegin(), maskImageIt.GoToBegin(), outputImageIt.GoToBegin(); !rescaledImageIt.IsAtEnd(); ++rescaledImageIt, ++maskImageIt, ++outputImageIt)
    {
      if(maskImageIt.Get() > 0)
      {
        float inputValue = rescaledImageIt.Get()  / (maxValue*1.0);
        
        //from pixel value to bin index
        unsigned int inputBin = (int) (inputValue * histograms[i].m_aCoefficient + histograms[i].m_bCoefficient);
        
        //from bin index to cdf of the current image
        float tempValue = histograms[i].m_normalizedCumulativeDistributionFunction[inputBin]; // tempValue in [0,sampleQuantification[
        
        //from cdf to inverse mean cdf
        T outputValue = meanNormalizedInverseCumulativeDistributionFunction[(int)(floor(tempValue))];
        outputImageIt.Set(outputValue);  
                       
      }
    }            
  }     

  
}

}
#endif // btkNLMTool_TXX
