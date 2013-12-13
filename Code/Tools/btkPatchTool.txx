/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 january 2013
< rousseau at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include "btkPatchTool.h"


namespace btk
{

template<typename T1, typename T2>
void PatchTool<T1,T2>::SetSpatialBandwidth(int s, Patch<T1> & patch)
{
  //std::cout<<"Computing spatial bandwidth (taking into account possible image anisotropy)\n";    
  float minVoxSz = patch.GetImageSpacing()[0];
  if(patch.GetImageSpacing()[1] < minVoxSz) minVoxSz = patch.GetImageSpacing()[1];
  if(patch.GetImageSpacing()[2] < minVoxSz) minVoxSz = patch.GetImageSpacing()[2];

  m_HalfSpatialBandwidth[0] = (int)(0.5 + s * minVoxSz / patch.GetImageSpacing()[0]);
  m_HalfSpatialBandwidth[1] = (int)(0.5 + s * minVoxSz / patch.GetImageSpacing()[1]);
  m_HalfSpatialBandwidth[2] = (int)(0.5 + s * minVoxSz / patch.GetImageSpacing()[2]);
  //std::cout<<"half spatialBandwidth : "<<m_HalfSpatialBandwidth[0]<<" "<<m_HalfSpatialBandwidth[1]<<" "<<m_HalfSpatialBandwidth[2]<<"\n";  
  
  m_FullSpatialBandwidth[0] = 2 * m_HalfSpatialBandwidth[0] + 1;
  m_FullSpatialBandwidth[1] = 2 * m_HalfSpatialBandwidth[1] + 1;
  m_FullSpatialBandwidth[2] = 2 * m_HalfSpatialBandwidth[2] + 1;
}

template<typename T1, typename T2>
void PatchTool<T1,T2>::AddPatchToImage(typename itkT1Image::IndexType & p, Patch<T2> & patch, itkT1ImagePointer & image, itkFloatImagePointer & weightImage, double weight)
{
  //this function add a patch value to the denoised image
  typename itkT1Image::RegionType imageRegion;
  typename itkT1Image::RegionType patchRegion;
  patch.ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkT1Iterator    imageIt( image, imageRegion);
  itkFloatIterator weightIt( weightImage, imageRegion);  
  itkT2Iterator    patchIt( patch.GetData(), patchRegion);

  for ( imageIt.GoToBegin(), patchIt.GoToBegin(), weightIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++patchIt, ++weightIt)
  {
    imageIt.Set( imageIt.Get() + patchIt.Get() );
    weightIt.Set( weightIt.Get() + weight );
  }

}

template<typename T1, typename T2>
void PatchTool<T1,T2>::PatchIntensityNormalizationUsingMeanAndVariance(Patch<T1> & inputPatch, Patch<T1> & refPatch, Patch<T2> & outputPatch)
{

  typename itkStatisticsT1ImageFilter::Pointer inputStat = itkStatisticsT1ImageFilter::New ();
  inputStat->SetInput(inputPatch.GetData());
  inputStat->Update();
  typename itkStatisticsT1ImageFilter::Pointer refStat = itkStatisticsT1ImageFilter::New ();
  refStat->SetInput(refPatch.GetData());
  refStat->Update();
  
  if( fabs(inputStat->GetMean() - inputStat->GetSigma()) > 0.0000001 )
  {
  
  	float a = (refStat->GetMean() - refStat->GetSigma()) / ( inputStat->GetMean() - inputStat->GetSigma() );
  	float b = refStat->GetMean() - a * inputStat->GetMean();
  
  	float minimum = refStat->GetMinimum();
  	float maximum = refStat->GetMaximum();
  
  	itkT1Iterator inputIt(inputPatch.GetData(), inputPatch.GetData()->GetLargestPossibleRegion());
  	itkT2Iterator outputIt(outputPatch.GetData(), outputPatch.GetData()->GetLargestPossibleRegion());

  	for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
  	{
  		float newValue = a*inputIt.Get()+b;
  		if(newValue < minimum)
  			newValue = minimum;
  		if(newValue > maximum)
  			newValue = maximum;
    	outputIt.Set( newValue );
  	}
	}
	else
	{
		itkT1Iterator inputIt(inputPatch.GetData(), inputPatch.GetData()->GetLargestPossibleRegion());
  	itkT2Iterator outputIt(outputPatch.GetData(), outputPatch.GetData()->GetLargestPossibleRegion());

  	for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
  	{
    	outputIt.Set( inputIt.Get() ); //what should we do? Impose new value or keep the old value ????
  	}
	}   
}

template<typename T1, typename T2>
void PatchTool<T1,T2>::CreatePatchFromAnotherPatch(Patch<T1> & inputPatch, Patch<T2> & outputPatch)
{
	outputPatch.SetImageSpacing(           inputPatch.GetImageSpacing()    );
	outputPatch.SetImageSize(              inputPatch.GetImageSize()       );
	outputPatch.SetImageRegion(            inputPatch.GetImageRegion()     );
	outputPatch.SetHalfPatchSize(          inputPatch.GetHalfPatchSize()   );
	outputPatch.SetFullPatchSize(          inputPatch.GetFullPatchSize()   );
	outputPatch.SetFullPatchRegion(        inputPatch.GetFullPatchRegion() );
	outputPatch.SetCentralPoint(           inputPatch.GetCentralPoint()    );
	outputPatch.SetCentralPointInImage(    inputPatch.GetCentralPointInImage()    );
	
	outputPatch.CreatePatch();
}

template<typename T1, typename T2>
void PatchTool<T1,T2>::ComputeSearchRegion(Patch<T1> & inputPatch, typename itkT1Image::RegionType & region)
{
  //create an appropiate search region around the current pixel p
  typename itkT1Image::RegionType::IndexType start;
  typename itkT1Image::RegionType::SizeType size;
  
  typename itkT1Image::RegionType::IndexType point = inputPatch.GetCentralPointInImage();
    
  for(unsigned int i=0; i!= point.GetIndexDimension(); i++){
    start[i] = point[i] - m_HalfSpatialBandwidth[i];
    size[i] = m_FullSpatialBandwidth[i];
        
    if(start[i] < 0){                  //if the starting index is outside the image (<0)
      size[i] += start[i];
      start[i] = 0;
    }
    else if(start[i] >= (int)m_ImageSize[i]){    //if the starting index is outside the image (>image size)
      start[i] = m_ImageSize[i]-1;
      size[i] = 0;
    }
    
    int d = (start[i] + size[i]) - m_ImageSize[i]; //if the region is not fully inside the image
    if(d>0){
//      size[i] = size[i] - d;
      if(static_cast< int >(size[i]) - d < 0)
      {
          size[i] = 0;
      }
      else
        size[i] = size[i] - d;

    }    
  }
  region.SetSize( size );
  region.SetIndex( start );

}		

template<typename T1, typename T2>
void PatchTool<T1,T2>::SetSelectionNeighbourPatchMethod(Patch<T1> & inputPatch)
{
  if(m_PatchSelectionMethod > 2) 
    m_PatchSelectionMethod = 0;
    
  switch(m_PatchSelectionMethod)
  {
    case 0: {
            //std::cout<<"patch selection method: chi 2"<<std::endl; 
            itk::Statistics::ChiSquareDistribution::Pointer chi2 = itk::Statistics::ChiSquareDistribution::New();
            int numberOfPointsToUseForPatchDistance = (2*inputPatch.GetHalfPatchSize()[0]+1) * (2*inputPatch.GetHalfPatchSize()[1]+1) * (2*inputPatch.GetHalfPatchSize()[2]+1);
            chi2->SetDegreesOfFreedom(numberOfPointsToUseForPatchDistance);
            m_ParamPatchSelectionMethod = exp(- chi2->EvaluateInverseCDF(m_ParamPatchSelectionMethod)/numberOfPointsToUseForPatchDistance);
            //std::cout<<m_ParamPatchSelectionMethod<<std::endl;
            break;
            }  
    case 1: //std::cout<<"patch selection method: percentage of nearest patches"<<std::endl; 
            break; 
    case 2: //std::cout<<"patch selection method: mean and variance"<<std::endl; 
            break; 
  }


}


template<typename T1, typename T2>
void PatchTool<T1,T2>::GetNeighbourPatches(Patch<T1> & inputPatch, std::vector< Patch<T2> > & neighbourPatches, itkT2ImagePointer & image)
{
  neighbourPatches.clear();
  
  //Set the search region for neighbours
  typename itkT2Image::RegionType searchRegion;
  this->ComputeSearchRegion(inputPatch, searchRegion);
  			
  //Go through the neighbourhood with a region iterator
  itkT2IteratorWithIndex it(image, searchRegion);
  typename itkT2IteratorWithIndex::IndexType neighbourPixelIndex;  
  
  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
                              
    neighbourPixelIndex = it.GetIndex();
              
    btk::Patch<T2> p;
    p.Initialize(image, inputPatch.GetHalfPatchSize() );
    p.ComputePatch(neighbourPixelIndex, image);
  
  	neighbourPatches.push_back(p);
  }
}

template<typename T1, typename T2>
void PatchTool<T1,T2>::ComputeNeighbourWeights(Patch<T1> & inputPatch, std::vector< Patch<T2> > & neighbourPatches, std::vector<float> & weights, float & smoothing)
{
  weights.resize( neighbourPatches.size() );
  switch(m_PatchSelectionMethod)
  {
    case 0: {//CHI 2 selection
            for(unsigned int i=0; i < neighbourPatches.size(); i++)
            {
                          
              float w = exp( - ComputeL2NormBetweenPatches(inputPatch, neighbourPatches[i]) / smoothing);
              //std::cout<<"w:"<<w<<", smoothing : "<<smoothing<<" "<<ComputeL2NormBetweenPatches(inputPatch, neighbourPatches[i])<<" "<<m_ParamPatchSelectionMethod<<std::endl;
                            
              if( w < m_ParamPatchSelectionMethod)
              	w = 0;
              weights[i] = w;
            }		
            break;
            }
    case 1: {//Nearest patches
            // TO BE DONE
            for(unsigned int i=0; i < neighbourPatches.size(); i++)
                weights[i] = exp( - ComputeL2NormBetweenPatches(inputPatch, neighbourPatches[i]) / smoothing);
            break;
            }
    case 2: {//Selection using mean and variance
    		for(unsigned int i=0; i < neighbourPatches.size(); i++)
    		{ 
    		  float w = 0;
    		  float meanRatio = 0;
    		  float stdDevRatio=0;
    		  
    		  if( neighbourPatches[i].GetMeanValue() == 0 )
    		  {
    		  	if(inputPatch.GetMeanValue() == 0)
    		  		meanRatio = 1;
    		  }
    		  else
    		  	meanRatio = inputPatch.GetMeanValue() / neighbourPatches[i].GetMeanValue();
    		  	
       		  if( neighbourPatches[i].GetStdDevValue() == 0 )
    		  {
    		  	if(inputPatch.GetStdDevValue() == 0)
    		  		stdDevRatio = 1;
    		  }
    		  else
    		  	stdDevRatio = inputPatch.GetStdDevValue() / neighbourPatches[i].GetStdDevValue(); 		  

              if( ((meanRatio > m_ParamPatchSelectionMethod) && (meanRatio < 1/m_ParamPatchSelectionMethod) && (stdDevRatio > 0.5)) || (stdDevRatio < 2) )
			    w = exp( - ComputeL2NormBetweenPatches(inputPatch, neighbourPatches[i]) / smoothing);

			  weights[i] = w;
    		}
            break;
            }
  }                  
}

template<typename T1, typename T2>
double PatchTool<T1,T2>::ComputeL2NormBetweenPatches(Patch<T1> & p, Patch<T2> & q)
{
  double diff=0;
  double dist = 0;
  itkConstT1Iterator itp( p.GetData(), p.GetData()->GetLargestPossibleRegion() );
  itkConstT2Iterator itq( q.GetData(), q.GetData()->GetLargestPossibleRegion() );

  for(itp.GoToBegin(), itq.GoToBegin(); !itp.IsAtEnd(); ++itp, ++itq){
    diff = itp.Get() - itq.Get();
    dist += diff*diff;	
  }
  return dist;   
}



} // namespace btk

