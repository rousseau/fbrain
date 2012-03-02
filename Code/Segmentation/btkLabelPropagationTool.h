/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

20 May 2011
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

#ifndef LabelFusionTool_H
#define LabelFusionTool_H


#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>

template <typename T>
class LabelFusionTool
{
 public:
  typedef typename itk::Image< T, 3> itkTImage;
  typedef typename itkTImage::Pointer itkTPointer;
  typedef typename itk::ImageDuplicator< itkTImage > itkTDuplicator;
  typedef typename itk::ImageRegionIterator< itkTImage > itkTIterator;
  typedef typename itk::ImageRegionConstIterator< itkTImage > itkTConstIterator;
  typedef typename itk::ImageRegionIteratorWithIndex< itkTImage > itkTIteratorWithIndex;
  typedef typename itk::ImageFileReader< itkTImage >  itkTReader;

  typedef typename itk::Image< float, 3> itkFloatImage;
  typedef typename itkFloatImage::Pointer itkFloatPointer;
  typedef typename itk::ImageRegionIterator< itkFloatImage > itkFloatIterator;
  typedef typename itk::ImageRegionConstIterator< itkFloatImage > itkFloatConstIterator;
  typedef typename itk::ImageRegionIteratorWithIndex< itkFloatImage > itkFloatIteratorWithIndex;

  typedef typename itk::Image< std::map<T, float>, 3> itkMapImage;
  typedef typename itkMapImage::Pointer itkMapPointer;
  typedef typename itk::ImageRegionIterator< itkMapImage > itkMapIterator;

  itkTPointer     m_inputImage;
  itkTPointer     m_outputImage;
  itkMapPointer   m_labelFusionImage;
  itkTPointer     m_maskImage;
  itkFloatPointer m_meanImage;
  itkFloatPointer m_varianceImage;  
  itkFloatPointer m_weightImage;

  std::vector<itkTPointer> m_anatomicalImages;
  std::vector<itkTPointer> m_labelImages;
  std::vector<itkFloatPointer> m_meanAnatomicalImages;
  std::vector<itkFloatPointer> m_varianceAnatomicalImages;

  float m_padding;
  float m_rangeBandwidth;
  int   m_blockwise;
  int   m_optimized;
  int   m_centralPointStrategy;
  float m_lowerMeanThreshold;
  float m_lowerVarianceThreshold;
  int   m_aggregation;
  int   m_normalization;
  int   m_minLabel;

  typename itkTImage::SizeType m_halfPatchSize;          //half of the patch size
  typename itkTImage::SizeType m_fullPatchSize;          //patch size  : 2 * halfPatchSize + 1
  typename itkTImage::SizeType m_halfSpatialBandwidth;   //equivalent to the half size of the volume search area in non-local means
  typename itkTImage::SizeType m_fullSpatialBandwidth;   //spatial bandwidth : 2 * halfSpatialBandwidth + 1

  itkTPointer m_hwnImage;
  itkTPointer m_hwvsImage;


  //Image information (size, spacing etc.)
  typename itkTImage::SpacingType m_spacing;
  typename itkTImage::SizeType    m_size;
  typename itkTImage::RegionType  m_region;

  void ReadInput(std::string input_file);
  void ReadAnatomicalImages(std::vector<std::string> & input_file);
  void ReadLabelImages(std::vector<std::string> & input_file);
  void SetPaddingValue(float padding);
  void SetMaskImage(itkTPointer maskImage);
  void SetPatchSize(int h);
  void SetSpatialBandwidth(int s);
  void SetSmoothing(float beta);
  void SetCentralPointStrategy(int s);
  void SetBlockwiseStrategy(int b);
  void SetAggregationStrategy(int a);
  void SetOptimizationStrategy(int o);
  void SetLowerThresholds(float m, float v);
  void SetNormalizationStrategy(int n);
  void SetMinLabel(int l);

  void CreatePatch(itkTPointer & patch);
  void CreatePatch(itkFloatPointer & patch);
  void CreatePatch(itkMapPointer & patch);
  void GetPatch(typename itkTImage::IndexType p, itkTPointer & patch, itkTPointer & image);
  void GetNormalizedPatch(typename itkTImage::IndexType p, itkFloatPointer & patch, itkTPointer & image, float mean, float stddev);
  void Get2Patches(typename itkTImage::IndexType p, itkTPointer & anatomicalPatch, itkTPointer & anatomicalImage, itkTPointer & labelPatch, itkTPointer & labelImage);
  void ComputePatchRegion(typename itkTImage::IndexType p, typename itkTImage::RegionType & imageRegion, typename itkTImage::RegionType & patchRegion);
  void ComputeSearchRegion(typename itkTImage::IndexType p, typename itkTImage::RegionType & region);
  void InitImage(itkFloatPointer & image);
  bool CheckSpeed(typename itkTImage::IndexType p, typename itkTImage::IndexType q, int i);

  void ComputeOutput();
  void ComputeOutput_groupwise();
  void ComputeDenoisedOutput();
  void ComputeHROutput();
  double GetDenoisedPatch(typename itkTImage::IndexType p, itkFloatPointer & patch);
  double GetHRPatch(typename itkTImage::IndexType p, itkFloatPointer & patch);
  void AddPatchToImage(typename itkTImage::IndexType p, itkFloatPointer & patch, itkFloatPointer & image, itkFloatPointer & weightImage, double weight);
  double GetLabelPatch(typename itkTImage::IndexType p, itkTPointer & patch);
  void GetFuzzyLabelPatch(typename itkTImage::IndexType p, itkTPointer & patch, std::vector< std::map<T, float> > & vectorMap);
  double GetLabelPatchUsingNormalization(typename itkTImage::IndexType p, itkTPointer & patch);
  void AddLabelPatchToLabelImage(typename itkTImage::IndexType p, itkTPointer & patch, double weight);
  void AddFuzzyLabelPatchToLabelImage(typename itkTImage::IndexType p, itkTPointer & patch, std::vector< std::map<T, float> > & vectorMap);
  double PatchDistance(itkTPointer & p,itkTPointer & q);
  double PatchDistance(itkFloatPointer & p,itkFloatPointer & q);
  void GetOutput(itkTPointer & outputImage);
  void GetWeightImage(itkTPointer & outputImage);
	void GetFuzzyWeightImage(itkFloatPointer & outputImage, int label);

  void ReadOneImage(std::string input_file, itkTPointer & image);

};


template <typename T>
void LabelFusionTool<T>::ReadInput(std::string input_file)
{
  typename itkTReader::Pointer reader = itkTReader::New();
  reader->SetFileName( input_file );
  reader->Update();
  m_inputImage = reader->GetOutput();

  //compute characteristics of the input image
  m_region  = m_inputImage->GetLargestPossibleRegion();
  m_size    = m_region.GetSize();
  m_spacing = m_inputImage->GetSpacing();
  
  //duplicate the input image into the output images to keep all header information
  typename itkTDuplicator::Pointer duplicator = itkTDuplicator::New();
  duplicator->SetInputImage( m_inputImage );
  duplicator->Update();
  m_outputImage = duplicator->GetOutput();
  m_outputImage->FillBuffer(0); 

  typename itk::CastImageFilter< itkTImage, itkFloatImage >::Pointer castFilter = itk::CastImageFilter< itkTImage, itkFloatImage >::New();
  castFilter->SetInput(m_outputImage);
  castFilter->Update();
  m_weightImage = castFilter->GetOutput();

  m_labelFusionImage = itkMapImage::New();
  m_labelFusionImage->SetRegions(m_inputImage->GetLargestPossibleRegion());
  m_labelFusionImage->SetSpacing( m_inputImage->GetSpacing() );
  m_labelFusionImage->SetOrigin( m_inputImage->GetOrigin() );
  m_labelFusionImage->SetDirection( m_inputImage->GetDirection() );
  m_labelFusionImage->Allocate();
}

template <typename T>
void LabelFusionTool<T>::ReadAnatomicalImages(std::vector<std::string> & input_file)
{
  m_anatomicalImages.resize(input_file.size());

  typename itkTImage::RegionType region;
  typename itkTImage::SizeType size;
  typename itkTImage::SpacingType spacing;

  for(unsigned int i=0;i<input_file.size();i++){
    //std::cout<<"Reading Anatomical Input Image : "<<input_file[i]<<"\n";

    //we have to declare the itk reader inside the loop (otherwise, if image sizes are different, it crashes)
    typename itkTReader::Pointer reader = itkTReader::New();
    reader->SetFileName( input_file[i]  );
    reader->Update();
    m_anatomicalImages[i] = reader->GetOutput();

    //Print image information (size and spacing)
    region = m_anatomicalImages[i]->GetLargestPossibleRegion();
    size = region.GetSize();
    //std::cout<<"Image size is : "<<size[0]<<" "<<size[1]<<" "<<size[2]<<"\n";
    spacing = m_anatomicalImages[i]->GetSpacing();
    //std::cout<<"Spacing of inputImage (in mm) : "<<spacing[0]<<", "<<spacing[1]<<", "<<spacing[2]<<"\n";
    //std::cout<<"------------------------------\n";
  }

  //Check if spacing and image size are not different
  for(unsigned int i=0;i<input_file.size();i++){
    region = m_anatomicalImages[i]->GetLargestPossibleRegion();
    size = region.GetSize();
    spacing = m_anatomicalImages[i]->GetSpacing();

    typename itkTImage::IndexType q;
    for(unsigned int j=0; j!= q.GetIndexDimension(); j++)
      if( (size[j] != m_size[j]) || (spacing[j] != m_spacing[j]) ){
        std::cout<<"*************************************************************************************\n";
        std::cout<<"WARNING : the size or the spacing of the anatomical image #"<<i+1<<" are incorrect wrt the input image\n";
        std::cout<<"*************************************************************************************\n";      
        break;
      }
  }

}

template <typename T>
void LabelFusionTool<T>::ReadLabelImages(std::vector<std::string> & input_file)
{
  m_labelImages.resize(input_file.size());

  typename itkTImage::RegionType region;
  typename itkTImage::SizeType size;
  typename itkTImage::SpacingType spacing;

  for(unsigned int i=0;i<input_file.size();i++){
    //std::cout<<"Reading Label Input Image : "<<input_file[i]<<"\n";

    //we have to declare the itk reader inside the loop (otherwise, if image sizes are different, it crashes)
    typename itkTReader::Pointer reader = itkTReader::New();
    reader->SetFileName( input_file[i]  );
    reader->Update();
    m_labelImages[i] = reader->GetOutput();

    //Print image information (size and spacing)
    region = m_labelImages[i]->GetLargestPossibleRegion();
    size = region.GetSize();
    //std::cout<<"Image size is : "<<size[0]<<" "<<size[1]<<" "<<size[2]<<"\n";
    spacing = m_labelImages[i]->GetSpacing();
    //std::cout<<"Spacing of inputImage (in mm) : "<<spacing[0]<<", "<<spacing[1]<<", "<<spacing[2]<<"\n";
    //std::cout<<"------------------------------\n";
  }

  //Check if spacing and image size are not different
  for(unsigned int i=0;i<input_file.size();i++){
    region = m_labelImages[i]->GetLargestPossibleRegion();
    size = region.GetSize();
    spacing = m_labelImages[i]->GetSpacing();

    typename itkTImage::IndexType q;
    for(unsigned int j=0; j!= q.GetIndexDimension(); j++)
      if( (size[j] != m_size[j]) || (spacing[j] != m_spacing[j]) ){
        std::cout<<"*************************************************************************************\n";
        std::cout<<"WARNING : the size or the spacing of the label image #"<<i+1<<" are incorrect wrt the input image\n";
        std::cout<<"*************************************************************************************\n";      
        break;
      }
  }

}

template <typename T>
void LabelFusionTool<T>::SetPaddingValue(float padding)
{
  m_padding = padding;
  
  std::cout<<"Creating the mask image using the padding value ("<<padding<<")\n";
  typename itkTDuplicator::Pointer duplicator = itkTDuplicator::New();
  duplicator->SetInputImage( m_inputImage );
  duplicator->Update();
  m_maskImage = duplicator->GetOutput();

  double count = 0;
  itkTIterator iterator( m_maskImage, m_maskImage->GetRequestedRegion() );
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator){
    if(iterator.Get() <= padding)
      iterator.Set(0);
    else{
      iterator.Set(1);
      count++;
    }
  }
  std::cout<<"Percentage of points to be processed : "<<(int)(count / m_maskImage->GetLargestPossibleRegion().GetNumberOfPixels() * 100.0) <<"\n";
  
}

template <typename T>
void LabelFusionTool<T>::SetMaskImage(itkTPointer maskImage)
{
  m_maskImage = maskImage; 
  
  //check the size and spacing of the mask wrt the input image
  typename itkTImage::SpacingType spacing;
  typename itkTImage::SizeType    size;
  typename itkTImage::RegionType  region;
  
  region  = m_maskImage->GetLargestPossibleRegion();
  size    = m_region.GetSize();
  spacing = m_maskImage->GetSpacing(); 
  
  typename itkTImage::IndexType q;
  for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
    if( (size[i] != m_size[i]) || (spacing[i] != m_spacing[i]) ){
      std::cout<<"*************************************************************************************\n";
      std::cout<<"WARNING : the size or the spacing of the mask image are incorrect wrt the input image\n";
      std::cout<<"*************************************************************************************\n";      
      break;
    }
  
  double count = 0;
  itkTIterator iterator( m_maskImage, m_maskImage->GetRequestedRegion() );
  for(iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator){
    if(iterator.Get() > 0)
      count++;
  }
  std::cout<<"Percentage of points to be processed : "<<(int)(count / m_maskImage->GetLargestPossibleRegion().GetNumberOfPixels() * 100.0) <<"\n";
  
  
}

template <typename T>
void LabelFusionTool<T>::SetPatchSize(int h)
{
  std::cout<<"Computing patch size (taking into account possible image anisotropy)\n";    
  float minVoxSz = m_spacing[0];
  if(m_spacing[1] < minVoxSz) minVoxSz = m_spacing[1];
  if(m_spacing[2] < minVoxSz) minVoxSz = m_spacing[2];

  m_halfPatchSize[0] = (int)(0.5 + h * minVoxSz / m_spacing[0]);
  m_halfPatchSize[1] = (int)(0.5 + h * minVoxSz / m_spacing[1]);
  m_halfPatchSize[2] = (int)(0.5 + h * minVoxSz / m_spacing[2]);
  std::cout<<"patchSize : "<<m_halfPatchSize[0]<<" "<<m_halfPatchSize[1]<<" "<<m_halfPatchSize[2]<<"\n";
  
  m_fullPatchSize[0] = 2 * m_halfPatchSize[0] + 1;
  m_fullPatchSize[1] = 2 * m_halfPatchSize[1] + 1;
  m_fullPatchSize[2] = 2 * m_halfPatchSize[2] + 1;
}

template <typename T>
void LabelFusionTool<T>::SetSpatialBandwidth(int s)
{
  std::cout<<"Computing spatial bandwidth (taking into account possible image anisotropy)\n";    
  float minVoxSz = m_spacing[0];
  if(m_spacing[1] < minVoxSz) minVoxSz = m_spacing[1];
  if(m_spacing[2] < minVoxSz) minVoxSz = m_spacing[2];

  m_halfSpatialBandwidth[0] = (int)(0.5 + s * minVoxSz / m_spacing[0]);
  m_halfSpatialBandwidth[1] = (int)(0.5 + s * minVoxSz / m_spacing[1]);
  m_halfSpatialBandwidth[2] = (int)(0.5 + s * minVoxSz / m_spacing[2]);
  std::cout<<"spatialBandwidth : "<<m_halfSpatialBandwidth[0]<<" "<<m_halfSpatialBandwidth[1]<<" "<<m_halfSpatialBandwidth[2]<<"\n";  
  
  m_fullSpatialBandwidth[0] = 2 * m_halfSpatialBandwidth[0] + 1;
  m_fullSpatialBandwidth[1] = 2 * m_halfSpatialBandwidth[1] + 1;
  m_fullSpatialBandwidth[2] = 2 * m_halfSpatialBandwidth[2] + 1;
}

template <typename T>
void LabelFusionTool<T>::SetSmoothing(float beta)
{
  //this function should be rewritten using a convolution-based approach
  std::cout<<"Computing the range bandwidth (corresponding to the smoothing parameter for the NLM algorithm).\n";
  std::cout<<"Implicit assumption : noise level is equal for all the images\n";
  unsigned int count = 0;
  double sigma2 = 0;
  int x,y,z;
  typename itkTImage::IndexType pixelIndex;	
  float ei = 0;
  std::vector<float> vecei;
  T value = 0;
  //since we have use to use a neighborhood around the current voxel, we neglect the border to avoid slow tests.
  for(z=1;z<(int)m_size[2]-1;z++)
    for(y=1;y<(int)m_size[1]-1;y++)
      for(x=1;x<(int)m_size[0]-1;x++){
	pixelIndex[0] = x;
	pixelIndex[1] = y;
	pixelIndex[2] = z;
	value = m_inputImage->GetPixel(pixelIndex);
        if( m_maskImage->GetPixel(pixelIndex) > 0){
	  //pourrait se faire avec une convolution !
	  pixelIndex[0] = x+1;	  pixelIndex[1] = y;	  pixelIndex[2] = z;
	  ei = m_inputImage->GetPixel(pixelIndex);
	  pixelIndex[0] = x-1;	  pixelIndex[1] = y;	  pixelIndex[2] = z;
	  ei += m_inputImage->GetPixel(pixelIndex);
	  pixelIndex[0] = x;	          pixelIndex[1] = y+1;	  pixelIndex[2] = z;
	  ei += m_inputImage->GetPixel(pixelIndex);
	  pixelIndex[0] = x;  	  pixelIndex[1] = y-1;	  pixelIndex[2] = z;
	  ei += m_inputImage->GetPixel(pixelIndex);
	  pixelIndex[0] = x;  	  pixelIndex[1] = y;	  pixelIndex[2] = z+1;
	  ei += m_inputImage->GetPixel(pixelIndex);
	  pixelIndex[0] = x;  	  pixelIndex[1] = y;	  pixelIndex[2] = z-1;
	  ei += m_inputImage->GetPixel(pixelIndex);
	  pixelIndex[0] = x;            pixelIndex[1] = y;	  pixelIndex[2] = z;
	  ei = sqrt(6.0/7.0)*(m_inputImage->GetPixel(pixelIndex) -ei/6.0);
	    
	  sigma2 += ei*ei;
	  count ++;
          if(fabs(ei>0))
            vecei.push_back(fabs(ei));
        }
      }
  //Estimation of sigma with MAD
  std::sort(vecei.begin(), vecei.end());
  float med = vecei[(int)(vecei.size()/2)];
  for(unsigned int i=0; i<vecei.size(); i++)
    vecei[i] = fabs(vecei[i] - med);
  std::sort(vecei.begin(), vecei.end());
    
  sigma2 = 1.4826 * vecei[(int)(vecei.size()/2)];
  sigma2 = sigma2 * sigma2;
  std::cout<<"Number of points for estimation : "<<count<<"\n";
  std::cout<<"Estimated Global Sigma : "<<sqrt(sigma2)<<"\n";
  float NLMsmooth = 2 * beta * sigma2 * (2*m_halfPatchSize[0]+1) * (2*m_halfPatchSize[1]+1) * (2*m_halfPatchSize[2]+1);
  std::cout<<"Global smoothing parameter h : "<<sqrt(NLMsmooth)<<"\n";
  m_rangeBandwidth = NLMsmooth;   
}

template <typename T>
void LabelFusionTool<T>::SetCentralPointStrategy(int s)
{
  m_centralPointStrategy = s;
}

template <typename T>
void LabelFusionTool<T>::SetBlockwiseStrategy(int b)
{
  m_blockwise = b;
}

template <typename T>
void LabelFusionTool<T>::SetAggregationStrategy(int a)
{
  m_aggregation = a;
}

template <typename T>
void LabelFusionTool<T>::SetOptimizationStrategy(int o)
{
  m_optimized = o;

  //if(o==1){
    std::cout<<"Optimized mode requirements: Computing Mean and Variance images\n";
    InitImage(m_meanImage);
    InitImage(m_varianceImage);

    int x,y,z;
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)m_size[2]; z++)
      for(y=0; y < (int)m_size[1]; y++)
        for(x=0; x < (int)m_size[0]; x++){
	  typename itkTImage::IndexType p;
	  p[0] = x;
	  p[1] = y;
	  p[2] = z;

  	  if( m_maskImage->GetPixel(p) > 0 ){	  
	    itkTPointer patch = itkTImage::New();
            CreatePatch(patch);
            GetPatch(p,patch,m_inputImage);
            itkTIterator itp( patch, patch->GetRequestedRegion() );
            double m = 0;
            double m2= 0;
            for(itp.GoToBegin(); !itp.IsAtEnd(); ++itp){
              m += itp.Get();
              m2+= (itp.Get() * itp.Get());	
            }
            int n = patch->GetLargestPossibleRegion().GetNumberOfPixels();
            float mean = m / n;
            float variance = (m2 / n) - (mean * mean) ;
            
            m_meanImage->SetPixel( p, mean );	  	  
            m_varianceImage->SetPixel( p, variance );	  	  
	  }
      }

  m_meanAnatomicalImages.resize(m_anatomicalImages.size());
  m_varianceAnatomicalImages.resize(m_anatomicalImages.size());

  for(unsigned int i=0; i < m_anatomicalImages.size(); i++){
    InitImage(m_meanAnatomicalImages[i]);
    InitImage(m_varianceAnatomicalImages[i]);

    int x,y,z;
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)m_size[2]; z++)
      for(y=0; y < (int)m_size[1]; y++)
        for(x=0; x < (int)m_size[0]; x++){
	  typename itkTImage::IndexType p;
	  p[0] = x;
	  p[1] = y;
	  p[2] = z;

  	  if( m_maskImage->GetPixel(p) > 0 ){	  
	    itkTPointer patch = itkTImage::New();
            CreatePatch(patch);
            GetPatch(p,patch,m_anatomicalImages[i]);
            itkTIterator itp( patch, patch->GetRequestedRegion() );
            double m = 0;
            double m2= 0;
            for(itp.GoToBegin(); !itp.IsAtEnd(); ++itp){
              m += itp.Get();
              m2+= (itp.Get() * itp.Get());	
            }
            int n = patch->GetLargestPossibleRegion().GetNumberOfPixels();
            float mean = m / n;
            float variance = (m2 / n) - (mean * mean) ;
            
            m_meanAnatomicalImages[i]->SetPixel( p, mean );	  	  
            m_varianceAnatomicalImages[i]->SetPixel( p, variance );	  	  
	  }
      }
  }

  //}
}

template <typename T>
void LabelFusionTool<T>::SetLowerThresholds(float m, float v)
{
  m_lowerMeanThreshold = m;
  m_lowerVarianceThreshold = v;
}

template <typename T>
void LabelFusionTool<T>::SetNormalizationStrategy(int n)
{
  m_normalization = n;
}

template <typename T>
void LabelFusionTool<T>::SetMinLabel(int l)
{
  m_minLabel = l;
}

template <typename T>
void LabelFusionTool<T>::CreatePatch(itkTPointer & patch)
{
  //resize and allocate the patch and set the estimate to 0
  typename itkTImage::SizeType size = m_fullPatchSize;
  typename itkTImage::RegionType region;
  region.SetSize(size);
  patch->SetRegions(region);
  patch->Allocate();
  patch->FillBuffer(0);  
}


template <typename T>
void LabelFusionTool<T>::CreatePatch(itkFloatPointer & patch)
{
  //resize and allocate the patch and set the estimate to 0
  typename itkFloatImage::SizeType size = m_fullPatchSize;
  typename itkFloatImage::RegionType region;
  region.SetSize(size);
  patch->SetRegions(region);
  patch->Allocate();
  patch->FillBuffer(0.0);  
}

template <typename T>
void LabelFusionTool<T>::CreatePatch(itkMapPointer & patch)
{
  //resize and allocate the patch and set the estimate to 0
  typename itkFloatImage::SizeType size = m_fullPatchSize;
  typename itkFloatImage::RegionType region;
  region.SetSize(size);
  patch->SetRegions(region);
  patch->Allocate();
}

template <typename T>
void LabelFusionTool<T>::GetPatch(typename itkTImage::IndexType p, itkTPointer & patch, itkTPointer & image)
{
  //this function is equivalent to a cropping function
  //WARNING : no check about the size of the input patch !
  patch->FillBuffer(0);
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkTConstIterator inputIt( image, imageRegion);
  itkTIterator outputIt( patch, patchRegion);

  for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
    outputIt.Set( inputIt.Get() );

}

template <typename T>
void LabelFusionTool<T>::GetNormalizedPatch(typename itkTImage::IndexType p, itkFloatPointer & patch, itkTPointer & image, float mean, float stddev)
{
  //this function is equivalent to a cropping function
  //WARNING : no check about the size of the input patch !
  patch->FillBuffer(0);

  if(stddev > 0){
    typename itkTImage::RegionType imageRegion;
    typename itkTImage::RegionType patchRegion;
    ComputePatchRegion(p,imageRegion,patchRegion);
  
    itkTConstIterator inputIt( image, imageRegion);
    itkFloatIterator outputIt( patch, patchRegion);

    for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
      outputIt.Set( (inputIt.Get() - mean) / stddev );
  }
}


template <typename T>
void LabelFusionTool<T>::Get2Patches(typename itkTImage::IndexType p, itkTPointer & anatomicalPatch, itkTPointer & anatomicalImage, itkTPointer & labelPatch, itkTPointer & labelImage)
{
  //this function is equivalent to a cropping function
  //WARNING : no check about the size of the input patch !
  anatomicalPatch->FillBuffer(0);
  labelPatch->FillBuffer(0);
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkTConstIterator inputAnatomicalIt( anatomicalImage, imageRegion);
  itkTConstIterator inputLabelIt( labelImage, imageRegion);
  itkTIterator anatomicalIt( anatomicalPatch, patchRegion);
  itkTIterator labelIt( labelPatch, patchRegion);

  for ( inputAnatomicalIt.GoToBegin(), inputLabelIt.GoToBegin(), anatomicalIt.GoToBegin(), labelIt.GoToBegin(); !inputAnatomicalIt.IsAtEnd(); ++inputAnatomicalIt, ++inputLabelIt, ++anatomicalIt, ++labelIt){
    anatomicalIt.Set( inputAnatomicalIt.Get() );
    labelIt.Set( inputLabelIt.Get() );
  }

}

template <typename T>
void LabelFusionTool<T>::ComputePatchRegion(typename itkTImage::IndexType p, typename itkTImage::RegionType & imageRegion, typename itkTImage::RegionType & patchRegion)
{
  //create an appropriate patch region around the current pixel p, and consider an offset for pixels close to boundaries.
  //the patch region is defined in the image coordinate system.
  //the offset allows to do the link between the patch coordinate system and the one of the image.
  typename itkTImage::RegionType::IndexType start;
  typename itkTImage::RegionType::IndexType offset;
  typename itkTImage::RegionType::SizeType size;

  for(unsigned int i=0; i!= p.GetIndexDimension(); i++){
    start[i] = p[i] - m_halfPatchSize[i];
    size[i] = m_fullPatchSize[i];
    offset[i] = 0;
    
    if(start[i] < 0){                  //if the starting index is outside the image (<0)
      size[i] += start[i];
      offset[i] = -start[i];
      start[i] = 0;
    }
    else if(start[i] >= (int)m_size[i]){    //if the starting index is outside the image (>image size)
      start[i] = m_size[i]-1;
      size[i] = 0;
    }
    
    int d = (start[i] + size[i]) - m_size[i]; //if the region is not fully inside the image
    if(d>0){
      size[i] = size[i] - d;
      if(size[i] < 0) size[i] = 0;
    }    
  }  
  imageRegion.SetSize( size );
  imageRegion.SetIndex( start );    
  patchRegion.SetSize( size );
  patchRegion.SetIndex( offset );    
}

template <typename T>
void LabelFusionTool<T>::ComputeSearchRegion(typename itkTImage::IndexType p, typename itkTImage::RegionType & region)
{
  //create an appropiate search region around the current pixel p
  typename itkTImage::RegionType::IndexType start;
  typename itkTImage::RegionType::SizeType size;
  
  for(unsigned int i=0; i!= p.GetIndexDimension(); i++){
    start[i] = p[i] - m_halfSpatialBandwidth[i];
    size[i] = m_fullSpatialBandwidth[i];
        
    if(start[i] < 0){                  //if the starting index is outside the image (<0)
      size[i] += start[i];
      start[i] = 0;
    }
    else if(start[i] >= (int)m_size[i]){    //if the starting index is outside the image (>image size)
      start[i] = m_size[i]-1;
      size[i] = 0;
    }
    
    int d = (start[i] + size[i]) - m_size[i]; //if the region is not fully inside the image
    if(d>0){
      size[i] = size[i] - d;
      if(size[i] < 0) size[i] = 0;
    }    
  }
  region.SetSize( size );
  region.SetIndex( start );
}

template <typename T>
void LabelFusionTool<T>::InitImage(itkFloatPointer & image)
{
  image = itkFloatImage::New();
  image->SetRegions(m_inputImage->GetLargestPossibleRegion());
  image->SetSpacing( m_inputImage->GetSpacing() );
  image->SetOrigin( m_inputImage->GetOrigin() );
  image->SetDirection( m_inputImage->GetDirection() );
  image->Allocate();
  image->FillBuffer(0); 
}

template <typename T>
bool LabelFusionTool<T>::CheckSpeed(typename itkTImage::IndexType p, typename itkTImage::IndexType q, int i)
{
  bool goForIt = true;

  float mSpeed = 0;
  if(m_meanAnatomicalImages[i]->GetPixel(q) ==0){
    if(m_meanImage->GetPixel(p) == 0)
      mSpeed = 1;
    else
      mSpeed = 0;
  }
  else
    mSpeed = m_meanImage->GetPixel(p) / m_meanAnatomicalImages[i]->GetPixel(q);

  if( (mSpeed < m_lowerMeanThreshold) || (mSpeed > 1/m_lowerMeanThreshold) )
    goForIt = false;

  float vSpeed = 0;
  if(m_varianceAnatomicalImages[i]->GetPixel(q) ==0){
    if(m_varianceImage->GetPixel(p) == 0)
      vSpeed = 1;
    else
      vSpeed = 0;
  }
  else
    vSpeed = m_varianceImage->GetPixel(p) / m_varianceAnatomicalImages[i]->GetPixel(q);

  if( (vSpeed < m_lowerVarianceThreshold) || (vSpeed > 1/m_lowerVarianceThreshold) )
    goForIt = false;

  return goForIt;
}

template <typename T>
void LabelFusionTool<T>::ComputeOutput()
{
  int x,y,z;
  itkTIterator maskImageIt( m_maskImage, m_maskImage->GetLargestPossibleRegion());
  itkTIterator outputImageIt( m_outputImage, m_outputImage->GetLargestPossibleRegion());
  itkFloatIterator weightImageIt( m_weightImage, m_weightImage->GetLargestPossibleRegion());

	  
  if(m_blockwise == 0){
    std::cout<<"pointwise approach\n";
    typename itkTImage::IndexType q;
    for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
      q[i] = m_halfPatchSize[i];
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)m_size[2]; z++)
      for(y=0; y < (int)m_size[1]; y++)
        for(x=0; x < (int)m_size[0]; x++){
          typename itkTImage::IndexType p;
          p[0] = x;
          p[1] = y;
          p[2] = z;

          if( m_maskImage->GetPixel(p) > 0 ){	  
            itkTPointer patch = itkTImage::New();
            double wmax;
            if(m_normalization == 0) wmax = GetLabelPatch(p, patch);
            else wmax = GetLabelPatchUsingNormalization(p, patch);
            T label = patch->GetPixel(q); 
            if(label < m_minLabel){ label = m_minLabel; wmax = 0;}  //we threshold possible values. In this case, the associated weight is 0.

            double weight = 1.0; //images of the training set have the same weight (basic option. the other option could be to set weight equal to sum)
            if(m_aggregation == 0) weight = wmax;

            m_outputImage->SetPixel( p, label );	  	  
            m_weightImage->SetPixel( p, weight );	  	  
          }
      }
  }
  if(m_blockwise >= 1){
    typename itkTImage::IndexType q;
    for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
      q[i] = m_halfPatchSize[i];

    if(m_blockwise == 1){
      std::cout<<"blockwise approach\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        for(y=0; y < (int)m_size[1]; y++)
          for(x=0; x < (int)m_size[0]; x++){
	    typename itkTImage::IndexType p;
	    p[0] = x;
	    p[1] = y;
	    p[2] = z;

	    if( m_maskImage->GetPixel(p) > 0 ){	  
	      itkTPointer patch = itkTImage::New();
              double wmax;
              if(m_normalization == 0) wmax = GetLabelPatch(p, patch);
              else wmax = GetLabelPatchUsingNormalization(p, patch);
              double weight = 1.0; //images of the training set have the same weight (basic option. the other option could be to set weight equal to sum)
              if(m_aggregation == 0) weight = wmax;
              if(patch->GetPixel(q) != -1){ // label = -1 means no relevant example has been found
                #pragma omp critical
                AddLabelPatchToLabelImage(p, patch, weight);
              }
	    }
          }
    }
    if(m_blockwise == 2){
      std::cout<<"fast blockwise approach\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        if( z%(m_halfPatchSize[2]+1) == 0){
        for(y=0; y < (int)m_size[1]; y++)
          if( y%(m_halfPatchSize[1]+1) == 0){
          for(x=0; x < (int)m_size[0]; x++)
            if( x%(m_halfPatchSize[0]+1) == 0 ){
	      typename itkTImage::IndexType p;
	      p[0] = x;
	      p[1] = y;
	      p[2] = z;

	      if( m_maskImage->GetPixel(p) > 0 ){	  
  	        itkTPointer patch = itkTImage::New();
                double wmax;
                if(m_normalization == 0) wmax = GetLabelPatch(p, patch);
                else wmax = GetLabelPatchUsingNormalization(p, patch);
                double weight = 1.0; //images of the training set have the same weight (basic option. the other option could be to set weight equal to sum)
                if(m_aggregation == 0) weight = wmax;
                if(patch->GetPixel(q) != -1){ // label = -1 means no relevant example has been found
                  #pragma omp critical
                  AddLabelPatchToLabelImage(p, patch, weight);
                }
	      }
            }
          }
        }
    }

    std::cout<<"Find the highest weight and assign the final labels ...\n";
    itkMapIterator labelImageIt( m_labelFusionImage, m_labelFusionImage->GetLargestPossibleRegion());

    for(labelImageIt.GoToBegin(), maskImageIt.GoToBegin(), outputImageIt.GoToBegin(), weightImageIt.GoToBegin(); !labelImageIt.IsAtEnd(); ++labelImageIt, ++maskImageIt, ++outputImageIt, ++weightImageIt){
      if(maskImageIt.Get() > 0){
        std::map<T, float> map;
        typename std::map<T, float>::iterator mapIt;

        map = labelImageIt.Get();
        short label = -1;
        float wmax = 0;
        for(mapIt = map.begin (); mapIt != map.end (); ++mapIt){
          if( (*mapIt).second > wmax ){
            wmax = (*mapIt).second;
            label = (*mapIt).first;
          }
        }

        if(label < m_minLabel){ label = m_minLabel; wmax = 0;}  //we threshold possible values. In this case, the associated weight is 0.

        outputImageIt.Set(label);
        weightImageIt.Set(wmax);   //warning : the weights are un-normalized (i.e. divided by the number of estimates (which is not constant for the fastblock version) ).
      }
    }

  }

  //check whether some points are unlabeled (label==-1)
  float count = 0;
  for(maskImageIt.GoToBegin(), outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++maskImageIt, ++outputImageIt)
      if( (maskImageIt.Get() > 0) && (outputImageIt.Get() == -1) )
        count++;
  std::cout<<"Number of unlabeled points : "<<count<<"\n";

}

template <typename T>
void LabelFusionTool<T>::ComputeOutput_groupwise()
{
  int x,y,z;
  itkTIterator maskImageIt( m_maskImage, m_maskImage->GetLargestPossibleRegion());
  itkTIterator outputImageIt( m_outputImage, m_outputImage->GetLargestPossibleRegion());
  itkFloatIterator weightImageIt( m_weightImage, m_weightImage->GetLargestPossibleRegion());

	  
  if(m_blockwise == 0){
    std::cout<<"pointwise approach\n";
    typename itkTImage::IndexType q;
    for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
      q[i] = m_halfPatchSize[i];
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)m_size[2]; z++)
      for(y=0; y < (int)m_size[1]; y++)
        for(x=0; x < (int)m_size[0]; x++){
          typename itkTImage::IndexType p;
          p[0] = x;
          p[1] = y;
          p[2] = z;

          if( m_maskImage->GetPixel(p) > 0 ){	  
            itkTPointer patch = itkTImage::New();
            double wmax;
            if(m_normalization == 0) wmax = GetLabelPatch(p, patch);
            else wmax = GetLabelPatchUsingNormalization(p, patch);
            T label = patch->GetPixel(q); 
            if(label < m_minLabel){ label = m_minLabel; wmax = 0;}  //we threshold possible values. In this case, the associated weight is 0.
            m_outputImage->SetPixel( p, label );	  	  
            m_weightImage->SetPixel( p, wmax );	  	  
	  }
      }
  }
  if(m_blockwise >= 1){
    typename itkTImage::IndexType q;
    for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
      q[i] = m_halfPatchSize[i];

    if(m_blockwise == 1){
      std::cout<<"blockwise approach\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        for(y=0; y < (int)m_size[1]; y++)
          for(x=0; x < (int)m_size[0]; x++){
            typename itkTImage::IndexType p;
            p[0] = x;
            p[1] = y;
            p[2] = z;

            if( m_maskImage->GetPixel(p) > 0 ){	  
              itkTPointer patch = itkTImage::New();

              //vector (1D patch) to store the cumulative weight for each label
              std::vector< std::map<T, float> > vectorMap;
              GetFuzzyLabelPatch(p, patch, vectorMap);

              #pragma omp critical
              AddFuzzyLabelPatchToLabelImage(p, patch, vectorMap);

            }
          }
    }
    if(m_blockwise == 2){
      std::cout<<"fast blockwise approach\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        if( z%(m_halfPatchSize[2]+1) == 0){
        for(y=0; y < (int)m_size[1]; y++)
          if( y%(m_halfPatchSize[1]+1) == 0){
          for(x=0; x < (int)m_size[0]; x++)
            if( x%(m_halfPatchSize[0]+1) == 0 ){
              typename itkTImage::IndexType p;
              p[0] = x;
              p[1] = y;
              p[2] = z;

              if( m_maskImage->GetPixel(p) > 0 ){	  
                itkTPointer patch = itkTImage::New();

                //vector (1D patch) to store the cumulative weight for each label
                std::vector< std::map<T, float> > vectorMap;
                GetFuzzyLabelPatch(p, patch, vectorMap);

                #pragma omp critical
                AddFuzzyLabelPatchToLabelImage(p, patch, vectorMap);

              }
            }
          }
        }
    }

    std::cout<<"Find the highest weight and assign the final labels ...\n";
    itkMapIterator labelImageIt( m_labelFusionImage, m_labelFusionImage->GetLargestPossibleRegion());

    for(labelImageIt.GoToBegin(), maskImageIt.GoToBegin(), outputImageIt.GoToBegin(), weightImageIt.GoToBegin(); !labelImageIt.IsAtEnd(); ++labelImageIt, ++maskImageIt, ++outputImageIt, ++weightImageIt){
      if(maskImageIt.Get() > 0){
        std::map<T, float> map;
        typename std::map<T, float>::iterator mapIt;

        map = labelImageIt.Get();
        short label = -1;
        float wmax = 0;
        for(mapIt = map.begin (); mapIt != map.end (); ++mapIt){
          if( (*mapIt).second > wmax ){
            wmax = (*mapIt).second;
            label = (*mapIt).first;
          }
        }

        if(label < m_minLabel){ label = m_minLabel; wmax = 0;}  //we threshold possible values. In this case, the associated weight is 0.

        outputImageIt.Set(label);
        weightImageIt.Set(wmax);   //warning : the weights are un-normalized (i.e. divided by the number of estimates (which is not constant for the fastblock version) ).
      }
    }

  }

  //check whether some points are unlabeled (label==-1)
  float count = 0;
  for(maskImageIt.GoToBegin(), outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++maskImageIt, ++outputImageIt)
      if( (maskImageIt.Get() > 0) && (outputImageIt.Get() == -1) )
        count++;
  std::cout<<"Number of unlabeled points : "<<count<<"\n";

}


template <typename T>
void LabelFusionTool<T>::ComputeDenoisedOutput()
{
  std::cout<<"Compute the denoised image using NLM algorithm\n";
  
  itkFloatPointer denoisedImage = itkFloatImage::New();
  denoisedImage->SetRegions(m_inputImage->GetLargestPossibleRegion());
  denoisedImage->SetSpacing( m_inputImage->GetSpacing() );
  denoisedImage->SetOrigin( m_inputImage->GetOrigin() );
  denoisedImage->SetDirection( m_inputImage->GetDirection() );
  denoisedImage->Allocate();
  denoisedImage->FillBuffer(0); 
  
  itkFloatIterator denoisedIt( denoisedImage, denoisedImage->GetLargestPossibleRegion() );
  itkTIterator outputIt( m_outputImage, m_outputImage->GetLargestPossibleRegion());

  int x,y,z;
	  
  if(m_blockwise == 0){
    std::cout<<"pointwise denoising\n";
    typename itkTImage::IndexType q;
    for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
      q[i] = m_halfPatchSize[i];
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)m_size[2]; z++)
      for(y=0; y < (int)m_size[1]; y++)
        for(x=0; x < (int)m_size[0]; x++){
          typename itkTImage::IndexType p;
          p[0] = x;
          p[1] = y;
          p[2] = z;

          if( m_maskImage->GetPixel(p) > 0 ){	  
            itkFloatPointer patch = itkFloatImage::New();
            double sum = GetDenoisedPatch(p, patch);
            denoisedImage->SetPixel( p, patch->GetPixel(q) );
            m_weightImage->SetPixel( p, sum );	  	  	  	  
          }
      }
  }
  if(m_blockwise >= 1){

    if(m_blockwise == 1){
      std::cout<<"blockwise denoising\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        for(y=0; y < (int)m_size[1]; y++)
          for(x=0; x < (int)m_size[0]; x++){
            typename itkTImage::IndexType p;
            p[0] = x;
            p[1] = y;
            p[2] = z;

            if( m_maskImage->GetPixel(p) > 0 ){	  
              itkFloatPointer patch = itkFloatImage::New();
              double sum = GetDenoisedPatch(p, patch);             
              double weight = 1.0;
              if(m_aggregation == 0) weight = sum;
              #pragma omp critical
              AddPatchToImage(p, patch, denoisedImage, m_weightImage, weight);
            }
          }
    }
    if(m_blockwise == 2){
      std::cout<<"fast blockwise denoising\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        if( z%(m_halfPatchSize[2]+1) == 0){
        for(y=0; y < (int)m_size[1]; y++)
          if( y%(m_halfPatchSize[1]+1) == 0){
          for(x=0; x < (int)m_size[0]; x++)
            if( x%(m_halfPatchSize[0]+1) == 0 ){
              typename itkTImage::IndexType p;
              p[0] = x;
              p[1] = y;
              p[2] = z;

              if( m_maskImage->GetPixel(p) > 0 ){	  
                itkFloatPointer patch = itkFloatImage::New();
                double sum = GetDenoisedPatch(p, patch);
                double weight = 1.0;
                if(m_aggregation == 0) weight = sum;
                #pragma omp critical
                AddPatchToImage(p, patch, denoisedImage, m_weightImage, weight);
              }
            }
          }
        }
    }

    itkTIterator maskImageIt( m_maskImage, m_maskImage->GetLargestPossibleRegion() );
    itkFloatIterator weightIt( m_weightImage, m_weightImage->GetLargestPossibleRegion() );
    itkTIterator inputIt( m_inputImage, m_inputImage->GetLargestPossibleRegion() );
    //weight normalization
    for ( denoisedIt.GoToBegin(), weightIt.GoToBegin(), inputIt.GoToBegin(); !denoisedIt.IsAtEnd(); ++denoisedIt, ++weightIt, ++inputIt)
      if( weightIt.Get() > 0 )
        denoisedIt.Set( denoisedIt.Get() / weightIt.Get() );
      else
        denoisedIt.Set( inputIt.Get() );

    for(maskImageIt.GoToBegin(), denoisedIt.GoToBegin(); !denoisedIt.IsAtEnd(); ++maskImageIt, ++denoisedIt)
      if( (maskImageIt.Get() == 0) )
        denoisedIt.Set( 0 );
  
  }

  //Convert float data from denoisedImage to T data for m_outputImage
  for ( denoisedIt.GoToBegin(), outputIt.GoToBegin(); !denoisedIt.IsAtEnd(); ++denoisedIt, ++outputIt)
    outputIt.Set( (T)(denoisedIt.Get()) );

}

template <typename T>
void LabelFusionTool<T>::ComputeHROutput()
{
  std::cout<<"Compute an HR image using label propagation approach\n";
  
  itkFloatPointer HRImage = itkFloatImage::New();
  HRImage->SetRegions(m_inputImage->GetLargestPossibleRegion());
  HRImage->SetSpacing( m_inputImage->GetSpacing() );
  HRImage->SetOrigin( m_inputImage->GetOrigin() );
  HRImage->SetDirection( m_inputImage->GetDirection() );
  HRImage->Allocate();
  HRImage->FillBuffer(0); 
  
  itkFloatIterator HRIt( HRImage, HRImage->GetLargestPossibleRegion() );
  itkTIterator outputIt( m_outputImage, m_outputImage->GetLargestPossibleRegion());
  itkTIterator inputIt( m_inputImage, m_outputImage->GetLargestPossibleRegion());
  
  int x,y,z;
  
  if(m_blockwise == 0){
    std::cout<<"pointwise HR estimation\n";
    typename itkTImage::IndexType q;
    for(unsigned int i=0; i!= q.GetIndexDimension(); i++)
      q[i] = m_halfPatchSize[i];
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)m_size[2]; z++)
      for(y=0; y < (int)m_size[1]; y++)
        for(x=0; x < (int)m_size[0]; x++){
          typename itkTImage::IndexType p;
          p[0] = x;
          p[1] = y;
          p[2] = z;
          
          if( m_maskImage->GetPixel(p) > 0 ){	  
            itkFloatPointer patch = itkFloatImage::New();
            double sum = GetHRPatch(p, patch);
            HRImage->SetPixel( p, patch->GetPixel(q) );
            m_weightImage->SetPixel( p, sum );	  	  
            
          }
        }
  }
  if(m_blockwise >= 1){
    
    if(m_blockwise == 1){
      std::cout<<"blockwise HR estimation\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        for(y=0; y < (int)m_size[1]; y++)
          for(x=0; x < (int)m_size[0]; x++){
            typename itkTImage::IndexType p;
            p[0] = x;
            p[1] = y;
            p[2] = z;
            
            if( m_maskImage->GetPixel(p) > 0 ){	  
              itkFloatPointer patch = itkFloatImage::New();
              double sum = GetHRPatch(p, patch);             
              double weight = 1.0;
              if(m_aggregation == 0) weight = sum;
              #pragma omp critical
              AddPatchToImage(p, patch, HRImage, m_weightImage, weight);
            }
          }
    }
    if(m_blockwise == 2){
      std::cout<<"fast blockwise HR estimation\n";
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)m_size[2]; z++)
        if( z%(m_halfPatchSize[2]+1) == 0){
          for(y=0; y < (int)m_size[1]; y++)
            if( y%(m_halfPatchSize[1]+1) == 0){
              for(x=0; x < (int)m_size[0]; x++)
                if( x%(m_halfPatchSize[0]+1) == 0 ){
                  typename itkTImage::IndexType p;
                  p[0] = x;
                  p[1] = y;
                  p[2] = z;
                  
                  if( m_maskImage->GetPixel(p) > 0 ){	  
                    itkFloatPointer patch = itkFloatImage::New();
                    double sum = GetHRPatch(p, patch);
                    double weight = 1.0;
                    if(m_aggregation == 0) weight = sum;
                    #pragma omp critical
                    AddPatchToImage(p, patch, HRImage, m_weightImage, weight);
                  }
                }
            }
        }
    }
    
    itkTIterator maskImageIt( m_maskImage, m_maskImage->GetLargestPossibleRegion() );
    itkFloatIterator weightIt( m_weightImage, m_weightImage->GetLargestPossibleRegion() );
    //weight normalization    
    for ( HRIt.GoToBegin(), weightIt.GoToBegin(), inputIt.GoToBegin(); !HRIt.IsAtEnd(); ++HRIt, ++weightIt, ++inputIt)
      if( weightIt.Get() > 0 )
        HRIt.Set( HRIt.Get() / weightIt.Get() );
      else{
        switch(m_centralPointStrategy){
          case 0:
            HRIt.Set( 0 ); //It can be a specific value to explicitly say that no example has been found 
            break;
          case 1: 
            HRIt.Set( inputIt.Get() ); //Copy the low resolution value into the HR image
            break;
          default: 
            HRIt.Set( inputIt.Get() ); //Copy the low resolution value into the HR image
            break;
        }
      }
        
    for(maskImageIt.GoToBegin(), HRIt.GoToBegin(); !HRIt.IsAtEnd(); ++maskImageIt, ++HRIt)
      if( (maskImageIt.Get() == 0) )
        HRIt.Set( 0 );
    
  }
  
  //Convert float data from denoisedImage to T data for m_outputImage
  for ( HRIt.GoToBegin(), outputIt.GoToBegin(); !HRIt.IsAtEnd(); ++HRIt, ++outputIt)
    outputIt.Set( (T)(HRIt.Get()) );
  
}

template <typename T>
double LabelFusionTool<T>::GetDenoisedPatch(typename itkTImage::IndexType p, itkFloatPointer & patch)
{
  double wmax = 0; //maximum weight of patches
  double sum  = 0; //sum of weights (used for normalization purpose)
  
  //create the patch and set the estimate to 0
  CreatePatch(patch);
  itkFloatIterator patchIt(patch, patch->GetRequestedRegion());

  //get the value of the patch around the current pixel
  itkTPointer centralPatch = itkTImage::New();
  CreatePatch(centralPatch);
  GetPatch(p, centralPatch, m_inputImage);
  itkTIterator centralPatchIt(centralPatch, centralPatch->GetRequestedRegion());

  //itkFloatPointer normalizedCentralPatch = itkFloatImage::New();
  //CreatePatch(normalizedCentralPatch);
  //float mean = m_meanImage->GetPixel(p);
  //float stddev = sqrt( m_varianceImage->GetPixel(p) );
  //GetNormalizedPatch(p, normalizedCentralPatch, m_inputImage, mean, stddev);
  //itkFloatIterator normalizedCentralPatchIt(normalizedCentralPatch, normalizedCentralPatch->GetRequestedRegion());


  //set the search region around the current pixel
  typename itkTImage::RegionType searchRegion;
  ComputeSearchRegion(p,searchRegion);
  
  //create the patch for pixels in the neighbourhood of the current pixel
  itkTPointer neighbourPatch = itkTImage::New();
  //itkFloatPointer neighbourPatch = itkFloatImage::New();
  CreatePatch(neighbourPatch);
  //itkFloatIterator neighbourPatchIt(neighbourPatch, neighbourPatch->GetRequestedRegion());
  itkTIterator neighbourPatchIt(neighbourPatch, neighbourPatch->GetRequestedRegion());
  
  //go through the neighbourhood with a region iterator
  itkTIteratorWithIndex it( m_inputImage, searchRegion);
  typename itkTImage::IndexType neighbourPixelIndex;	

  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
    neighbourPixelIndex = it.GetIndex();

    for(unsigned int i=0; i < m_anatomicalImages.size(); i++){

      bool goForIt = true;
      if(m_optimized == 1)
        goForIt = CheckSpeed(p, neighbourPixelIndex, i);

      if(goForIt == true){
      
        //float meanNeighbour = m_meanAnatomicalImages[i]->GetPixel(neighbourPixelIndex);
        //float stddevNeighbour = sqrt( m_varianceAnatomicalImages[i]->GetPixel(neighbourPixelIndex) );

        GetPatch(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i]);
        //GetNormalizedPatch(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i], meanNeighbour, stddevNeighbour);

        double weight = exp( - PatchDistance(centralPatch, neighbourPatch) / m_rangeBandwidth);
        //double weight = exp( - PatchDistance(normalizedCentralPatch, neighbourPatch) / m_rangeBandwidth);

        if(weight>wmax)
          if( (p[0] != neighbourPixelIndex[0]) && (p[1] != neighbourPixelIndex[1]) && (p[2] != neighbourPixelIndex[2]) ) //has to be modify (not appropriate if the input image is in the anatomical input images.
	    wmax = weight;
      
        sum += weight;
    
        //Add this patch to the current estimate using the computed weight
        for(patchIt.GoToBegin(), neighbourPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++neighbourPatchIt){
          patchIt.Set( patchIt.Get() + neighbourPatchIt.Get() * weight );
          //patchIt.Set( patchIt.Get() + (neighbourPatchIt.Get() * stddev + mean ) * weight );
        }
      }
    }
  }

  
  //consider now the special case of the central patch
  switch(m_centralPointStrategy){
    case 0:                                        //remove the central patch to the estimated patch
      for(patchIt.GoToBegin(), centralPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++centralPatchIt)
        patchIt.Set( patchIt.Get() -1.0 * centralPatchIt.Get() );
      sum -= 1.0;
      break;
    case 1: break;                                 //nothing to do
    case -1:
      for(patchIt.GoToBegin(), centralPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++centralPatchIt)
        patchIt.Set( patchIt.Get()  + (wmax -1.0) * centralPatchIt.Get() );
      sum += (wmax - 1.0);
      break;
    default:                                       //as in case -1
      for(patchIt.GoToBegin(), centralPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++centralPatchIt)
        patchIt.Set( patchIt.Get()  + (wmax -1.0) * centralPatchIt.Get() );
      sum += (wmax - 1.0);
      break;
  }  
    
  if(sum>0.0001)
    //Normalization of the denoised patch 
    for(patchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt)
      patchIt.Set( patchIt.Get() / sum );
  else
    //copy the central patch to the denoised patch
    for(patchIt.GoToBegin(), centralPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++centralPatchIt)
      patchIt.Set( centralPatchIt.Get() );  
    
  return sum;
}

template <typename T>
double LabelFusionTool<T>::GetHRPatch(typename itkTImage::IndexType p, itkFloatPointer & patch)
{
  double wmax = 0; //maximum weight of patches
  double sum  = 0; //sum of weights (used for normalization purpose)
  
  //create the patch and set the estimate to 0
  CreatePatch(patch);
  itkFloatIterator patchIt(patch, patch->GetRequestedRegion());

  //get the value of the patch around the current pixel
  itkTPointer centralPatch = itkTImage::New();
  CreatePatch(centralPatch);
  GetPatch(p, centralPatch, m_inputImage);
  itkTIterator centralPatchIt(centralPatch, centralPatch->GetRequestedRegion());  
  
  //set the search region around the current pixel
  typename itkTImage::RegionType searchRegion;
  ComputeSearchRegion(p,searchRegion);
  
  //create the patch for pixels in the neighbourhood of the current pixel in the anatomical images
  itkTPointer neighbourPatch = itkTImage::New();
  CreatePatch(neighbourPatch);
  itkTIterator neighbourPatchIt(neighbourPatch, neighbourPatch->GetRequestedRegion());

  //create the patch for pixels in the neighbourhood of the current pixel in the HR images
  itkTPointer neighbourHRPatch = itkTImage::New();
  CreatePatch(neighbourHRPatch);
  itkTIterator neighbourHRPatchIt(neighbourHRPatch, neighbourHRPatch->GetRequestedRegion());
    
  //go through the neighbourhood with a region iterator
  itkTIteratorWithIndex it( m_inputImage, searchRegion);
  typename itkTImage::IndexType neighbourPixelIndex;	
  
  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
    neighbourPixelIndex = it.GetIndex();
    
    for(unsigned int i=0; i < m_anatomicalImages.size(); i++){
      
      bool goForIt = true;
      if(m_optimized == 1)
        goForIt = CheckSpeed(p, neighbourPixelIndex, i);
      
      if(goForIt == true){
        
        GetPatch(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i]);
        GetPatch(neighbourPixelIndex, neighbourHRPatch, m_labelImages[i]);
        
        double weight = exp( - PatchDistance(centralPatch, neighbourPatch) / m_rangeBandwidth);
        
        if(weight>wmax)
          if( (p[0] != neighbourPixelIndex[0]) && (p[1] != neighbourPixelIndex[1]) && (p[2] != neighbourPixelIndex[2]) ) //has to be modify (not appropriate if the input image is in the anatomical input images.
            wmax = weight;
        
        sum += weight;
        
        //Add this patch to the current estimate using the computed weight
        for(patchIt.GoToBegin(), neighbourHRPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++neighbourHRPatchIt){
          patchIt.Set( patchIt.Get() + neighbourHRPatchIt.Get() * weight );
        }
      }
    }
  }
  
  
  //In this case, no consideration of the central patch
  
  if(sum>0.0001)
    //Normalization of the denoised patch 
    for(patchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt){
      patchIt.Set( patchIt.Get() / sum );
    }
  else{
    switch(m_centralPointStrategy){
      case 0: 
        sum = 0;
        for(patchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt)
          patchIt.Set( 0 ); 
        break;
      case 1: 
        //copy the central patch to the current patch
        for(patchIt.GoToBegin(), centralPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++centralPatchIt)
          patchIt.Set( centralPatchIt.Get() ); 
        break;
      default: 
        //copy the central patch to the current patch
        for(patchIt.GoToBegin(), centralPatchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++centralPatchIt)
          patchIt.Set( centralPatchIt.Get() ); 
        break;                
    }
  }  
  
  return sum;
}


template <typename T>
void LabelFusionTool<T>::AddPatchToImage(typename itkTImage::IndexType p, itkFloatPointer & patch, itkFloatPointer & image, itkFloatPointer & weightImage, double weight)
{
  //this function add a patch value to the denoised image
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkFloatIterator imageIt( image, imageRegion);
  itkFloatIterator weightIt( weightImage, imageRegion);
  itkFloatIterator patchIt( patch, patchRegion);

  for ( imageIt.GoToBegin(), patchIt.GoToBegin(), weightIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++patchIt, ++weightIt){
    double previous_weight = weightIt.Get();
    double new_weight = previous_weight + weight;
    //trick to avoid limits overflow
    float new_value = imageIt.Get();
    if(new_weight > 0)
      new_value = imageIt.Get() * (previous_weight / new_weight) + patchIt.Get() * (weight / new_weight);
    
    //imageIt.Set(new_value);
    //weightIt.Set(new_weight);
    imageIt.Set( imageIt.Get() + weight * patchIt.Get() );
    weightIt.Set( weightIt.Get() + weight );
  }
  
}

template <typename T>
double LabelFusionTool<T>:: GetLabelPatch(typename itkTImage::IndexType p, itkTPointer & patch)
{
  double sum  = 0; //sum of weights (used for normalization purpose)
  double wmax = 0;

  //create the patch and set the estimate to 0
  CreatePatch(patch);
  itkTIterator patchIt(patch, patch->GetRequestedRegion());

  //get the (intensity) value of the patch around the current pixel
  itkTPointer centralPatch = itkTImage::New();
  CreatePatch(centralPatch);
  GetPatch(p, centralPatch, m_inputImage);

  //set the search region around the current pixel
  typename itkTImage::RegionType searchRegion;
  ComputeSearchRegion(p,searchRegion);

  //create the (intensity) patch for pixels in the neighbourhood of the current pixel
  itkTPointer neighbourPatch = itkTImage::New();
  CreatePatch(neighbourPatch);

  //create the (label) patch for pixels in the neighbourhood of the current pixel
  itkTPointer neighbourLabelPatch = itkTImage::New();
  CreatePatch(neighbourLabelPatch);
  itkTIterator neighbourPatchIt(neighbourLabelPatch, neighbourLabelPatch->GetLargestPossibleRegion());

  //vector (1D patch) to store the cumulative weight for each label
  std::vector< std::map<T, float> > vectorMap( patch->GetLargestPossibleRegion().GetNumberOfPixels() );
  typename std::map<T, float>::iterator mapIt;
  unsigned int k = 0;

  //go through the neighbourhood with a region iterator
  itkTIteratorWithIndex it( m_inputImage, searchRegion);
  typename itkTImage::IndexType neighbourPixelIndex;

  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
    neighbourPixelIndex = it.GetIndex();

    for(unsigned int i=0; i < m_anatomicalImages.size(); i++){

      bool goForIt = true;
      if(m_optimized == 1)
        goForIt = CheckSpeed(p, neighbourPixelIndex, i);

      if(goForIt == true){

        Get2Patches(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i], neighbourLabelPatch, m_labelImages[i]); 
        //GetPatch(neighbourPixelIndex, neighbourLabelPatch, m_labelImages[i]);
        //GetPatch(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i]);

        double weight = exp( - PatchDistance(centralPatch, neighbourPatch) / m_rangeBandwidth);
        sum += weight;

        //Add this label patch to the current estimate using the computed weight
        k = 0;
        for(neighbourPatchIt.GoToBegin(); !neighbourPatchIt.IsAtEnd(); ++neighbourPatchIt, ++k){
          vectorMap[k][neighbourPatchIt.Get()] += weight;
        }
      }
    }
  }

  k = 0;
  for(patchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++k){
    short label = -1;
    wmax = 0;
    for(mapIt = vectorMap[k].begin (); mapIt != vectorMap[k].end (); ++mapIt){
      if( (*mapIt).second > wmax ){
        wmax = (*mapIt).second;
        label = (*mapIt).first;
      }
    }
    patchIt.Set( label );
  }


  return wmax;
}

template <typename T>
void LabelFusionTool<T>:: GetFuzzyLabelPatch(typename itkTImage::IndexType p, itkTPointer & patch, std::vector< std::map<T, float> > & vectorMap)
{
  //create the patch and set the estimate to 0
  CreatePatch(patch);
  itkTIterator patchIt(patch, patch->GetRequestedRegion());

  //get the (intensity) value of the patch around the current pixel
  itkTPointer centralPatch = itkTImage::New();
  CreatePatch(centralPatch);
  GetPatch(p, centralPatch, m_inputImage);

  //set the search region around the current pixel
  typename itkTImage::RegionType searchRegion;
  ComputeSearchRegion(p,searchRegion);

  //create the (intensity) patch for pixels in the neighbourhood of the current pixel
  itkTPointer neighbourPatch = itkTImage::New();
  CreatePatch(neighbourPatch);

  //create the (label) patch for pixels in the neighbourhood of the current pixel
  itkTPointer neighbourLabelPatch = itkTImage::New();
  CreatePatch(neighbourLabelPatch);
  itkTIterator neighbourPatchIt(neighbourLabelPatch, neighbourLabelPatch->GetLargestPossibleRegion());

  vectorMap.resize( patch->GetLargestPossibleRegion().GetNumberOfPixels() );
  typename std::map<T, float>::iterator mapIt;
  unsigned int k = 0;

  //go through the neighbourhood with a region iterator
  itkTIteratorWithIndex it( m_inputImage, searchRegion);
  typename itkTImage::IndexType neighbourPixelIndex;

  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
    neighbourPixelIndex = it.GetIndex();

    for(unsigned int i=0; i < m_anatomicalImages.size(); i++){

      bool goForIt = true;
      if(m_optimized == 1)
        goForIt = CheckSpeed(p, neighbourPixelIndex, i);

      if(goForIt == true){

        Get2Patches(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i], neighbourLabelPatch, m_labelImages[i]); 
        //GetPatch(neighbourPixelIndex, neighbourLabelPatch, m_labelImages[i]);
        //GetPatch(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i]);
        double weight = exp( - PatchDistance(centralPatch, neighbourPatch) / m_rangeBandwidth);

        //Add this label patch to the current estimate using the computed weight
        k = 0;
        for(neighbourPatchIt.GoToBegin(); !neighbourPatchIt.IsAtEnd(); ++neighbourPatchIt, ++k){
          vectorMap[k][neighbourPatchIt.Get()] += weight;
        }
      }
    }
  }
}


template <typename T>
double LabelFusionTool<T>:: GetLabelPatchUsingNormalization(typename itkTImage::IndexType p, itkTPointer & patch)
{
  double sum  = 0; //sum of weights (used for normalization purpose)
  double wmax = 0;

  //create the patch and set the estimate to 0
  CreatePatch(patch);
  itkTIterator patchIt(patch, patch->GetRequestedRegion());

  //get the (intensity) value of the patch around the current pixel
  itkFloatPointer normalizedCentralPatch = itkFloatImage::New();
  CreatePatch(normalizedCentralPatch);
  float mean = m_meanImage->GetPixel(p);
  float stddev = sqrt( m_varianceImage->GetPixel(p) );
  GetNormalizedPatch(p, normalizedCentralPatch, m_inputImage, mean, stddev);

  //set the search region around the current pixel
  typename itkTImage::RegionType searchRegion;
  ComputeSearchRegion(p,searchRegion);

  //create the (intensity) patch for pixels in the neighbourhood of the current pixel
  itkFloatPointer neighbourPatch = itkFloatImage::New();
  CreatePatch(neighbourPatch);

  //create the (label) patch for pixels in the neighbourhood of the current pixel
  itkTPointer neighbourLabelPatch = itkTImage::New();
  CreatePatch(neighbourLabelPatch);
  itkTIterator neighbourPatchIt(neighbourLabelPatch, neighbourLabelPatch->GetLargestPossibleRegion());

  //vector (1D patch) to store the cumulative weight for each label
  std::vector< std::map<T, float> > vectorMap( patch->GetLargestPossibleRegion().GetNumberOfPixels() );
  typename std::map<T, float>::iterator mapIt;
  unsigned int k = 0;

  //go through the neighbourhood with a region iterator
  itkTIteratorWithIndex it( m_inputImage, searchRegion);
  typename itkTImage::IndexType neighbourPixelIndex;

  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
    neighbourPixelIndex = it.GetIndex();

    for(unsigned int i=0; i < m_anatomicalImages.size(); i++){

      bool goForIt = true;
      if(m_optimized == 1)
        goForIt = CheckSpeed(p, neighbourPixelIndex, i);

      if(goForIt == true){

        float meanNeighbour = m_meanAnatomicalImages[i]->GetPixel(neighbourPixelIndex);
        float stddevNeighbour = sqrt( m_varianceAnatomicalImages[i]->GetPixel(neighbourPixelIndex) );

        GetNormalizedPatch(neighbourPixelIndex, neighbourPatch, m_anatomicalImages[i], meanNeighbour, stddevNeighbour);
        GetPatch(neighbourPixelIndex, neighbourLabelPatch, m_labelImages[i]);

        double weight = exp( - PatchDistance(normalizedCentralPatch, neighbourPatch) / m_rangeBandwidth);
        sum += weight;

        //Add this label patch to the current estimate using the computed weight
        k = 0;
        for(neighbourPatchIt.GoToBegin(); !neighbourPatchIt.IsAtEnd(); ++neighbourPatchIt, ++k){
          vectorMap[k][neighbourPatchIt.Get()] += weight;
        }
      }
    }
  }

  k = 0;
  for(patchIt.GoToBegin(); !patchIt.IsAtEnd(); ++patchIt, ++k){
    short label = -1;
    wmax = 0;
    for(mapIt = vectorMap[k].begin (); mapIt != vectorMap[k].end (); ++mapIt){
      if( (*mapIt).second > wmax ){
        wmax = (*mapIt).second;
        label = (*mapIt).first;
      }
    }
    patchIt.Set( label );
  }


  return wmax;
}

template <typename T>
void LabelFusionTool<T>::AddLabelPatchToLabelImage(typename itkTImage::IndexType p, itkTPointer & patch, double weight)
{
  //this function add a label patch value to the cumulative label image
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkMapIterator imageIt( m_labelFusionImage, imageRegion);
  itkTIterator patchIt( patch, patchRegion);
  std::map<T, float> map;
  typename std::map<T, float>::iterator mapIt;

  for ( imageIt.GoToBegin(), patchIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++patchIt){
    map = imageIt.Get();
    map[patchIt.Get()] += weight;
    imageIt.Set( map );
  }
  
}

template <typename T>
void LabelFusionTool<T>::AddFuzzyLabelPatchToLabelImage(typename itkTImage::IndexType p, itkTPointer & patch, std::vector< std::map<T, float> > & vectorMap)
{
  //this function add a fuzzy label patch value to the cumulative label image
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkMapIterator imageIt( m_labelFusionImage, imageRegion);
  itkTIterator patchIt( patch, patchRegion);
  std::map<T, float> map;
  typename std::map<T, float>::iterator mapIt;

  unsigned int k = 0;

  for ( imageIt.GoToBegin(), patchIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++patchIt, ++k){
    map = imageIt.Get();
    for(mapIt = vectorMap[k].begin (); mapIt != vectorMap[k].end (); ++mapIt)
      map[(*mapIt).first] += (*mapIt).second;

    imageIt.Set( map );
  }
  
}


template <typename T>
double LabelFusionTool<T>::PatchDistance(itkTPointer & p,itkTPointer & q)
{
  double diff=0;
  double dist = 0;
  itkTConstIterator itp( p, p->GetLargestPossibleRegion() );
  itkTConstIterator itq( q, q->GetLargestPossibleRegion() );

  for(itp.GoToBegin(), itq.GoToBegin(); !itp.IsAtEnd(); ++itp, ++itq){
    diff = itp.Get() - itq.Get();
    dist += diff*diff;	
  }
  return dist;
}

template <typename T>
double LabelFusionTool<T>::PatchDistance(itkFloatPointer & p,itkFloatPointer & q)
{
  double diff=0;
  double dist = 0;
  itkFloatConstIterator itp( p, p->GetLargestPossibleRegion() );
  itkFloatConstIterator itq( q, q->GetLargestPossibleRegion() );

  for(itp.GoToBegin(), itq.GoToBegin(); !itp.IsAtEnd(); ++itp, ++itq){
    diff = itp.Get() - itq.Get();
    dist += diff*diff;	
  }
  return dist;
}

template <typename T>
void LabelFusionTool<T>::GetOutput(itkTPointer & outputImage)
{
  outputImage = m_outputImage; 
}

template <typename T>
void LabelFusionTool<T>::GetWeightImage(itkTPointer & outputImage)
{
  typename itk::RescaleIntensityImageFilter<itkFloatImage,itkTImage>::Pointer rescaleFilter = itk::RescaleIntensityImageFilter<itkFloatImage,itkTImage>::New();
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(10000);
  rescaleFilter->SetInput(m_weightImage);
  rescaleFilter->UpdateLargestPossibleRegion();
    
  outputImage = rescaleFilter->GetOutput();
}

template <typename T>
void LabelFusionTool<T>::GetFuzzyWeightImage(itkFloatPointer & outputImage, int label)
{
  itkMapIterator labelImageIt( m_labelFusionImage, m_labelFusionImage->GetLargestPossibleRegion());
  itkFloatIterator outputImageIt( outputImage, outputImage->GetLargestPossibleRegion());

  for(labelImageIt.GoToBegin(), outputImageIt.GoToBegin(); !labelImageIt.IsAtEnd(); ++labelImageIt, ++outputImageIt){
    std::map<T, float> map;
    typename std::map<T, float>::iterator mapIt;

    map = labelImageIt.Get();

    for(mapIt = map.begin (); mapIt != map.end (); ++mapIt){
			if( (*mapIt).first == label ){
				outputImageIt.Set( outputImageIt.Get() + (*mapIt).second );
			}
		}
  }

}


template <typename T>
void LabelFusionTool<T>::ReadOneImage(std::string input_file, itkTPointer & image)
{
  typename itkTReader::Pointer reader = itkTReader::New();
  reader->SetFileName( input_file );
  reader->Update();
  image = reader->GetOutput();

  typename itkTImage::RegionType region;
  typename itkTImage::SizeType size;
  typename itkTImage::SpacingType spacing;

  //Print image information (size and spacing)
  region = image->GetLargestPossibleRegion();
  size = region.GetSize();
  std::cout<<"\n ************** debug ******************* \n Image size is : "<<size[0]<<" "<<size[1]<<" "<<size[2]<<"\n";
  spacing = image->GetSpacing();
  //std::cout<<"Spacing of inputImage (in mm) : "<<spacing[0]<<", "<<spacing[1]<<", "<<spacing[2]<<"\n";
  //std::cout<<"------------------------------\n";}

}

#endif

