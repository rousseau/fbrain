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

#include "btkPatch.h"


namespace btk
{

template<typename T>
void Patch<T>::Initialize(itkTImagePointer & image, typename itkTImage::SizeType h)
{
  SetInputImage(image);
  SetHalfPatchSize(h);
  
  m_FullPatchSize[0] = 2 * m_HalfPatchSize[0] + 1;
  m_FullPatchSize[1] = 2 * m_HalfPatchSize[1] + 1;
  m_FullPatchSize[2] = 2 * m_HalfPatchSize[2] + 1;

  m_FullPatchRegion.SetSize(m_FullPatchSize);
  
  m_CentralPoint[0] = m_HalfPatchSize[0];
  m_CentralPoint[1] = m_HalfPatchSize[1];
  m_CentralPoint[2] = m_HalfPatchSize[2];
    
  CreatePatch();  
}


template<typename T>
void Patch<T>::Initialize(itkTImagePointer & image, int h)
{
  SetInputImage(image);
  SetPatchSize(h);
  CreatePatch();  
}

template<typename T>
void Patch<T>::ComputePatch(typename itkTImage::IndexType p, itkTImagePointer & image)
{
  //Set the center point of the patch in the image coordinate system
  SetCentralPointInImage(p);
  
  //Set the patch values with the corresponding image values
  m_Data->FillBuffer(0);
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);
  
  itkConstIterator  inputIt( image, imageRegion);
  itkTIterator      outputIt(m_Data, patchRegion);

  for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
    outputIt.Set( inputIt.Get() ); 
}



template<typename T>
void Patch<T>::SetInputImage(itkTImagePointer & image)
{
  //compute characteristics of the input image
  m_ImageRegion  = image->GetLargestPossibleRegion();
  m_ImageSize    = m_ImageRegion.GetSize();
  m_ImageSpacing = image->GetSpacing();  
}

template<typename T>
void Patch<T>::SetPatchSize(int h)
{
  //std::cout<<"Compute patch size (taking into account possible image anisotropy)\n";    
  float minVoxSz = m_ImageSpacing[0];
  if(m_ImageSpacing[1] < minVoxSz) minVoxSz = m_ImageSpacing[1];
  if(m_ImageSpacing[2] < minVoxSz) minVoxSz = m_ImageSpacing[2];

  m_HalfPatchSize[0] = (int)(0.5 + h * minVoxSz / m_ImageSpacing[0]);
  m_HalfPatchSize[1] = (int)(0.5 + h * minVoxSz / m_ImageSpacing[1]);
  m_HalfPatchSize[2] = (int)(0.5 + h * minVoxSz / m_ImageSpacing[2]);

  m_FullPatchSize[0] = 2 * m_HalfPatchSize[0] + 1;
  m_FullPatchSize[1] = 2 * m_HalfPatchSize[1] + 1;
  m_FullPatchSize[2] = 2 * m_HalfPatchSize[2] + 1;

  m_FullPatchRegion.SetSize(m_FullPatchSize);
  
  m_CentralPoint[0] = m_HalfPatchSize[0];
  m_CentralPoint[1] = m_HalfPatchSize[1];
  m_CentralPoint[2] = m_HalfPatchSize[2];
}

template<typename T>
void Patch<T>::CreatePatch()
{
  m_Data = itkTImage::New();
  //resize and allocate the patch and set the estimate to 0
  m_Data->SetRegions(m_FullPatchRegion);
  m_Data->Allocate();
  m_Data->FillBuffer(0);  
}


template<typename T>
void Patch<T>::ComputePatchRegion(typename itkTImage::IndexType p, typename itkTImage::RegionType & imageRegion, typename itkTImage::RegionType & patchRegion)
{
  //create an appropriate patch region around the current pixel p, and consider an offset for pixels close to boundaries.
  //the patch region is defined in the image coordinate system.
  //the offset allows to do the link between the patch coordinate system and the one of the image.
  typename itkTImage::RegionType::IndexType start;
  typename itkTImage::RegionType::IndexType offset;
  typename itkTImage::RegionType::SizeType size;

  for(unsigned int i=0; i!= p.GetIndexDimension(); i++){
    start[i] = p[i] - m_HalfPatchSize[i];
    size[i] = m_FullPatchSize[i];
    offset[i] = 0;
    
    if(start[i] < 0){                  //if the starting index is outside the image (<0)
      size[i] += start[i];
      offset[i] = -start[i];
      start[i] = 0;
    }
    else if(start[i] >= (int)m_ImageSize[i]){    //if the starting index is outside the image (>image size)
      start[i] = m_ImageSize[i]-1;
      size[i] = 0;
    }
    
    int d = (start[i] + size[i]) - m_ImageSize[i]; //if the region is not fully inside the image
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

template<typename T>
T Patch<T>::GetCentralValue()
{
  return m_Data->GetPixel(m_CentralPoint); 
}

template<typename T>
void Patch<T>::ComputeMeanAndStdDevValues()
{
  //Possible enhancement: take into account a mask image
  double m = 0;
  double m2= 0;
  itkConstIterator itp( m_Data, m_Data->GetRequestedRegion() );

  for(itp.GoToBegin(); !itp.IsAtEnd(); ++itp)
  {
    m += itp.Get();
    m2+= (itp.Get() * itp.Get());
  }
  int n = m_Data->GetLargestPossibleRegion().GetNumberOfPixels();
  m_MeanValue   = m / n;
  m_StdDevValue = (m2 / n) - (m_MeanValue * m_MeanValue) ;                        
}



} // namespace btk

