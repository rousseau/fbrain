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

#include "btkPatch2.h"


namespace btk
{

template<typename T>
void Patch2<T>::CreateEmptyPatch(itkTImagePointer & image, typename itkTImage::SizeType & fullPatchSize)
{
  //purpose : resize and allocate the patch and set the estimate to 0
  m_Data = itkTImage::New();

  //Compute patch size
  typename itkTImage::RegionType  fullPatchRegion;        //ITK region corresponding to a patch
  fullPatchRegion.SetSize(fullPatchSize);

  m_Data->SetRegions(fullPatchRegion);
  m_Data->Allocate();

  m_Data->FillBuffer(0);
}

template<typename T>
void Patch2<T>::ComputePatch(itkTImagePointer & image, typename itkTImage::IndexType & p)
{
  //Set the center point of the patch in the image coordinate system
  SetCentralPointInImage(p);

  //Set the patch values with the corresponding image values
  typename itkTImage::RegionType imageRegion;
  typename itkTImage::RegionType patchRegion;
  ComputePatchRegion(image,p, imageRegion,patchRegion);

  itkConstIterator  inputIt( image, imageRegion);
  itkTIterator      outputIt(m_Data, patchRegion);

  for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt, ++outputIt)
    outputIt.Set( inputIt.Get() );
}

template<typename T>
void Patch2<T>::ComputePatchRegion(itkTImagePointer & image, typename itkTImage::IndexType & p, typename itkTImage::RegionType & imageRegion, typename itkTImage::RegionType & patchRegion)
{
  //create an appropriate patch region around the current pixel p, and consider an offset for pixels close to boundaries.
  //the patch region is defined in the image coordinate system.
  //the offset allows to do the link between the patch coordinate system and the one of the image.
  typename itkTImage::RegionType::IndexType start;
  typename itkTImage::RegionType::IndexType offset;
  typename itkTImage::RegionType::SizeType size;

  typename itkTImage::SizeType    fullPatchSize = m_Data->GetLargestPossibleRegion().GetSize();          //patch size  : 2 * halfPatchSize + 1
  typename itkTImage::SizeType    imageSize     = image->GetLargestPossibleRegion().GetSize();

  for(unsigned int i=0; i!= p.GetIndexDimension(); i++){
    start[i] = p[i] - (fullPatchSize[i]-1)/2;
    size[i] = fullPatchSize[i];
    offset[i] = 0;

    if(start[i] < 0){                  //if the starting index is outside the image (<0)
      size[i] += start[i];
      offset[i] = -start[i];
      start[i] = 0;
    }
    else if(start[i] >= (int)imageSize[i]){    //if the starting index is outside the image (>image size)
      start[i] = imageSize[i]-1;
      size[i] = 0;
    }

    int d = (start[i] + size[i]) - imageSize[i]; //if the region is not fully inside the image
    if(d>0){
      if(static_cast< int >(size[i]) - d < 0)
      {
          size[i] = 0;
      }
      else
        size[i] = size[i] - d;

    }
  }
  imageRegion.SetSize( size );
  imageRegion.SetIndex( start );
  patchRegion.SetSize( size );
  patchRegion.SetIndex( offset );
}

template<typename T>
void Patch2<T>::ComputeMeanAndVariance(float & mean, float & variance)
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
  mean   = m / n;
  variance = (m2 / n) - (mean * mean) ;
}

template<typename T>
T Patch2<T>::GetDataAtCenter()
{
  typename itkTImage::IndexType p;
  for(unsigned int i=0; i!= p.GetIndexDimension(); i++)
   p[i] = (m_Data->GetLargestPossibleRegion().GetSize()[i]-1)/2;

  return m_Data->GetPixel(p);
}

} // namespace btk

