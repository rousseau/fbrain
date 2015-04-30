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

#include "btkPatchTool2.h"


namespace btk
{


template<typename T>
void PatchTool2<T>::ComputePatchSize(itkTImagePointer & image, int & h)
{
  std::cout<<"Compute patch size (taking into account possible image anisotropy)\n";
  typename itkTImage::SpacingType imageSpacing = image->GetSpacing();
  float minVoxSz = imageSpacing[0];
  if(imageSpacing[1] < minVoxSz) minVoxSz = imageSpacing[1];
  if(imageSpacing[2] < minVoxSz) minVoxSz = imageSpacing[2];

  m_HalfPatchSize[0] = (int)(0.5 + h * minVoxSz / imageSpacing[0]);
  m_HalfPatchSize[1] = (int)(0.5 + h * minVoxSz / imageSpacing[1]);
  m_HalfPatchSize[2] = (int)(0.5 + h * minVoxSz / imageSpacing[2]);

  m_FullPatchSize[0] = 2 * m_HalfPatchSize[0] + 1;
  m_FullPatchSize[1] = 2 * m_HalfPatchSize[1] + 1;
  m_FullPatchSize[2] = 2 * m_HalfPatchSize[2] + 1;
  std::cout<<"half patch size : "<<m_HalfPatchSize[0]<<" "<<m_HalfPatchSize[1]<<" "<<m_HalfPatchSize[2]<<std::endl;
  std::cout<<"full patch size : "<<m_FullPatchSize[0]<<" "<<m_FullPatchSize[1]<<" "<<m_FullPatchSize[2]<<std::endl;

  m_ImageSize = image->GetLargestPossibleRegion().GetSize();
  std::cout<<"Input Image Size : "<<m_ImageSize[0]<<" "<<m_ImageSize[1]<<" "<<m_ImageSize[2]<<std::endl;
}

template<typename T>
void PatchTool2<T>::ComputeSpatialBandwidth(itkTImagePointer & image, int & h)
{
  std::cout<<"Computing spatial bandwidth (taking into account possible image anisotropy)\n";
  typename itkTImage::SpacingType imageSpacing = image->GetSpacing();
  float minVoxSz = imageSpacing[0];
  if(imageSpacing[1] < minVoxSz) minVoxSz = imageSpacing[1];
  if(imageSpacing[2] < minVoxSz) minVoxSz = imageSpacing[2];

  m_HalfSpatialBandwidth[0] = (int)(0.5 + h * minVoxSz / imageSpacing[0]);
  m_HalfSpatialBandwidth[1] = (int)(0.5 + h * minVoxSz / imageSpacing[1]);
  m_HalfSpatialBandwidth[2] = (int)(0.5 + h * minVoxSz / imageSpacing[2]);
  std::cout<<"half spatialBandwidth : "<<m_HalfSpatialBandwidth[0]<<" "<<m_HalfSpatialBandwidth[1]<<" "<<m_HalfSpatialBandwidth[2]<<"\n";

  m_FullSpatialBandwidth[0] = 2 * m_HalfSpatialBandwidth[0] + 1;
  m_FullSpatialBandwidth[1] = 2 * m_HalfSpatialBandwidth[1] + 1;
  m_FullSpatialBandwidth[2] = 2 * m_HalfSpatialBandwidth[2] + 1;
}

template<typename T>
void PatchTool2<T>::ComputeSearchRegion(typename itkTImage::IndexType & point, itkTImagePointer & image, typename itkTImage::RegionType & region)
{
  //create an appropiate search region around the current pixel p
  typename itkTImage::RegionType::IndexType start;
  typename itkTImage::RegionType::SizeType size;

  typename itkTImage::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

  for(unsigned int i=0; i!= point.GetIndexDimension(); i++){
    start[i] = point[i] - m_HalfSpatialBandwidth[i];
    size[i] = m_FullSpatialBandwidth[i];

    if(start[i] < 0){                  //if the starting index is outside the image (<0)
      size[i] += start[i];
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
  region.SetSize( size );
  region.SetIndex( start );

}

template<typename T>
double PatchTool2<T>::ComputeL2NormBetweenPatches(Patch2<T> & p, Patch2<T> & q)
{
  double diff=0;
  double dist = 0;
  itkConstTIterator itp( p.GetData(), p.GetData()->GetLargestPossibleRegion() );
  itkConstTIterator itq( q.GetData(), q.GetData()->GetLargestPossibleRegion() );

  for(itp.GoToBegin(), itq.GoToBegin(); !itp.IsAtEnd(); ++itp, ++itq){
    diff = itp.Get() - itq.Get();
    dist += diff*diff;
  }
  return dist;
}

template<typename T>
void PatchTool2<T>::ComputeMeanAndVarianceImage(itkTImagePointer & image, itkTImagePointer & maskImage, itkFloatImagePointer & meanImage, itkFloatImagePointer & varianceImage)
{
  std::cout<<"Computing mean and variance images ...\n";
  int x,y,z;
  typename itkTImage::SizeType imageSize = image->GetLargestPossibleRegion().GetSize();

  #pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)imageSize[2]; z++)
    for(y=0; y < (int)imageSize[1]; y++)
      for(x=0; x < (int)imageSize[0]; x++){
         typename itkTImage::IndexType p;
         p[0] = x;
         p[1] = y;
         p[2] = z;

         if( maskImage->GetPixel(p) > 0 ){
           btk::Patch2<T> patch = btk::Patch2<T>(image, m_FullPatchSize, p);
           float mean = 0;
           float variance = 0;
           patch.ComputeMeanAndVariance(mean, variance);

           meanImage->SetPixel( p, mean );
           varianceImage->SetPixel( p, variance );
         }
      }
}

template<typename T>
void PatchTool2<T>::GetNeighboursUsingMeanAndVariance(typename itkTImage::IndexType & point, itkTImagePointer & image, float & mean, float & variance, itkFloatImagePointer & meanImage, itkFloatImagePointer & varianceImage, std::vector< typename itkTImage::IndexType > & neighbours)
{
    neighbours.clear();

    //Set the search region for neighbours
    typename itkTImage::RegionType searchRegion;
    this->ComputeSearchRegion(point, image, searchRegion);

    //Go through the neighbourhood with a region iterator
    itkTIteratorWithIndex it(image, searchRegion);
    itkFloatIterator itMean(meanImage, searchRegion);
    itkFloatIterator itVariance(varianceImage, searchRegion);
    typename itkTIteratorWithIndex::IndexType neighbourPixelIndex;

    for(it.GoToBegin(), itMean.GoToBegin(), itVariance.GoToBegin(); !it.IsAtEnd(); ++it, ++itMean, ++itVariance){

      neighbourPixelIndex = it.GetIndex();

      double meanRatio = 0;
      double varianceRatio = 0;

      float meanNeighbour     = itMean.Get();
      float varianceNeighbour = itVariance.Get();

      if( meanNeighbour == 0 )
      {
        if(mean == 0)
            meanRatio = 1;
      }
      else
        meanRatio = mean / meanNeighbour;

      if( varianceNeighbour == 0 )
      {
        if(variance == 0)
            varianceRatio = 1;
      }
      else
        varianceRatio = variance / varianceNeighbour;

      if( (meanRatio > 0.95) && (meanRatio < 1/0.95) && (varianceRatio > 0.5) && (varianceRatio < 2) )
      {
          neighbours.push_back(neighbourPixelIndex);
      }
    }
}

template<typename T>
void PatchTool2<T>::GetNeighbours(typename itkTImage::IndexType & point, itkTImagePointer & image, std::vector< typename itkTImage::IndexType > & neighbours)
{
  neighbours.clear();

  //Set the search region for neighbours
  typename itkTImage::RegionType searchRegion;
  this->ComputeSearchRegion(point, image, searchRegion);

  //Go through the neighbourhood with a region iterator
  itkTIteratorWithIndex it(image, searchRegion);
  typename itkTIteratorWithIndex::IndexType neighbourPixelIndex;

  for(it.GoToBegin(); !it.IsAtEnd(); ++it){
    neighbours.push_back( it.GetIndex() );
  }
}

template<typename T>
void PatchTool2<T>::RemoveCentralPointInNeighborhood(typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours)
{
  for (unsigned int i = 0; i < neighbours.size(); i++)
    if( (neighbours[i][0] == point[0]) && (neighbours[i][1] == point[1]) && (neighbours[i][2] == point[2]) )
    {
      neighbours.erase(neighbours.begin()+i);
      break;
    }
}

template<typename T>
void PatchTool2<T>::ComputeNeighbourWeights(itkTImagePointer & image, typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours, std::vector<double> & weights, double & sumOfWeights, double & smoothing)
{
  btk::Patch2<T> patch = btk::Patch2<T>(image, m_FullPatchSize, point);
  btk::Patch2<T> neighbourPatch = btk::Patch2<T>(image, m_FullPatchSize);

  sumOfWeights = 0;
  double neighbourWeight = 0;

  for(unsigned int n = 0; n < neighbours.size(); n++)
  {
    neighbourPatch.ComputePatch(image, neighbours[n]);
    neighbourWeight = exp( - ComputeL2NormBetweenPatches(patch, neighbourPatch ) / smoothing);
    weights.push_back(neighbourWeight);
    sumOfWeights += neighbourWeight;
  }
}

template<typename T>
void PatchTool2<T>::ComputeWeightedMeanOfPatches(itkTImagePointer & image, typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours, Patch2<float> & outputPatch, double & outputWeight, double & maxWeight, double & smoothing)
{
  btk::Patch2<T> patch = btk::Patch2<T>(image, m_FullPatchSize, point);
  outputPatch.GetData()->FillBuffer(0);

  btk::Patch2<T> neighbourPatch = btk::Patch2<T>(image, m_FullPatchSize);

  itkFloatIterator itp( outputPatch.GetData(), outputPatch.GetData()->GetLargestPossibleRegion() );
  itkTIterator itq( neighbourPatch.GetData(), neighbourPatch.GetData()->GetLargestPossibleRegion() );

  double neighbourWeight = 0;
  double sum = 0;
  maxWeight = 0;

  for(unsigned int n = 0; n < neighbours.size(); n++)
  {
    neighbourPatch.ComputePatch(image, neighbours[n]);
    neighbourWeight = exp( - ComputeL2NormBetweenPatches(patch, neighbourPatch ) / smoothing);
    sum += neighbourWeight;

    if(neighbourWeight > maxWeight)
        maxWeight = neighbourWeight;

    for(itp.GoToBegin(), itq.GoToBegin(); !itp.IsAtEnd(); ++itp, ++itq)
      itp.Set( itp.Get() + neighbourWeight * itq.Get() );
  }

  outputWeight = sum;
}

template<typename T>
void PatchTool2<T>::ComputeWeightedMeanAtPatchCenter(itkTImagePointer & image, typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours, double & outputValue, double & outputWeight, double & maxWeight, double & smoothing)
{
  btk::Patch2<T> patch = btk::Patch2<T>(image, m_FullPatchSize, point);
  outputValue = 0.0;

  btk::Patch2<T> neighbourPatch = btk::Patch2<T>(image, m_FullPatchSize);

  itkTIterator itq( neighbourPatch.GetData(), neighbourPatch.GetData()->GetLargestPossibleRegion() );

  double neighbourWeight = 0;
  double sum = 0;
  maxWeight = 0;

  for(unsigned int n = 0; n < neighbours.size(); n++)
  {
    neighbourPatch.ComputePatch(image, neighbours[n]);
    neighbourWeight = exp( - ComputeL2NormBetweenPatches(patch, neighbourPatch ) / smoothing);
    sum += neighbourWeight;

    if(neighbourWeight > maxWeight)
        maxWeight = neighbourWeight;

    outputValue += (neighbourWeight * neighbourPatch.GetDataAtCenter() );
  }

  outputWeight = sum;
}

template <typename T>
void PatchTool2<T>::ComputePatchRegion(typename itkTImage::IndexType & p, typename itkTImage::RegionType & imageRegion, typename itkTImage::RegionType & patchRegion)
{
  //create an appropriate patch region around the current pixel p, and consider an offset for pixels close to boundaries.
  //the patch region is defined in the image coordinate system.
  //the offset allows to do the link between the patch coordinate system and the one of the image.
  typename itkTImage::RegionType::IndexType start;
  typename itkTImage::RegionType::IndexType offset;
  typename itkTImage::RegionType::SizeType size;

  for(unsigned int i=0; i!= p.GetIndexDimension(); i++)
  {
    start[i] = p[i] - m_HalfPatchSize[i];
    size[i] = m_FullPatchSize[i];
    offset[i] = 0;

    if(start[i] < 0)
    {                  //if the starting index is outside the image (<0)
      size[i] += start[i];
      offset[i] = -start[i];
      start[i] = 0;
    }
    else if(start[i] >= (int)m_ImageSize[i])
    {
      //if the starting index is outside the image (>image size)
      start[i] = m_ImageSize[i]-1;
      size[i] = 0;
    }

    int d = (start[i] + size[i]) - m_ImageSize[i]; //if the region is not fully inside the image
    if(d>0)
    {
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

template <typename T>
void PatchTool2<T>::AddPatchToImage(typename itkTImage::IndexType & p, Patch2<float> & patch, itkFloatImagePointer & image, itkFloatImagePointer & weightImage, double & weight)
{
  //this function add a patch value to the denoised image
  typename itkFloatImage::RegionType imageRegion;
  typename itkFloatImage::RegionType patchRegion;
  ComputePatchRegion(p,imageRegion,patchRegion);

  itkFloatIterator imageIt( image, imageRegion);
  itkFloatIterator weightIt( weightImage, imageRegion);
  itkFloatIterator patchIt( patch.GetData(), patchRegion);

  for ( imageIt.GoToBegin(), patchIt.GoToBegin(), weightIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt, ++patchIt, ++weightIt)
  {
    imageIt.Set( imageIt.Get() + patchIt.Get() );
    weightIt.Set( weightIt.Get() + weight );
  }
}

template <typename T>
void PatchTool2<T>::ComputePatchSimilarityOverMask(itkTImagePointer & image, itkTImagePointer & maskImage, std::vector< typename itkTImage::IndexType> & seeds, itkFloatImagePointer & outputImage, double & smoothing)
{
  typename itkTImage::RegionType  region  = image->GetLargestPossibleRegion();
  typename itkTImage::SizeType    size    = region.GetSize();

  //Loop over seeds
  for(unsigned int s=0; s<seeds.size(); s++)
  {
    btk::Patch2<T> patch = btk::Patch2<T>(image, m_FullPatchSize, seeds[s]);

    int x,y,z;
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)size[2]; z++)
      for(y=0; y < (int)size[1]; y++)
        for(x=0; x < (int)size[0]; x++)
        {
          itkFloatImage::IndexType p;
          p[0] = x;
          p[1] = y;
          p[2] = z;
          if( maskImage->GetPixel(p) > 0 )
          {
            btk::Patch2<T> patch2 = btk::Patch2<T>(image, m_FullPatchSize, p);
            float neighbourWeight = exp( - ComputeL2NormBetweenPatches(patch, patch2 ) / smoothing);
            outputImage->SetPixel(p,neighbourWeight);
          }
        }
  }
}

} // namespace btk

