/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 january 2013
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

#ifndef BTK_PATCHTOOL2_H
#define BTK_PATCHTOOL2_H

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkStatisticsImageFilter.h"
#include "itkChiSquareDistribution.h"


#include "btkPatch2.h"

namespace btk
{

  /**
  * @class PatchTool
  * @brief 3D patch handling
  */
  template<typename T>
  class PatchTool2
  {
    public:

    //Defining ITK stuff for the input image and the patch image
    typedef typename itk::Image< T, 3>                              itkTImage;
    typedef typename itkTImage::Pointer                             itkTImagePointer;
    typedef typename itk::ImageRegionIterator< itkTImage >          itkTIterator;
    typedef typename itk::ImageRegionConstIterator< itkTImage >     itkConstTIterator;
    typedef typename itk::ImageRegionIteratorWithIndex< itkTImage > itkTIteratorWithIndex;

    typedef typename itk::Image< float, 3>                              itkFloatImage;
    typedef typename itkFloatImage::Pointer                             itkFloatImagePointer;
    typedef typename itk::ImageRegionIterator< itkFloatImage >          itkFloatIterator;

    typedef typename itk::StatisticsImageFilter<itkTImage> itkStatisticsTImageFilter;
    


    void ComputePatchSize(itkTImagePointer & image, int & h);
    void ComputeSpatialBandwidth(itkTImagePointer & image, int & h);
    void ComputeSearchRegion(typename itkTImage::IndexType & point, itkTImagePointer & image, typename itkTImage::RegionType & region);
    double ComputeL2NormBetweenPatches(Patch2<T> & p, Patch2<T> & q);

    void ComputeMeanAndVarianceImage(itkTImagePointer & image, itkTImagePointer & maskImage, itkFloatImagePointer & meanImage, itkFloatImagePointer & varianceImage);

    //Using vector of patches is really too slow (because of the operator New:: in ITK image)
    //So we have to store only the coordinates of the patch center into vectors...
    void GetNeighboursUsingMeanAndVariance(typename itkTImage::IndexType & point, itkTImagePointer & image, float & mean, float & variance, itkFloatImagePointer & meanImage, itkFloatImagePointer & varianceImage, std::vector< typename itkTImage::IndexType > & neighbours);
    void GetNeighbours(typename itkTImage::IndexType & point, itkTImagePointer & image, std::vector< typename itkTImage::IndexType > & neighbours);
    void RemoveCentralPointInNeighborhood(typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours);

    void ComputeNeighbourWeights(itkTImagePointer & image, typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours, std::vector<double> & weights, double & sumOfWeights, double & smoothing);
    void ComputeWeightedMeanOfPatches(itkTImagePointer & image, typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours, Patch2<float> & outputPatch, double & outputWeight, double & maxWeight, double & smoothing);
    void ComputeWeightedMeanAtPatchCenter(itkTImagePointer & image, typename itkTImage::IndexType & point, std::vector< typename itkTImage::IndexType > & neighbours, double & outputValue, double & outputWeight, double & maxWeight, double & smoothing);

    void ComputePatchRegion(typename itkTImage::IndexType & p, typename itkTImage::RegionType & imageRegion, typename itkTImage::RegionType & patchRegion);
    void AddPatchToImage(typename itkTImage::IndexType & p, Patch2<float> & patch, itkFloatImagePointer & image, itkFloatImagePointer & weightImage, double & weight);

    void ComputePatchSimilarityOverMask(itkTImagePointer & image, itkTImagePointer & maskImage, std::vector<typename itkTImage::IndexType> & seeds, itkFloatImagePointer & outputImage, double & smoothing);

    //TO BE IMPLEMENTED : CHI2, Gestion du patch central, CorrelationBetweenNeighbourhood
    //-----------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------
    /*
    void SetSpatialBandwidth(int s, Patch<T1> & patch);

    void AddPatchToImage(typename itkT1Image::IndexType & p, Patch<T2> & patch, itkT1ImagePointer & image, itkFloatImagePointer & weightImage, double weight);

    void PatchIntensityNormalizationUsingMeanAndVariance(Patch<T1> & inputPatch, Patch<T1> & refPatch, Patch<T2> & outputPatch);

	void CreatePatchFromAnotherPatch(Patch<T1> & inputPatch, Patch<T2> & outputPatch);
		
	void ComputeSearchRegion(Patch<T1> & inputPatch,  typename itkT1Image::RegionType & region);	
		
	void SetSelectionNeighbourPatchMethod(Patch<T1> & inputPatch);
		
	void GetNeighbourPatches(Patch<T1> & inputPatch, std::vector< Patch<T2> > & neighbourPatches, itkT2ImagePointer & image);

	void ComputeNeighbourWeights(Patch<T1> & inputPatch, std::vector< Patch<T2> > & neighbourPatches, std::vector<float> & weights, float & smoothing);
	
	double ComputeL2NormBetweenPatches(Patch<T1> & p, Patch<T2> & q);
    */

    btkSetMacro(HalfPatchSize,typename itkTImage::SizeType);
    btkGetMacro(HalfPatchSize,typename itkTImage::SizeType);
    btkSetMacro(FullPatchSize,typename itkTImage::SizeType);
    btkGetMacro(FullPatchSize,typename itkTImage::SizeType);
    btkGetMacro(HalfSpatialBandwidth,typename itkTImage::SizeType);
    btkGetMacro(FullSpatialBandwidth,typename itkTImage::SizeType);
   	btkGetMacro(PatchSelectionMethod, int);
	btkSetMacro(PatchSelectionMethod, int);
	btkSetMacro(ParamPatchSelectionMethod, float);
    btkSetMacro(ImageSize, typename itkTImage::SizeType);

  private:
    typename itkTImage::SizeType    m_HalfPatchSize;          //half of the patch size
    typename itkTImage::SizeType    m_FullPatchSize;          //patch size  : 2 * halfPatchSize + 1
    typename itkTImage::SizeType    m_HalfSpatialBandwidth;   //equivalent to the half size of the volume search area in non-local means
    typename itkTImage::SizeType    m_FullSpatialBandwidth;   //spatial bandwidth : 2 * halfSpatialBandwidth + 1
    int                             m_PatchSelectionMethod;
    float                           m_ParamPatchSelectionMethod;
    typename itkTImage::SizeType    m_ImageSize;
  	
  };

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkPatchTool2.txx"
#endif

#endif // BTK_PATCHTOOL2_H
