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

#ifndef BTK_PATCHTOOL_H
#define BTK_PATCHTOOL_H

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkStatisticsImageFilter.h"
#include "itkChiSquareDistribution.h"


#include "btkPatch.h"

namespace btk
{

  /**
  * @class PatchTool
  * @brief 3D patch handling
  */
  template<typename T1, typename T2>
  class PatchTool
  {
    public:

    //Defining ITK stuff for the input image and the patch image
    typedef typename itk::Image< float, 3>                         itkFloatImage;
    typedef typename itkFloatImage::Pointer                        itkFloatImagePointer;
    typedef typename itk::ImageRegionIterator< itkFloatImage >     itkFloatIterator;

    typedef typename itk::Image< T1, 3>                           itkT1Image;
    typedef typename itkT1Image::Pointer                          itkT1ImagePointer;
    typedef typename itk::ImageRegionIterator< itkT1Image >       itkT1Iterator;
    typedef typename itk::ImageRegionConstIterator< itkT1Image >  itkConstT1Iterator;
  	typedef typename itk::ImageRegionIteratorWithIndex< itkT1Image > itkT1IteratorWithIndex;

    typedef typename itk::Image< T2, 3>                           itkT2Image;
    typedef typename itkT2Image::Pointer                          itkT2ImagePointer;
    typedef typename itk::ImageRegionIterator< itkT2Image >       itkT2Iterator;
    typedef typename itk::ImageRegionConstIterator< itkT2Image >  itkConstT2Iterator;
  	typedef typename itk::ImageRegionIteratorWithIndex< itkT2Image > itkT2IteratorWithIndex;

    typedef typename itk::StatisticsImageFilter<itkT1Image> itkStatisticsT1ImageFilter;
    
    void SetSpatialBandwidth(int s, Patch<T1> & patch);

    void AddPatchToImage(typename itkT1Image::IndexType & p, Patch<T2> & patch, itkT1ImagePointer & image, itkFloatImagePointer & weightImage, double weight);

    void PatchIntensityNormalizationUsingMeanAndVariance(Patch<T1> & inputPatch, Patch<T1> & refPatch, Patch<T2> & outputPatch);

	void CreatePatchFromAnotherPatch(Patch<T1> & inputPatch, Patch<T2> & outputPatch);
		
	void ComputeSearchRegion(Patch<T1> & inputPatch,  typename itkT1Image::RegionType & region);	
		
	void SetSelectionNeighbourPatchMethod(Patch<T1> & inputPatch);
		
	void GetNeighbourPatches(Patch<T1> & inputPatch, std::vector< Patch<T2> > & neighbourPatches, itkT2ImagePointer & image);

	void ComputeNeighbourWeights(Patch<T1> & inputPatch, std::vector< Patch<T2> > & neighbourPatches, std::vector<float> & weights, float & smoothing);
	
	double ComputeL2NormBetweenPatches(Patch<T1> & p, Patch<T2> & q);
	
   	btkGetMacro(HalfSpatialBandwidth,typename itkT1Image::SizeType);
   	btkGetMacro(FullSpatialBandwidth,typename itkT1Image::SizeType);
   	btkGetMacro(PatchSelectionMethod, int);
	btkSetMacro(PatchSelectionMethod, int);
	btkSetMacro(ParamPatchSelectionMethod, float);
	btkSetMacro(ImageSize, typename itkT1Image::SizeType);

  private:
	typename itkT1Image::SizeType    m_HalfSpatialBandwidth;   //equivalent to the half size of the volume search area in non-local means
  	typename itkT1Image::SizeType    m_FullSpatialBandwidth;   //spatial bandwidth : 2 * halfSpatialBandwidth + 1
  	int                              m_PatchSelectionMethod;
  	float                            m_ParamPatchSelectionMethod;
  	typename itkT1Image::SizeType    m_ImageSize;
  	
  };

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkPatchTool.txx"
#endif

#endif // BTK_PATCHTOOL_H
