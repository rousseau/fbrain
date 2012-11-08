/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 24/10/2012
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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

#ifndef BTKLOWTOHIGHRESOLUTIONFILTER_HXX
#define BTKLOWTOHIGHRESOLUTIONFILTER_HXX

/* ITK */
#include "itkObject.h"
#include "itkExceptionObject.h"
#include "itkImage.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"

/* STD */
#include "vector"

/* BTK */
#include "btkMacro.h"

namespace btk
{
/** \class LowToHighResolutionFilter
 * \brief Class for obtaining an isovoxel image from low resolution images, and previously
 * computed transformations.
 *
 *\ingroup Reconstruction
 */
template <typename TImage>
class LowToHighResolutionFilter : public itk::Object
{
    public:
    /** typedefs : */
        typedef LowToHighResolutionFilter Self;
        typedef itk::Object SuperClass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        typedef TImage ImageType;
        typedef itk::MatrixOffsetTransformBase< double , 3> TransformType;
        typedef itk::LinearInterpolateImageFunction< ImageType, double> Interpolator;
        typedef itk::Image<unsigned char, 3>    MaskType;
        typedef itk::LinearInterpolateImageFunction< MaskType, double> MaskInterpolator;

        typedef itk::ImageRegionIteratorWithIndex< ImageType> Iterator;
        typedef itk::ImageRegionIteratorWithIndex<MaskType> MaskIterator;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(LowToHighResolutionFilter, itk::Object);

        /** Update Method */
        virtual void Update()throw(itk::ExceptionObject);

        /** Resample each low resolution image on the HighResolution space */

        std::vector<typename ImageType::Pointer>& ResampleImages() throw(itk::ExceptionObject);


        /** GET/SET */

        btkSetMacro(Images, std::vector<typename ImageType::Pointer>);

        btkSetMacro(Masks,std::vector<typename MaskType::Pointer>);

        btkSetMacro(Transforms,std::vector<typename TransformType::Pointer>);
        btkGetMacro(InverseTransforms,std::vector<typename TransformType::Pointer>);

        btkSetMacro(ReferenceImage,typename ImageType::Pointer);
        btkSetMacro(ReferenceMask, typename MaskType::Pointer);

        btkSetMacro(UseReference,bool);
        btkGetMacro(UseReference,bool);

        btkSetMacro(ReferenceImageNum, unsigned int);
        btkGetMacro(ReferenceImageNum, unsigned int);



        btkGetMacro(Output, typename ImageType::Pointer);

        btkGetMacro(MaskCombination, typename MaskType::Pointer );




    protected:
        /** Constructor */
        LowToHighResolutionFilter();
        /** Destructor */
        virtual ~LowToHighResolutionFilter(){};

        /** Initialize method */
        virtual void Initialize() throw(itk::ExceptionObject);

    private:

        /** private Members  */
        std::vector<typename ImageType::Pointer> m_Images; /** Inputs Images */
        std::vector<typename MaskType::Pointer> m_Masks; /** Inputs masks */
        std::vector<typename TransformType::Pointer> m_Transforms; /** Inputs Transform */
        std::vector<typename TransformType::Pointer> m_InverseTransforms; /** Inputs Transform */
        typename ImageType::Pointer m_ReferenceImage; /** Reference Image */
        typename MaskType::Pointer m_ReferenceMask; /** Reference Mask */
        unsigned int m_ReferenceImageNum;/** Reference image num, in array of inputs images */
        bool m_UseReference;

        typename ImageType::Pointer m_Output; /** Output Image */
        typename ImageType::Pointer m_HighResolutionImage; /** HighResolution Image */
        typename MaskType::Pointer m_MaskCombination; /** Combination of input masks */
        std::vector<typename ImageType::Pointer> m_ResampledImages; /** Resampled Low resolution images */

        std::vector<typename Interpolator::Pointer> m_Interpolators; /** Linear interpolators on images */


        typename ImageType::RegionType m_ReferenceRegion;


        unsigned int m_NumberOfImages; /** Number of inputs images */
        double m_Margin; /** Margin for Highresolution , default = 0 */





};
}
#ifndef ITK_MANUAL_INSTANTIATION
#include "btkLowToHighResolutionFilter.txx"
#endif

#endif // BTKLOWTOHIGHRESOLUTIONFILTER_HXX
