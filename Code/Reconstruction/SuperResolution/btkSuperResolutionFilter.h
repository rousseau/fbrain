/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 16/03/2012
  Author(s): Schweitzer Marc (marc.schweitzer@unistra.fr)

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

#ifndef __BTK_SUPERRESOLUTIONFILTER_h__
#define __BTK_SUPERRESOLUTIONFILTER_h__

/* ITK */
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"
#include "itkBSplineInterpolationWeightFunction.h"
#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkTransformFactory.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMatrixOffsetTransformBase.h"

/* VNL */
#include "vnl/vnl_sparse_matrix.h"

/* BTK */
#include "btkSliceBySliceTransform.h"
#include "btkCreateHRMaskFilter.h"
#include "btkSimulateLRImageFilter.h"
#include "btkMacro.h"


// Filter used for SuperResolution pipeline :
#include "btkBiasCorrectionFilter.h"
#include "btkHighResolutionReconstructionFilter.h"
#include "btkMotionCorrectionFilter.h"
#include "btkPSFEstimationFilter.h"
#include "btkSliceRejectionFilter.h"




namespace btk
{
class SuperResolutionFilter
{
public:
    typedef float PixelType;
    typedef itk::Image< PixelType, 3>         itkImage;
    typedef itk::Image< PixelType, 2>         SliceType;

    typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
    typedef itk::ImageRegionIterator< itkImage > itkIterator;
    typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;
    typedef itk::ContinuousIndex<double,3>     itkContinuousIndex;


    typedef itk::ResampleImageFilter<itkImage, itkImage>                    itkResampleFilter;
    typedef itk::IdentityTransform<double, 3>                               itkIdentityTransform;
    typedef itk::LinearInterpolateImageFunction<itkImage, double>           itkLinearInterpolator;
    typedef itk::BSplineInterpolateImageFunction<itkImage, double, double>  itkBSplineInterpolator;
    typedef itk::SubtractImageFilter <itkImage, itkImage >                  itkSubtractImageFilter;
    typedef itk::AddImageFilter <itkImage, itkImage >                       itkAddImageFilter;
    typedef itk::BinaryThresholdImageFilter <itkImage, itkImage>            itkBinaryThresholdImageFilter;
    typedef itk::AbsoluteValueDifferenceImageFilter <itkImage, itkImage, itkImage>  itkAbsoluteValueDifferenceImageFilter;
    typedef itk::StatisticsImageFilter<itkImage>    itkStatisticsImageFilter;

    typedef itk::Image< unsigned char, 3 >     itkImageMask;
    typedef itk::ImageMaskSpatialObject< 3 >   itkMask;

    typedef itk::AffineTransform<double,3>     itkAffineDeformation;

    typedef btk::SliceBySliceTransform<double,3> SbSTransformType;
    typedef SbSTransformType::Pointer            SbSTransformPointer;

    typedef itk::Euler3DTransform<double> EulerTransformType;
    typedef EulerTransformType::Pointer EulerTransformPointerType;
    typedef std::vector< std::vector< EulerTransformPointerType > > EulerTransformArrayType;

    typedef itk::MatrixOffsetTransformBase<double,3,3> MatrixTransformType;

    typedef itk::Transform<double, 3> TransformType;

    enum TRANSFORMATION_TYPE
    {
        AFFINE = 0,
        EULER_3D,
        SLICE_BY_SLICE
    };

    // Example of reconstruction type :
    enum RECONSTRUCTION_TYPE
    {
        ALGO1 = 0,
        ALGO2,
        ALGO3
    };


    /** Constructor with no parameters */
    SuperResolutionFilter();

    /** Constructors with parameters */
    SuperResolutionFilter(int loop, float beta, int nlm);
    SuperResolutionFilter(int loop, float beta, int nlm, TRANSFORMATION_TYPE transfoType);
    // TODO : Add the others constructors with wanted parameters
    //...


    /** Destructor */
    ~SuperResolutionFilter();

    /** Method called for compute the SuperResolution pipeline */
    void Update();

    /** Method who return the High Resolution image  */
    itkImage::Pointer GetOutput();



    // GETTER/SETTER :

    btkGetMacro(TransformsLR,std::vector< TransformType::Pointer >);
    btkSetMacro(TransformsLR,std::vector< TransformType::Pointer >);

    btkGetMacro(ImagesLR,std::vector< itkImage::Pointer > );
    btkSetMacro(ImagesLR,std::vector< itkImage::Pointer > );


    btkGetMacro(ImagesMaskLR,std::vector< itkImageMask::Pointer >  );
    btkSetMacro(ImagesMaskLR,std::vector< itkImageMask::Pointer >  );

    btkGetMacro(MasksLR,std::vector< itkMask::Pointer >);
    btkSetMacro(MasksLR,std::vector< itkMask::Pointer >);

    btkGetMacro(ImageHR, itkImage::Pointer);
    btkSetMacro(ImageHR, itkImage::Pointer);

    btkGetMacro(ImageMaskHR, itkImageMask::Pointer);
    btkSetMacro(ImageMaskHR, itkImageMask::Pointer);

    btkGetMacro(ReferenceImage, itkImage::Pointer);
    btkSetMacro(ReferenceImage, itkImage::Pointer);

    btkGetMacro(Nlm, int);
    btkSetMacro(Nlm,int);

    btkGetMacro(Beta, float);
    btkSetMacro(Beta,float);

    btkGetMacro(Loop, int);
    btkSetMacro(Loop,int);







protected:

private:

    // !! Since use of btkMacro for get/set use a Capital letter after the m_ (ex : m_MyVariable) !!

    int                       m_Psftype; // 0: 3D interpolated boxcar, 1: 3D oversampled boxcar
    std::vector<itkImage::Pointer>   m_PSF;
    int                       m_InterpolationOrderPSF;
    int                       m_InterpolationOrderIBP;
    vnl_sparse_matrix<float>  m_H;
    vnl_vector<float>         m_Y;
    vnl_vector<float>         m_X;
    float                     m_PaddingValue;
    std::vector<unsigned int> m_Offset;
    TRANSFORMATION_TYPE m_TransformationType;
    RECONSTRUCTION_TYPE m_ReconstructionType;


    std::vector< itkImage::Pointer >         m_ImagesLR;
    std::vector< TransformType::Pointer  >   m_TransformsLR;
    std::vector< TransformType::Pointer >    m_InverseTransformsLR;
    std::vector< itkImageMask::Pointer >     m_ImagesMaskLR;
    std::vector< itkMask::Pointer >          m_MasksLR;
    itkImage::Pointer                        m_ImageHR;
    itkImage::Pointer                        m_OutputHRImage;
    itkImageMask::Pointer                    m_ImageMaskHR;
    itkMask::Pointer                         m_MaskHR;
    itkImage::Pointer                        m_ReferenceImage;

    int m_Nlm;
    float m_Beta;
    int m_Loop;

    //Image information (size, spacing etc.)
    itkImage::SpacingType m_Spacing;
    itkImage::SizeType    m_Size;
    itkImage::RegionType  m_Region;

    btk::PSFEstimationFilter * m_PSFEstimationFilter;
    btk::MotionCorrectionFilter * m_MotionCorrectionFilter;
    btk::BiasCorrectionFilter *   m_BiasCorrectionFilter;
    btk::HighResolutionReconstructionFilter * m_HighResolutionReconstructionFilter;
    btk::SliceRejectionFilter * m_SliceRejectionFilter;




};//end class

}//end namespace

#endif
