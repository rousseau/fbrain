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
#include "itkTransformFileReader.h"

/* VNL */
#include "vnl/vnl_sparse_matrix.h"

/* BTK */
#include "btkSliceBySliceTransform.h"
#include "btkCreateHRMaskFilter.h"
#include "btkSimulateLRImageFilter.h"
#include "btkMacro.h"
//Typedefs :
#include "btkSuperResolutionType.h"


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

    btkGetMacro(TransformsLR,std::vector< itkTransformBase::Pointer >);// To remove
    btkSetMacro(TransformsLR,std::vector< itkTransformBase::Pointer >);


    btkGetMacro(TransformsLRSbS,std::vector< btkSliceBySliceTransform::Pointer >);
    btkSetMacro(TransformsLRSbS,std::vector< btkSliceBySliceTransform::Pointer >);

    btkGetMacro(TransformsLRAffine,std::vector< itkAffineTransform::Pointer >);
    btkSetMacro(TransformsLRAffine,std::vector< itkAffineTransform::Pointer >);


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

    btkGetMacro(MedianIBP, int);
    btkSetMacro(MedianIBP,int);

    btkGetMacro(Psftype, int);
    btkSetMacro(Psftype,int);

    btkGetMacro(InterpolationOrderPSF, int);
    btkSetMacro(InterpolationOrderPSF,int);

    btkGetMacro(InterpolationOrderIBP, int);
    btkSetMacro(InterpolationOrderIBP,int);


    btkGetMacro(TransformationType, TRANSFORMATION_TYPE);
    btkSetMacro(TransformationType, TRANSFORMATION_TYPE);

    btkGetMacro(ReconstructionType,RECONSTRUCTION_TYPE);
    btkSetMacro(ReconstructionType,RECONSTRUCTION_TYPE);



    void SetParameters(int Nlm, float Beta, int Loop, int MedianIBP, int PsfType, int InterpolationOrderIBP, int InterpolationOrderPSF);
    void SetDefaultParameters();

protected:

    void InverseTransforms();
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
    std::vector< itkTransformBase::Pointer  >   m_TransformsLR;
    std::vector< itkTransformBase::Pointer >    m_InverseTransformsLR;

    std::vector< btkSliceBySliceTransform::Pointer  >   m_TransformsLRSbS;
    std::vector< btkSliceBySliceTransform::Pointer >    m_InverseTransformsLRSbS;

    std::vector< itkAffineTransform::Pointer  >   m_TransformsLRAffine;
    std::vector< itkAffineTransform::Pointer >    m_InverseTransformsLRAffine;

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
    int m_MedianIBP;

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
