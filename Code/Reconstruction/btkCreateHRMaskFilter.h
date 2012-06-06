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

#ifndef __BTK_CREATEHRMASKFILTER_H__
#define __BTK_CREATEHRMASKFILTER_H__

#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkTransformFactory.h"
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

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMatrixOffsetTransformBase.h"

#include "btkSliceBySliceTransform.h"
#include "btkMacro.h"
#include "btkSuperResolutionType.h"


namespace btk
{
class CreateHRMaskFilter
{
public:

    typedef itk::Euler3DTransform<double> EulerTransformType;
    typedef EulerTransformType::Pointer EulerTransformPointerType;
    typedef std::vector< std::vector< EulerTransformPointerType > > EulerTransformArrayType;



public :
    CreateHRMaskFilter();
    ~CreateHRMaskFilter();
    void Update();

    itkImage::Pointer GetOutput()
    {
        return m_MaskHRImage;
    }


    btkGetMacro(InputLRImages,std::vector< itkImage::Pointer >);
    btkSetMacro(InputLRImages,std::vector< itkImage::Pointer >);

    btkGetMacro(Transforms, std::vector< itkTransformBase::Pointer >);
    btkSetMacro(Transforms, std::vector< itkTransformBase::Pointer >);

//    btkGetMacro(Transforms, std::vector< itkAffineTransform::Pointer >);
//    btkSetMacro(Transforms, std::vector< itkAffineTransform::Pointer >);

    btkGetMacro(TransformsAffine, std::vector< itkAffineTransform::Pointer >);
    btkSetMacro(TransformsAffine, std::vector< itkAffineTransform::Pointer >);

    btkGetMacro(TransformsSbS, std::vector< btkEulerSliceBySliceTransform::Pointer >);
    btkSetMacro(TransformsSbS, std::vector< btkEulerSliceBySliceTransform::Pointer >);

    btkGetMacro(HRImage,itkImage::Pointer);
    btkSetMacro(HRImage, itkImage::Pointer);

    btkGetMacro(TransformType,TRANSFORMATION_TYPE);
    btkSetMacro(TransformType,TRANSFORMATION_TYPE);







protected:
private:

    std::vector< itkTransformBase::Pointer > m_Transforms;
//    std::vector< itkAffineTransform::Pointer > m_Transforms;

    std::vector< itkAffineTransform::Pointer > m_TransformsAffine;
    std::vector< btkEulerSliceBySliceTransform::Pointer > m_TransformsSbS;

    std::vector< itkImage::Pointer > m_InputLRImages;
    itkImage::Pointer m_MaskHRImage;
    itkImage::Pointer m_HRImage;

    TRANSFORMATION_TYPE m_TransformType;



};
}
#endif
