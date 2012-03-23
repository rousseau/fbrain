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


namespace btk
{
class CreateHRMaskFilter
{
public:
    // Typedef
    typedef float PixelType;
    typedef itk::Image< PixelType, 3>          itkImage;
    typedef itkImage::Pointer                  itkPointer;
    typedef itk::ImageDuplicator< itkImage >   itkDuplicator;
    typedef itk::Image< unsigned char, 3 >     itkImageMask;
    typedef itk::ImageMaskSpatialObject< 3 >   itkMask;
    typedef itkMask::Pointer                   itkMaskPointer;

    typedef itk::AffineTransform<double,3>     itkAffineDeformation;

    typedef btk::SliceBySliceTransform<double,3> SbSTransformType;
    typedef SbSTransformType::Pointer            SbSTransformPointer;


    typedef itk::Transform<double,3> TransformType;

    typedef itk::Euler3DTransform<double> EulerTransformType;
    typedef EulerTransformType::Pointer EulerTransformPointerType;
    typedef std::vector< std::vector< EulerTransformPointerType > > EulerTransformArrayType;

    typedef itk::ResampleImageFilter<itkImage, itkImage>                    itkResampleFilter;
    typedef itk::IdentityTransform<double, 3>                               itkIdentityTransform;
    typedef itk::LinearInterpolateImageFunction<itkImage, double>           itkLinearInterpolator;
    typedef itk::BSplineInterpolateImageFunction<itkImage, double, double>  itkBSplineInterpolator;
    typedef itk::SubtractImageFilter <itkImage, itkImage >                  itkSubtractImageFilter;
    typedef itk::AddImageFilter <itkImage, itkImage >                       itkAddImageFilter;
    typedef itk::BinaryThresholdImageFilter <itkImage, itkImage>            itkBinaryThresholdImageFilter;
    typedef itk::AbsoluteValueDifferenceImageFilter <itkImage, itkImage, itkImage>  itkAbsoluteValueDifferenceImageFilter;
    typedef itk::StatisticsImageFilter<itkImage>                            itkStatisticsImageFilter;

public :
    CreateHRMaskFilter();
    ~CreateHRMaskFilter();
    void Update();

    itkPointer GetOutput()
    {
        return m_maskHRImage;
    }

    void SetInputs(std::vector< itkPointer > inputs)
    {
        m_inputLRImages = inputs;
    }

    std::vector< itkPointer > GetInputs()
    {
        return m_inputLRImages;
    }

    void SetTransforms(std::vector< TransformType::Pointer > transforms)
    {
        m_transforms = transforms;
    }

    std::vector< TransformType::Pointer > GetTransforms()
    {
        return m_transforms;
    }



protected:
private:

    std::vector< TransformType::Pointer > m_transforms;
    std::vector< itkPointer > m_inputLRImages;
    itkPointer m_maskHRImage;
    itkPointer m_HRImage;


};
}
#endif
