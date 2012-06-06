#ifndef __BTK_SUPERRESOLUTIONTYPE_H__
#define __BTK_SUPERRESOLUTIONTYPE_H__

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


#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMatrixOffsetTransformBase.h"


/* BTK */
#include "btkLowToHighImageResolutionMethod.h"
#include "btkSliceBySliceTransform.h"
#include "btkSliceBySliceTransformBase.h"
#include "btkAffineSliceBySliceTransform.h"
#include "btkEulerSliceBySliceTransform.h"

namespace btk
{



typedef float PixelType;
typedef double TScalarType;
const unsigned int Dimension = 3;
typedef itk::Image< PixelType, Dimension>         itkImage;
typedef itk::Image< PixelType, 2>         itkSlice;

typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
typedef itk::ImageRegionIterator< itkImage > itkIterator;
typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;
typedef itk::ContinuousIndex<TScalarType,Dimension>     itkContinuousIndex;


typedef itk::ResampleImageFilter<itkImage, itkImage>                    itkResampleFilter;
typedef itk::ImageToImageFilter< itkImage, itkImage >                   itkFilter;
typedef itk::IdentityTransform<TScalarType, Dimension>                               itkIdentityTransform;
typedef itk::LinearInterpolateImageFunction<itkImage, TScalarType>           itkLinearInterpolator;
typedef itk::BSplineInterpolateImageFunction<itkImage, TScalarType, TScalarType>  itkBSplineInterpolator;
typedef itk::SubtractImageFilter <itkImage, itkImage >                  itkSubtractImageFilter;
typedef itk::AddImageFilter <itkImage, itkImage >                       itkAddImageFilter;
typedef itk::BinaryThresholdImageFilter <itkImage, itkImage>            itkBinaryThresholdImageFilter;
typedef itk::AbsoluteValueDifferenceImageFilter <itkImage, itkImage, itkImage>  itkAbsoluteValueDifferenceImageFilter;
typedef itk::StatisticsImageFilter<itkImage>    itkStatisticsImageFilter;

typedef itk::Image< unsigned char, Dimension >     itkImageMask;
typedef itk::ImageMaskSpatialObject< Dimension >   itkMask;

typedef itkImage::RegionType               itkRegion;

typedef itk::AffineTransform<TScalarType,Dimension>     itkAffineTransform;

typedef itk::Euler3DTransform<TScalarType>      itkEulerTransform;

typedef btk::EulerSliceBySliceTransform<TScalarType,Dimension> btkEulerSliceBySliceTransform;

typedef btk::AffineSliceBySliceTransform< TScalarType, Dimension> btkAffineSliceBySliceTransform;

typedef btk::SliceBySliceTransform<TScalarType,Dimension> btkSliceBySliceTransform;

typedef btk::SliceBySliceTransformBase<TScalarType,Dimension> btkSliceBySliceTransformBase;


typedef itk::Euler3DTransform<TScalarType> itkEulerTransform;


typedef itk::MatrixOffsetTransformBase<TScalarType,Dimension,Dimension> itkMatrixTransform;

typedef itk::Transform<TScalarType, Dimension> itkTransformBase;

typedef btk::LowToHighImageResolutionMethod<itkImage, itkEulerTransform> LowToHighResFilterRigid;
typedef btk::LowToHighImageResolutionMethod<itkImage, itkAffineTransform> LowToHighResFilterAffine;
enum TRANSFORMATION_TYPE
{
    AFFINE = 0,
    EULER_3D,
    SLICE_BY_SLICE_AFFINE,
    SLICE_BY_SLICE_EULER,
    SLICE_BY_SLICE,
    ANTS
};

// Example of reconstruction type :
enum RECONSTRUCTION_TYPE
{
    SR = 0,
    IBP
};
enum PSF_TYPE
{
    BOXCAR = 0,
    GAUSSIAN
};
}




#endif // BTKSUPERRESOLUTIONTYPE_H
