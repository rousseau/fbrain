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
#include "btkSliceBySliceTransform.h"


namespace btk
{



typedef float PixelType;
typedef itk::Image< PixelType, 3>         itkImage;
typedef itk::Image< PixelType, 2>         itkSlice;

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

typedef itk::AffineTransform<double,3>     itkAffineTransform;

typedef btk::SliceBySliceTransform<double,3> btkSliceBySliceTransform;


typedef itk::Euler3DTransform<double> itkEulerTransform;


typedef itk::MatrixOffsetTransformBase<double,3,3> itkMatrixTransform;

typedef itk::Transform<double, 3> itkTransformBase;

enum TRANSFORMATION_TYPE
{
    AFFINE = 0,
    EULER_3D,
    SLICE_BY_SLICE,
    ANTS
};

// Example of reconstruction type :
enum RECONSTRUCTION_TYPE
{
    IBP = 0,
    SR
};
}




#endif // BTKSUPERRESOLUTIONTYPE_H
