/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 15/10/2013
  Author(s): François Rousseau
  
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

#include "btkPandoraBoxImageFilters.h"

namespace btk
{

void PandoraBoxImageFilters::DiscreteGaussianFiltering(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, float variance)
{
  itkGaussianFilter::Pointer filter = itkGaussianFilter::New();
  filter->SetInput(inputImage);

  itkGaussianFilter::ArrayType var;
  var[0] = variance;
  var[1] = variance;
  var[2] = variance;
  filter->SetVariance(var);

  filter->SetUseImageSpacingOff();
  filter->Update();

  outputImage = filter->GetOutput();
}

void PandoraBoxImageFilters::ResampleImageUsingSpacing(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, itkFloatImage::SpacingType & itkSpacing, int interpolationOrder)
{
  //Warning : interpolation order 0 (nearest neighbors) or 1 (linear) are fine. 3 (spline) may provide null image (ITK bug?).

  itkResampleFilter::Pointer resample = itkResampleFilter::New();

  //parameters for interpolation (identity transform and bspline interpolator)
  itkIdentityTransform::Pointer transform = itkIdentityTransform::New();
  transform->SetIdentity();
  itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
  bsInterpolator->SetSplineOrder(interpolationOrder);

  resample->SetTransform(transform);
  resample->SetInterpolator(bsInterpolator);

  itkFloatImage::SizeType    itkSize;
  itkFloatImage::SizeType    itkInputSize    = inputImage->GetLargestPossibleRegion().GetSize();
  itkFloatImage::SpacingType itkInputSpacing = inputImage->GetSpacing();

  for(unsigned int i=0; i<3; i++)
    itkSize[i] = itkInputSize[i] * itkInputSpacing[i] / itkSpacing[i] ;

  resample->SetSize(itkSize);
  resample->SetOutputSpacing(itkSpacing);
  resample->SetDefaultPixelValue(0.0);
  resample->SetInput( inputImage );

  //Set the new origin correctly:
  itkContinuousIndex outputIndex;
  outputIndex[0] = (-0.5 * inputImage->GetSpacing()[0] + 0.5 * resample->GetOutputSpacing()[0]) / inputImage->GetSpacing()[0] ;
  outputIndex[1] = (-0.5 * inputImage->GetSpacing()[1] + 0.5 * resample->GetOutputSpacing()[1]) / inputImage->GetSpacing()[1] ;
  outputIndex[2] = (-0.5 * inputImage->GetSpacing()[2] + 0.5 * resample->GetOutputSpacing()[2]) / inputImage->GetSpacing()[2] ;

  itkFloatImage::PointType outputPoint;
  inputImage->TransformContinuousIndexToPhysicalPoint(outputIndex,outputPoint);

  resample->SetOutputOrigin( outputPoint );
  resample->SetOutputDirection( inputImage->GetDirection() );
  resample->Update();
  outputImage = resample->GetOutput();
}

void PandoraBoxImageFilters::ResampleImageUsingSize(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, itkFloatImage::SizeType & itkSize, int interpolationOrder)
{
  itkResampleFilter::Pointer resample = itkResampleFilter::New();

  //parameters for interpolation (identity transform and bspline interpolator)
  itkIdentityTransform::Pointer transform = itkIdentityTransform::New();
  transform->SetIdentity();
  itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
  bsInterpolator->SetSplineOrder(interpolationOrder);

  resample->SetTransform(transform);
  resample->SetInterpolator(bsInterpolator);

  itkFloatImage::SpacingType itkSpacing;
  itkFloatImage::SizeType    itkInputSize    = inputImage->GetLargestPossibleRegion().GetSize();
  itkFloatImage::SpacingType itkInputSpacing = inputImage->GetSpacing();

  for(uint i=0; i<3; i++)
    itkSpacing[i] = itkInputSize[i] * itkInputSpacing[i] / itkSize[i] ;

  resample->SetSize(itkSize);
  resample->SetOutputSpacing(itkSpacing);
  resample->SetDefaultPixelValue(0.0);
  resample->SetInput( inputImage );

  //Set the new origin correctly:
  itkContinuousIndex outputIndex;
  outputIndex[0] = (-0.5 * inputImage->GetSpacing()[0] + 0.5 * resample->GetOutputSpacing()[0]) / inputImage->GetSpacing()[0] ;
  outputIndex[1] = (-0.5 * inputImage->GetSpacing()[1] + 0.5 * resample->GetOutputSpacing()[1]) / inputImage->GetSpacing()[1] ;
  outputIndex[2] = (-0.5 * inputImage->GetSpacing()[2] + 0.5 * resample->GetOutputSpacing()[2]) / inputImage->GetSpacing()[2] ;

  itkFloatImage::PointType outputPoint;
  inputImage->TransformContinuousIndexToPhysicalPoint(outputIndex,outputPoint);

  resample->SetOutputOrigin( outputPoint );
  resample->SetOutputDirection( inputImage->GetDirection() );
  resample->Update();
  outputImage = resample->GetOutput();
}

void PandoraBoxImageFilters::ResampleImageUsingReference(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, itkFloatImagePointer & referenceImage, int interpolationOrder)
{
  itkResampleFilter::Pointer resample = itkResampleFilter::New();

  //parameters for interpolation (identity transform and bspline interpolator)
  itkIdentityTransform::Pointer transform = itkIdentityTransform::New();
  transform->SetIdentity();
  itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
  bsInterpolator->SetSplineOrder(interpolationOrder);

  resample->SetTransform(transform);
  resample->SetInterpolator(bsInterpolator);
  resample->UseReferenceImageOn();
  resample->SetReferenceImage( referenceImage );
  resample->SetDefaultPixelValue(0.0);
  resample->SetInput( inputImage );
  resample->Update();
  outputImage = resample->GetOutput();
}

void PandoraBoxImageFilters::ResampleImageUsingTransform(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, itkFloatImagePointer & referenceImage, int interpolationOrder, itkTransformPointer & inputTransform)
  {
    itkResampleFilter::Pointer resample = itkResampleFilter::New();
    
    //parameters for interpolation (bspline interpolator)
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(interpolationOrder);
    
    resample->SetTransform(inputTransform);
    resample->SetInterpolator(bsInterpolator);
    resample->UseReferenceImageOn();
    resample->SetReferenceImage( referenceImage );
    resample->SetDefaultPixelValue(0.0);
    resample->SetInput( inputImage );
    resample->Update();
    outputImage = resample->GetOutput();
  }
  
void PandoraBoxImageFilters::DownscaleImage(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, float factor)
{
  itkFloatImage::Pointer tmpImage;
  float variance = 6 * factor / 3.14159;   //PhD thesis P. Hellier 2000, page 35
  DiscreteGaussianFiltering(inputImage,tmpImage,variance);

  itkFloatImage::SpacingType spacing = inputImage->GetSpacing();
  float minSpacing = spacing[0];
  for(unsigned int i=1; i<3; i++)
    if(minSpacing > spacing[i] )
      minSpacing = spacing[i];

  for(unsigned int i=0; i<3; i++)
    //spacing[i] = spacing[i] * factor;
    spacing[i] = minSpacing * factor;

  ResampleImageUsingSpacing(tmpImage,outputImage,spacing,0);
}


void PandoraBoxImageFilters::ProbabilityImageNormalization(std::vector< itkFloatImagePointer > & inputImages, std::vector< itkFloatImagePointer > & outputImages)
{
  int x,y,z;
  itkFloatImage::SizeType    size    = inputImages[0]->GetLargestPossibleRegion().GetSize();

#pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)size[2]; z++)
  for(y=0; y < (int)size[1]; y++)
  for(x=0; x < (int)size[0]; x++)
  {
    itkFloatImage::IndexType index;
    index[0] = x;
    index[1] = y;
    index[2] = z;

    std::vector<double> valueVector(inputImages.size());
    double sumVector = 0.0;

    for(unsigned int l = 0; l < inputImages.size(); l++){
      valueVector[l] = inputImages[l]->GetPixel( index );
      if(valueVector[l] < 0)
        valueVector[l] = 0;
      sumVector += valueVector[l];
    }

    if(sumVector > 0)
      for(unsigned int l = 0; l < inputImages.size(); l++){
        valueVector[l] /= sumVector;
        outputImages[l]->SetPixel(index, valueVector[l]);
      }
    else
      for(unsigned int l = 0; l < inputImages.size(); l++)
          outputImages[l]->SetPixel(index, 0);
  }

}

void PandoraBoxImageFilters::GetLabelWithMaxProbabilityImage(std::vector< itkFloatImagePointer > & inputImages, itkShortImagePointer & outputImage)
{
  int x,y,z;
  itkFloatImage::SizeType    size    = inputImages[0]->GetLargestPossibleRegion().GetSize();

  #pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)size[2]; z++)
  for(y=0; y < (int)size[1]; y++)
  for(x=0; x < (int)size[0]; x++)
  {
    itkShortImage::IndexType index;
    index[0] = x;
    index[1] = y;
    index[2] = z;

    float maxProba = 0;
    short maxLabel = 0;
    for(unsigned int l = 0; l < inputImages.size(); l++)
      if( maxProba < inputImages[l]->GetPixel( index ) )
      {
        maxProba = inputImages[l]->GetPixel( index );
        maxLabel = l+1;
      }

    outputImage->SetPixel(index, maxLabel);
  }
}

//Display common image info
void PandoraBoxImageFilters::DisplayImageInfo(itkShortImagePointer & inputImage)
{
  itkShort2FloatImageCastFilter::Pointer castFilter = itkShort2FloatImageCastFilter::New();
  castFilter->SetInput(inputImage);
  castFilter->Update();
  itkFloatImagePointer tmpImage = castFilter->GetOutput();
  DisplayImageInfo(tmpImage);
}

void PandoraBoxImageFilters::DisplayImageInfo(itkFloatImagePointer & inputImage)
{
  std::cout<<"Image size : "<<inputImage->GetLargestPossibleRegion().GetSize()<<std::endl;
  std::cout<<"Image spacing : "<<inputImage->GetSpacing()<<std::endl;

  itkImageCalculatorFilter::Pointer calculatorFilter = itkImageCalculatorFilter::New();
  calculatorFilter->SetImage(inputImage);
  calculatorFilter->Compute();
  std::cout<<"Intensity Max : "<<calculatorFilter->GetMaximum()<<std::endl;
  std::cout<<"Intensity min : "<<calculatorFilter->GetMinimum()<<std::endl;

}







} // namespace btk
