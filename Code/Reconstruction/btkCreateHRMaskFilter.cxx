#include "btkCreateHRMaskFilter.h"

namespace btk
{

CreateHRMaskFilter::CreateHRMaskFilter()
{
}

CreateHRMaskFilter::~CreateHRMaskFilter()
{
}

void CreateHRMaskFilter::Update()
{

    std::cout<<"Create Mask HR Image by interpolating Masks of LR Images\n";

    //Initialize to 0 the output HR image
    m_maskHRImage->FillBuffer(0);

    //parameters for interpolation (bspline interpolator)
    int interpolationOrder = 0;
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(interpolationOrder);

    for(uint i=0; i<m_inputLRImages.size(); i++)
    {
       const char * className =  m_transforms[i]->GetNameOfClass();
       // TODO: Add a different method for each type of transformation

      //interpolate the LR mask
      itkResampleFilter::Pointer resample = itkResampleFilter::New();
      resample->SetTransform(m_transforms[i]); //testing
      resample->SetInterpolator(bsInterpolator);
      resample->UseReferenceImageOn();
      resample->SetReferenceImage(m_HRImage);
      resample->SetInput(m_inputLRImages[i]);

      //Add the interpolated mask
      itkAddImageFilter::Pointer addFilter = itkAddImageFilter::New ();
      addFilter->SetInput1(m_maskHRImage);
      addFilter->SetInput2(resample->GetOutput());
      addFilter->Update();

      m_maskHRImage = addFilter->GetOutput();
    }

    itkStatisticsImageFilter::Pointer statisticsImageFilter = itkStatisticsImageFilter::New ();
    statisticsImageFilter->SetInput(m_maskHRImage);
    statisticsImageFilter->Update();
    std::cout << "Stat of the HR mask: \nMean: " << statisticsImageFilter->GetMean() << std::endl;
    std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
    std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
    std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;

    std::cout<<"Binarize the HR mask image\n";
    itkBinaryThresholdImageFilter::Pointer thresholdFilter = itkBinaryThresholdImageFilter::New();
    thresholdFilter->SetInput( m_maskHRImage );
    thresholdFilter->SetLowerThreshold(0.5);
    thresholdFilter->SetUpperThreshold( statisticsImageFilter->GetMaximum() + 1);
    thresholdFilter->SetInsideValue(1.0);
    thresholdFilter->SetOutsideValue(0.0);
    thresholdFilter->Update();
    m_maskHRImage = thresholdFilter->GetOutput();
}

}//end namespace
