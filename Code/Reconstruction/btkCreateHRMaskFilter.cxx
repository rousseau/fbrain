#include "btkCreateHRMaskFilter.h"

#include "itkImageDuplicator.h"

namespace btk
{

CreateHRMaskFilter::CreateHRMaskFilter()
{
    m_MaskHRImage = NULL;
}
//-----------------------------------------------------------------------------------------------------------
CreateHRMaskFilter::~CreateHRMaskFilter()
{
}
//-----------------------------------------------------------------------------------------------------------
void CreateHRMaskFilter::Update()
{


    itkDuplicator::Pointer duplicator3 = itkDuplicator::New();
    duplicator3->SetInputImage( m_HRImage );
    duplicator3->Update();
    m_MaskHRImage = duplicator3->GetOutput();


    std::cout<<"Create Mask HR Image by interpolating Masks of LR Images\n";
    //Initialize to 0 the output HR image
    m_MaskHRImage->FillBuffer(0);

    //parameters for interpolation (bspline interpolator)
    int interpolationOrder = 0;
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(interpolationOrder);

    for(unsigned int i=0; i<m_InputLRImages.size(); i++)
    {

       // TODO: Check if SbS transform do the correct Transformation

      //interpolate the LR mask
      itkResampleFilter::Pointer resample = itkResampleFilter::New();
      if(m_TransformType ==  AFFINE)
      {
          resample->SetTransform(m_TransformsAffine[i]);
      }
      else
      {
          resample->SetTransform(m_TransformsSbS[i]);
      }

      resample->SetInterpolator(bsInterpolator);
      resample->UseReferenceImageOn();
      resample->SetReferenceImage(m_HRImage);
      resample->SetInput(m_InputLRImages[i]);

      //Add the interpolated mask
      itkAddImageFilter::Pointer addFilter = itkAddImageFilter::New ();
      addFilter->SetInput1(m_MaskHRImage);
      addFilter->SetInput2(resample->GetOutput());
      addFilter->Update();

      m_MaskHRImage = addFilter->GetOutput();
    }

    itkStatisticsImageFilter::Pointer statisticsImageFilter = itkStatisticsImageFilter::New ();
    statisticsImageFilter->SetInput(m_MaskHRImage);
    statisticsImageFilter->Update();
    std::cout << "Stat of the HR mask: \nMean: " << statisticsImageFilter->GetMean() << std::endl;
    std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
    std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
    std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;

    std::cout<<"Binarize the HR mask image\n";
    itkBinaryThresholdImageFilter::Pointer thresholdFilter = itkBinaryThresholdImageFilter::New();
    thresholdFilter->SetInput( m_MaskHRImage );
    thresholdFilter->SetLowerThreshold(0.5);
    thresholdFilter->SetUpperThreshold( statisticsImageFilter->GetMaximum() + 1);
    thresholdFilter->SetInsideValue(1.0);
    thresholdFilter->SetOutsideValue(0.0);
    thresholdFilter->Update();
    m_MaskHRImage = thresholdFilter->GetOutput();
}

}//end namespace
