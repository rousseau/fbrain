#ifndef __BTK_APPLYTRANSFORMTOIMAGEFILTER_TXX__
#define __BTK_APPLYTRANSFORMTOIMAGEFILTER_TXX__

#include "btkApplyTransformToImageFilter.h"


namespace btk
{
template<typename TImageIn, typename TImageOut>
ApplyTransformToImageFilter<TImageIn,TImageOut>::ApplyTransformToImageFilter()
{
    m_OutputImage = itkImageOut::New();
}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
ApplyTransformToImageFilter<TImageIn,TImageOut>::~ApplyTransformToImageFilter()
{
    m_Transform = NULL;
    m_OutputImage = NULL;
    m_InputImage = NULL;

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Initialize()
{
    if(!m_Transform)
    {
        btkException("Missing input transform !");
    }
    if(!m_InputImage)
    {
        btkException("Missing input image !");
    }
    if(!m_ReferenceImage)
    {
        btkException("Missing reference image !");
    }

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Update()
{
    try
    {
        this->Resample();
    }
    catch(itk::ExceptionObject & excp)
    {
        std::cerr << "Error when apply transformation" << std::endl;
        std::cerr << excp << std::endl;
        std::cout << "[FAILED]" << std::endl;
        throw excp;
    }



}
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Resample() throw(itk::ExceptionObject &)
{
    typename Interpolator::Pointer interpolator = Interpolator::New();
    interpolator->SetInputImage(m_ReferenceImage);

    m_OutputImage =  btk::ImageHelper<itkImageOut>::CreateNewImageFromPhysicalSpaceOf(m_InputImage.GetPointer());

    IteratorIn it(m_InputImage, m_InputImage->GetLargestPossibleRegion());
    IteratorOut itO(m_OutputImage, m_OutputImage->GetLargestPossibleRegion());


    for(it.GoToBegin(), itO.GoToBegin(); !it.IsAtEnd(), !itO.IsAtEnd(); ++it, ++itO)
    {
        typename itkImage::IndexType index = it.GetIndex();
        typename Interpolator::ContinuousIndexType Tindex;


        typename itkImage::PointType point, Tpoint;
        m_InputImage->TransformIndexToPhysicalPoint(index, point);

        Tpoint = m_Transform->TransformPoint(point);


        if(interpolator->IsInsideBuffer(Tpoint))
        {
            m_ReferenceImage->TransformPhysicalPointToContinuousIndex(Tpoint,Tindex);

            typename itkImage::PixelType value = interpolator->EvaluateAtContinuousIndex(Tindex);

            itO.Set(value);
        }
        else
        {
            itO.Set(0);
        }

    }
}
}

#endif
