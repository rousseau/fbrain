#ifndef __BTK_APPLYTRANSFORMTOIMAGEFILTER_TXX__
#define __BTK_APPLYTRANSFORMTOIMAGEFILTER_TXX__

#include "btkApplyTransformToImageFilter.h"


namespace btk
{
template<typename TImageIn, typename TImageOut>
ApplyTransformToImageFilter<TImageIn,TImageOut>::ApplyTransformToImageFilter()
{
    m_Resampler = Resampler::New();

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
ApplyTransformToImageFilter<TImageIn,TImageOut>::~ApplyTransformToImageFilter()
{
    m_Resampler = NULL;
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
        btkException("Missing input image !")
    }

    m_Resampler->SetInput(m_InputImage);
    m_Resampler->SetTransform(m_Transform);
    m_Resampler->SetOutputOrigin(m_InputImage->GetOrigin());
    m_Resampler->SetOutputDirection(m_InputImage->GetDirection());
    m_Resampler->SetOutputSpacing(m_InputImage->GetSpacing());
    m_Resampler->SetSize(m_InputImage->GetLargestPossibleRegion().GetSize());
}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Update()
{
    try
    {
        m_Resampler->Update();
    }
    catch(itk::ExceptionObject & excp)
    {
        std::cerr << "Error in the Warping" << std::endl;
        std::cerr << excp << std::endl;
        std::cout << "[FAILED]" << std::endl;
        throw excp;
    }

    this->m_OutputImage = m_Resampler->GetOutput();


}
}

#endif
