#ifndef __BTK_WARPTRANSFORMTOIMAGEFILTER_TXX__
#define __BTK_WARPTRANSFORMTOIMAGEFILTER_TXX__

#include "btkWarpTransformToImageFilter.h"


namespace btk
{
template<typename TImageIn, typename TImageOut>
WarpTransformToImageFilter<TImageIn,TImageOut>::WarpTransformToImageFilter()
{
    m_Resampler = Resampler::New();
    m_Transform = itkTransform::New();

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
WarpTransformToImageFilter<TImageIn,TImageOut>::~WarpTransformToImageFilter()
{
    m_Resampler = NULL;
    m_Transform = NULL;
    m_OutputImage = NULL;
    m_InputImage = NULL;

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void WarpTransformToImageFilter<TImageIn,TImageOut>::Initialize()
{
    m_Resampler->SetInput(m_InputImage);
    m_Resampler->SetTransform(m_Transform);
    m_Resampler->SetOutputOrigin(m_InputImage->GetOrigin());
    m_Resampler->SetOutputDirection(m_InputImage->GetDirection());
    m_Resampler->SetOutputSpacing(m_InputImage->GetSpacing());
    m_Resampler->SetSize(m_InputImage->GetLargestPossibleRegion().GetSize());
}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void WarpTransformToImageFilter<TImageIn,TImageOut>::Update()
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
