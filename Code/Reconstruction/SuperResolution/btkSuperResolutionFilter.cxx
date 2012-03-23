#include "btkSuperResolutionFilter.h"

namespace btk
{

SuperResolutionFilter::SuperResolutionFilter()
{
    m_PSFEstimationFilter = NULL;
    m_MotionCorrectionFilter = NULL;
    m_BiasCorrectionFilter = NULL;
    m_HighResolutionReconstructionFilter = NULL;

    btkCoutMacro(SuperResolutionFilter : Constructor );




}

SuperResolutionFilter::SuperResolutionFilter(int loop, float beta, int nlm)
{

}

SuperResolutionFilter::SuperResolutionFilter(int loop, float beta, int nlm, SuperResolutionFilter::TRANSFORMATION_TYPE transfoType)
{
}

SuperResolutionFilter::~SuperResolutionFilter()
{
    btkCoutMacro(SuperResolutionFilter : Destructor );

    if(m_PSFEstimationFilter != NULL)
    {
        delete m_PSFEstimationFilter;
        m_PSFEstimationFilter = NULL;

    }

    if(m_MotionCorrectionFilter != NULL)
    {
        delete m_MotionCorrectionFilter;
        m_MotionCorrectionFilter = NULL;
    }

    if(m_BiasCorrectionFilter != NULL)
    {
        delete m_BiasCorrectionFilter;
        m_BiasCorrectionFilter = NULL;
    }

    if(m_HighResolutionReconstructionFilter != NULL)
    {
        delete m_HighResolutionReconstructionFilter;
        m_HighResolutionReconstructionFilter = NULL;
    }
}

void SuperResolutionFilter::Update()
{
    btkCoutMacro(SuperResolutionFilter : Update Method );

    m_PSFEstimationFilter = new btk::PSFEstimationFilter();
    m_PSFEstimationFilter->Update();

    m_MotionCorrectionFilter = new btk::MotionCorrectionFilter();
    m_MotionCorrectionFilter->Update();

    m_BiasCorrectionFilter = new btk::BiasCorrectionFilter();
    m_BiasCorrectionFilter->Update();

    m_HighResolutionReconstructionFilter = new btk::HighResolutionReconstructionFilter();
    m_HighResolutionReconstructionFilter->Update();


}

SuperResolutionFilter::itkImage::Pointer SuperResolutionFilter::GetOutput()
{

    return m_ImageHR;
}



}
