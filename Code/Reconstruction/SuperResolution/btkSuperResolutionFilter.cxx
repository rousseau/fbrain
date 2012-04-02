#include "btkSuperResolutionFilter.h"

#include "btkHighResolutionIBPFilter.h"
#include "btkHighResolutionSRFilter.h"
#include "btkMotionCorrectionAffine3DFilter.h"
#include "btkMotionCorrectionSliceBySliceFilter.h"

namespace btk
{

SuperResolutionFilter::SuperResolutionFilter()
{
    m_PSFEstimationFilter = NULL;
    m_MotionCorrectionFilter = NULL;
    m_BiasCorrectionFilter = NULL;
    m_HighResolutionReconstructionFilter = NULL;
    m_SliceRejectionFilter = NULL;

    m_ReconstructionType = IBP;

    m_TransformationType = AFFINE;

    btkCoutMacro(SuperResolutionFilter : Constructor );




}
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter(int loop, float beta, int nlm)
{

}
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter(int loop, float beta, int nlm, TRANSFORMATION_TYPE transfoType)
{
}
//-----------------------------------------------------------------------------------------------------------
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

    if(m_SliceRejectionFilter != NULL)
    {
        delete m_SliceRejectionFilter;
        m_SliceRejectionFilter = NULL;
    }
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::Update()
{
    btkCoutMacro(SuperResolutionFilter : Update Method );

    this->InverseTransforms();

    m_PSFEstimationFilter = new btk::PSFEstimationFilter();
    m_PSFEstimationFilter->Update();


    if(m_TransformationType == AFFINE)
    {
        m_MotionCorrectionFilter = new btk::MotionCorrectionAffine3DFilter();
    }
    else
    {
        m_MotionCorrectionFilter = new btk::MotionCorrectionSliceBySliceFilter();
    }

    m_MotionCorrectionFilter->Update();

    m_BiasCorrectionFilter = new btk::BiasCorrectionFilter();
    m_BiasCorrectionFilter->Update();




    if(m_ReconstructionType == IBP)
    {

         m_HighResolutionReconstructionFilter = new btk::HighResolutionIBPFilter();

         dynamic_cast< btk::HighResolutionIBPFilter* >(m_HighResolutionReconstructionFilter)->SetNloops(m_Loop);
         dynamic_cast< btk::HighResolutionIBPFilter* >(m_HighResolutionReconstructionFilter)->SetNlm(m_Nlm);
         dynamic_cast< btk::HighResolutionIBPFilter* >(m_HighResolutionReconstructionFilter)->SetBeta(m_Beta);
         dynamic_cast< btk::HighResolutionIBPFilter* >(m_HighResolutionReconstructionFilter)->SetMedianIBP(m_MedianIBP);
         m_HighResolutionReconstructionFilter->SetInterpolationOrderIBP(m_InterpolationOrderIBP);

    }
    else
    {

        m_HighResolutionReconstructionFilter = new btk::HighResolutionSRFilter();     
        m_HighResolutionReconstructionFilter->SetInterpolationOrderPSF(m_InterpolationOrderPSF);
    }

     m_HighResolutionReconstructionFilter->SetImagesLR(m_ImagesLR);
     m_HighResolutionReconstructionFilter->SetImagesMaskLR(m_ImagesMaskLR);
     m_HighResolutionReconstructionFilter->SetImageHR(m_ImageHR);


     if(m_TransformationType == AFFINE)
     {
         m_HighResolutionReconstructionFilter->SetTransformsLRAffine(m_TransformsLRAffine);
         m_HighResolutionReconstructionFilter->SetTransformType(AFFINE);
         m_HighResolutionReconstructionFilter->SetInverseTransformsLRAffine(m_InverseTransformsLRAffine);
     }
     else
     {
         m_HighResolutionReconstructionFilter->SetTransformsLRSbS(m_TransformsLRSbS);
         m_HighResolutionReconstructionFilter->SetTransformType(SLICE_BY_SLICE);
         m_HighResolutionReconstructionFilter->SetInverseTransformsLRSbS(m_InverseTransformsLRSbS);
     }


     m_HighResolutionReconstructionFilter->Update();



    m_SliceRejectionFilter = new btk::SliceRejectionFilter();
    m_SliceRejectionFilter->Update();

    m_OutputHRImage = m_HighResolutionReconstructionFilter->GetOutput();

}
//-----------------------------------------------------------------------------------------------------------
itkImage::Pointer SuperResolutionFilter::GetOutput()
{

    return m_OutputHRImage;
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::SetDefaultParameters()
{
    m_Nlm = 0;
    m_Beta = 1;
    m_Loop = 50;
    m_MedianIBP = 0;
    m_Psftype = 1;
    m_InterpolationOrderPSF = 5;
    m_InterpolationOrderIBP = 5;
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::SetParameters(int Nlm, float Beta, int Loop, int MedianIBP, int PsfType, int InterpolationOrderIBP, int InterpolationOrderPSF)
{
    m_Nlm = Nlm;
    m_Beta = Beta;
    m_Loop = Loop;
    m_MedianIBP = MedianIBP;
    m_Psftype = PsfType;
    m_InterpolationOrderPSF = InterpolationOrderPSF;
    m_InterpolationOrderIBP = InterpolationOrderIBP;
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::InverseTransforms()
{


    m_InverseTransformsLRAffine.resize(m_TransformsLRAffine.size());
    m_InverseTransformsLRSbS.resize(m_TransformsLRSbS.size());
    for(int i = 0; i < m_ImagesLR.size(); i++)
    {
        if(m_TransformationType == AFFINE)
        {
            m_TransformsLRAffine[i]->GetInverse( m_InverseTransformsLRAffine[i]);
        }
        else
        {
            m_TransformsLRSbS[i]->GetInverse(m_InverseTransformsLRSbS[i]);
        }


    }

}
//-----------------------------------------------------------------------------------------------------------


}
