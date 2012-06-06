

#include "btkMotionCorrectionSliceBySliceEulerFilter.h"

namespace btk
{

MotionCorrectionSliceBySliceEulerFilter::MotionCorrectionSliceBySliceEulerFilter()
{
    btkCoutMacro(MotionCorrectionSliceBySliceEulerFilter : Constructor );
}
//-----------------------------------------------------------------------------------------------------------

 MotionCorrectionSliceBySliceEulerFilter::~MotionCorrectionSliceBySliceEulerFilter()
{
    btkCoutMacro(MotionCorrectionSliceBySliceEulerFilter : Destructor );
}
//-----------------------------------------------------------------------------------------------------------

void  MotionCorrectionSliceBySliceEulerFilter::Update()
{
    btkCoutMacro(MotionCorrectionSliceBySliceEulerFilter : Update Method );
    if(SuperClass::m_ImagesLR.empty() || SuperClass::m_ImagesMaskLR.empty() || SuperClass::m_ReferenceImage.IsNull() || SuperClass::m_TransformsLR.empty())
    {
        throw std::string("Missing input arguments (LR images, LR Masks, Ref images, LR Transfos) ");
    }

    this->Initialize();
    this->DoRegistration();

}
//-----------------------------------------------------------------------------------------------------------

void  MotionCorrectionSliceBySliceEulerFilter::Initialize()
{
    m_SliceBySliceRegistration.resize(SuperClass::m_ImagesLR.size());
    m_OutputTransformsLR.resize(SuperClass::m_TransformsLR.size());
}
//-----------------------------------------------------------------------------------------------------------

void  MotionCorrectionSliceBySliceEulerFilter::DoRegistration()
{

    unsigned int im = 0;
    #pragma omp parallel for private(im) schedule(dynamic)

    for(im = 0; im< SuperClass::m_ImagesLR.size(); im++)
    {
        m_SliceBySliceRegistration[im] = SliceBySliceRegistration::New();
        m_SliceBySliceRegistration[im]->SetFixedImage(SuperClass::m_ImagesLR[im]);
        m_SliceBySliceRegistration[im]->SetMovingImage(SuperClass::m_ReferenceImage);
        m_SliceBySliceRegistration[im]->SetImageMask(SuperClass::m_ImagesMaskLR[im]);
        btkSliceBySliceTransformBase* transfo = reinterpret_cast<btkSliceBySliceTransformBase*>(m_TransformsLR[im].GetPointer());
        m_SliceBySliceRegistration[im]->SetTransform(transfo);

        try
        {
            m_SliceBySliceRegistration[im]->StartRegistration();
        }
        catch(itk::ExceptionObject & err)
        {
            throw err;
        }

        m_OutputTransformsLR[im] =  m_SliceBySliceRegistration[im]->GetTransform();


    }
}
//-----------------------------------------------------------------------------------------------------------

}


