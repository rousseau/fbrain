

#include "btkMotionCorrection3DEulerFilter.h"

namespace btk
{

MotionCorrection3DEulerFilter::MotionCorrection3DEulerFilter()
{
    btkCoutMacro("MotionCorrection3DEulerFilter : Constructor");
    m_UseAffine = false;
    m_UseEuler = false;
}
//-----------------------------------------------------------------------------------------------------------

MotionCorrection3DEulerFilter::~MotionCorrection3DEulerFilter()
{
    btkCoutMacro("MotionCorrection3DEulerFilter : Destructor");
}
//-----------------------------------------------------------------------------------------------------------

void MotionCorrection3DEulerFilter::Update()
{
    btkCoutMacro("MotionCorrection3DEulerFilter : Update Method");
    if(SuperClass::m_ImagesLR.empty() || SuperClass::m_MasksLR.empty() || SuperClass::m_ReferenceImage.IsNull() || SuperClass::m_TransformsLR.empty())
    {
        std::cout<<"Images LR "<<SuperClass::m_ImagesLR.empty()<<std::endl;
        std::cout<<"Images masks "<<SuperClass::m_ImagesMaskLR.empty()<<std::endl;
        std::cout<<"Image Reference "<<SuperClass::m_ReferenceImage.IsNull()<<std::endl;
        std::cout<<"Transforms "<<SuperClass::m_TransformsLR.empty()<<std::endl;
        throw std::string("Missing input arguments (LR images, LR Masks, Ref images, LR Transfos) ");
    }

    this->Initialize();
    this->DoRegistration();




}
//-----------------------------------------------------------------------------------------------------------

void MotionCorrection3DEulerFilter::Initialize()
{
    m_Rigid3DRegistration.resize(SuperClass::m_ImagesLR.size());
    m_OutputTransformsLR.resize(m_TransformsLR.size());



}
//-----------------------------------------------------------------------------------------------------------

void MotionCorrection3DEulerFilter::DoRegistration()
{

    unsigned int im = 0;


    #pragma omp parallel for private(im) schedule(dynamic)

    for(im = 0; im< SuperClass::m_ImagesLR.size(); im++)
    {
        m_Rigid3DRegistration[im] = Rigid3DRegistration::New();
        m_Rigid3DRegistration[im]->SetFixedImage(SuperClass::m_ImagesLR[im]);
        m_Rigid3DRegistration[im]->SetMovingImage(SuperClass::m_ReferenceImage);
        m_Rigid3DRegistration[im]->SetFixedImageMask(SuperClass::m_ImagesMaskLR[im]);
        m_Rigid3DRegistration[im]->SetTransform(SuperClass::m_TransformsLR[im]);

        try
        {
            m_Rigid3DRegistration[im]->Update();
        }
        catch(itk::ExceptionObject & err)
        {
            throw err;
        }


        m_OutputTransformsLR[im] =  m_Rigid3DRegistration[im]->GetTransform();



    }
}
//-----------------------------------------------------------------------------------------------------------

}


