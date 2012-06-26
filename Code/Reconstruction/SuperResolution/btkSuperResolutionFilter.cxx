#include "btkSuperResolutionFilter.h"



namespace btk
{
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter()
{
    m_PSFEstimationFilter = NULL;
    m_MotionCorrectionFilter = NULL;
    m_BiasCorrectionFilter = NULL;
    m_HighResolutionReconstructionFilter = NULL;
    m_SliceRejectionFilter = NULL;

    m_ReconstructionType = SR;

    m_TransformationType = AFFINE;

    m_UseMotionCorrection = true;

    btkCoutMacro(SuperResolutionFilter : Constructor );

    m_PSFEstimationFilter = new btk::PSFEstimationFilter();
    m_SliceRejectionFilter = new btk::SliceRejectionFilter();
    m_BiasCorrectionFilter = new btk::BiasCorrectionFilter();




}
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter(int loop, float beta, int nlm)
{
    m_PSFEstimationFilter = new btk::PSFEstimationFilter();
    m_SliceRejectionFilter = new btk::SliceRejectionFilter();
    m_BiasCorrectionFilter = new btk::BiasCorrectionFilter();

}
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter(int loop, float beta, int nlm, TRANSFORMATION_TYPE transfoType)
{
    m_PSFEstimationFilter = new btk::PSFEstimationFilter();
    m_SliceRejectionFilter = new btk::SliceRejectionFilter();
    m_BiasCorrectionFilter = new btk::BiasCorrectionFilter();

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
int SuperResolutionFilter::Update()
{
    btkCoutMacro(SuperResolutionFilter : Update Method );

    if(m_ImagesLR.empty() || m_ImagesMaskLR.empty() )
    {
        throw std::string("No Low Resolution Images or Masks ! Super Resolution can not be performed !");
    }

    //

    this->Initialize();

    if(m_ImageHR.IsNull())
    {
        this->ComputeHRImage();
        std::string name("tmp_reconstruct3D.nii.gz");
        btk::ImageHelper< itkImage >::WriteImage(m_ImageHR,name);
    }

    if(m_UseMotionCorrection)
    {
        this->MotionCorrection();
    }

    this->InverseTransforms();


     m_HighResolutionReconstructionFilter = new btk::HighResolutionSRFilter;
     m_HighResolutionReconstructionFilter->SetTransformsLR(m_TransformsLR);
     std::vector<itkTransformBase*> transfos;
     transfos.resize(m_TransformsLR.size());

     for(int i =0; i<transfos.size(); i++)
     {
         switch(m_TransformationType)
         {

             case SLICE_BY_SLICE :

                 break;

             case SLICE_BY_SLICE_AFFINE:

                 break;

             case SLICE_BY_SLICE_EULER:


                 break;

             case EULER_3D:


                 break;

             case AFFINE:

                 break;

             default:
                 break;


         }
     }

     m_HighResolutionReconstructionFilter->SetPsfType(m_Psftype);
     m_HighResolutionReconstructionFilter->SetNloops(m_Loop);
     m_HighResolutionReconstructionFilter->SetImagesLR(m_ImagesLR);
     m_HighResolutionReconstructionFilter->SetMasksLR(m_MasksLR);
     m_HighResolutionReconstructionFilter->SetImagesMaskLR(m_ImagesMaskLR);
     m_HighResolutionReconstructionFilter->SetImageHR(m_ImageHR);


     m_HighResolutionReconstructionFilter->SetTransformType(m_TransformationType);



     std::cout<<"Performing Super-Resolution with the previous computed HR image..."<<std::endl;

     m_HighResolutionReconstructionFilter->Update();

     std::cout<<"Done !"<<std::endl;




    //m_OutputHRImage = m_ImageHR;//m_HighResolutionReconstructionFilter->GetOutput();
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
    m_InterpolationOrderPSF = 1;
    m_InterpolationOrderIBP = 5;
    m_IterMax = 5;
    m_ComputeRegistration = true;
    m_Epsilon = 1e-4;
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::SetParameters(int Nlm, float Beta, int Loop, int MedianIBP, int PsfType,float lambda,int iterMax, int InterpolationOrderIBP, int InterpolationOrderPSF)
{
    m_Nlm = Nlm;
    m_Beta = Beta;
    m_Loop = Loop;
    m_MedianIBP = MedianIBP;
    m_Psftype = PsfType;
    m_Lambda = lambda;
    m_InterpolationOrderPSF = InterpolationOrderPSF;
    m_InterpolationOrderIBP = InterpolationOrderIBP;
    m_IterMax = iterMax;
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::InverseTransforms()
{


}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::Initialize()
{
    //TODO: Analyse the inputs for knowing what we should compute !
    //ex: if transfos are set we don't have to compute registration...
    switch(m_TransformationType)
    {
        case SLICE_BY_SLICE :
            m_MotionCorrectionFilter = new btk::MotionCorrectionSliceBySliceAffineFilter();
            break;

        case SLICE_BY_SLICE_AFFINE:
            m_MotionCorrectionFilter = new btk::MotionCorrectionSliceBySliceAffineFilter();
            break;

        case SLICE_BY_SLICE_EULER:
            m_MotionCorrectionFilter = new btk::MotionCorrectionSliceBySliceEulerFilter();
            break;

        case EULER_3D:
            m_MotionCorrectionFilter = new btk::MotionCorrection3DEulerFilter();
            break;

        case AFFINE:
            m_MotionCorrectionFilter = new btk::MotionCorrection3DAffineFilter();
            break;

        default:
            btkException("Wrong type of Transformation !")
            break;
    }



    if(m_ImageHR.IsNotNull())
    {
        m_ComputeRegistration = false;
    }
    else
    {
        m_ComputeRegistration = true;
    }


    m_MasksLR.resize(m_ImagesLR.size());

    if(m_RoisLR.empty())
    {
        m_RoisLR.resize(m_ImagesLR.size());
        for(int i = 0; i < m_ImagesMaskLR.size(); i++ )
        {
            m_MasksLR[i] = itkMask::New();
            m_MasksLR[i]->SetImage(m_ImagesMaskLR[i]);
            m_RoisLR[i] = m_MasksLR[i]->GetAxisAlignedBoundingBoxRegion();
        }
    }

}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::ComputeHRImage()
{

    if(m_ComputeRegistration)
    {

        std::cout<<"Compute HR Image with the "<< m_ImagesLR.size() <<" LR Images "<<std::endl;

        int numberOfImages = m_ImagesLR.size();

        if(m_TransformationType == SLICE_BY_SLICE ||
           m_TransformationType == SLICE_BY_SLICE_EULER ||
           m_TransformationType == EULER_3D)
        {
            m_CreateRigidHighResolutionImage = btk::LowToHighResFilterRigid::New();
            m_CreateRigidHighResolutionImage -> SetNumberOfImages(numberOfImages);
            m_CreateRigidHighResolutionImage -> SetTargetImage( 0 );
            m_CreateRigidHighResolutionImage -> SetMargin( 0 );

            for(int i = 0; i < m_ImagesLR.size(); i++)
            {
                m_CreateRigidHighResolutionImage->SetImageArray(i, m_ImagesLR[i]);
                m_CreateRigidHighResolutionImage->SetRegionArray(i, m_RoisLR[i]);
                m_CreateRigidHighResolutionImage->SetImageMaskArray(i, m_ImagesMaskLR[i]);

            }

            try
            {
                m_CreateRigidHighResolutionImage->StartRegistration();
            }
            catch(itk::ExceptionObject & excpt)
            {
                throw excpt;
            }

            m_ImageHR = m_CreateRigidHighResolutionImage->GetHighResolutionImage();
            m_ImageMaskHR = m_CreateRigidHighResolutionImage->GetImageMaskCombination();
            m_TransformsLR.resize(m_ImagesLR.size());
            m_TransformsLRSbS.resize(m_ImagesLR.size());
        }
        else
        {
            m_CreateAffineHighResolutionImage = btk::LowToHighResFilterAffine::New();
            m_CreateAffineHighResolutionImage -> SetNumberOfImages(numberOfImages);
            m_CreateAffineHighResolutionImage -> SetTargetImage( 0 );
            m_CreateAffineHighResolutionImage -> SetMargin( 0 );

            for(int i = 0; i < m_ImagesLR.size(); i++)
            {
                m_CreateAffineHighResolutionImage->SetImageArray(i, m_ImagesLR[i]);
                m_CreateAffineHighResolutionImage->SetRegionArray(i, m_RoisLR[i]);
                m_CreateAffineHighResolutionImage->SetImageMaskArray(i, m_ImagesMaskLR[i]);

            }

            try
            {
                m_CreateAffineHighResolutionImage->StartRegistration();
            }
            catch(itk::ExceptionObject & excpt)
            {
                throw excpt;
            }

            std::cout<<" Done. "<<std::endl;
            m_ImageHR = m_CreateAffineHighResolutionImage->GetHighResolutionImage();
            m_ImageMaskHR = m_CreateAffineHighResolutionImage->GetImageMaskCombination();
            m_TransformsLR.resize(m_ImagesLR.size());
            m_TransformsLRSbS.resize(m_ImagesLR.size());
        }


        for(unsigned int i=0; i<m_ImagesLR.size(); i++)
        {

            //TODO : New method, use a vector of itkTransform
            // Add .GetPointer() for casting

            switch(m_TransformationType)
            {
                default:
                    break;

                case SLICE_BY_SLICE_EULER:
                    m_TransformsLR[i] = btkEulerSliceBySliceTransform::New();

                    dynamic_cast<btkEulerSliceBySliceTransform*>(m_TransformsLR[i].GetPointer())->SetImage(m_ImagesLR[i]);
                    //TODO: Inverse ?
                    dynamic_cast<btkEulerSliceBySliceTransform*>(m_TransformsLR[i].GetPointer())->Initialize(m_CreateRigidHighResolutionImage->GetInverseTransformArray(i));
                    break;

                case SLICE_BY_SLICE:
                    m_TransformsLR[i] = btkEulerSliceBySliceTransform::New();
                    dynamic_cast<btkEulerSliceBySliceTransform*>(m_TransformsLR[i].GetPointer())->SetImage(m_ImagesLR[i]);
                    //TODO: Inverse ?
                    dynamic_cast<btkEulerSliceBySliceTransform*>(m_TransformsLR[i].GetPointer())->Initialize(m_CreateRigidHighResolutionImage->GetInverseTransformArray(i));
                    break;


                case SLICE_BY_SLICE_AFFINE:
                    m_TransformsLR[i] = btkAffineSliceBySliceTransform::New();
                    dynamic_cast<btkAffineSliceBySliceTransform*>(m_TransformsLR[i].GetPointer())->SetImage(m_ImagesLR[i]);
                    //m_TransformsLR[i]->SetImage(m_ImagesLR[i]);
                    //TODO: Inverse ?
                    dynamic_cast<btkAffineSliceBySliceTransform*>(m_TransformsLR[i].GetPointer())->Initialize(m_CreateAffineHighResolutionImage->GetInverseTransformArray(i));
                    break;


                case EULER_3D:
                    m_TransformsLR[i] = itkEulerTransform::New();
                    m_TransformsLR[i] = m_CreateRigidHighResolutionImage->GetTransformArray(i);
                    break;

                case AFFINE:
                    m_TransformsLR[i] = itkAffineTransform::New();
                    m_TransformsLR[i] = m_CreateAffineHighResolutionImage->GetTransformArray(i);
                    break;

            }

        }

    }
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::MotionCorrection()
{

    //TODO: Simply this
    itkImage::Pointer HRold = itkImage::New();
    itkImage::Pointer HRini = m_ImageHR;

    float previousMetric = 0.0;
    float currentMetric = 0.0;

    m_MotionCorrectionFilter->SetReferenceImage(m_ImageHR);
    m_MotionCorrectionFilter->SetImagesLR(m_ImagesLR);
    m_MotionCorrectionFilter->SetMasksLR(m_MasksLR);
    m_MotionCorrectionFilter->SetImageMaskHR(m_ImageMaskHR);
    m_MotionCorrectionFilter->SetImagesMaskLR(m_ImagesMaskLR);
    m_MotionCorrectionFilter->SetTransformsLR(m_TransformsLR);

    for(unsigned int it=1; it<= m_IterMax; it++)
    {

        std::cout << "Iteration : " << it <<" / "<<m_IterMax<< std::endl; std::cout.flush();


        m_MotionCorrectionFilter->Update();

        // Inject images

        std::cout << "Injecting images ... "; std::cout.flush();

        m_ResampleByInjectionFilter = ResamplerByInjectionType::New();

        for(unsigned int i = 0; i< m_ImagesLR.size(); i++)
        {
            m_ResampleByInjectionFilter->AddInput(m_ImagesLR[i]);
            m_ResampleByInjectionFilter->AddRegion(m_RoisLR[i]);

            //NOTE : Since ResampleByInjection use only SliceBySliceTransform we should do that this way :
            if(m_TransformationType == AFFINE)
            {
                btkAffineSliceBySliceTransform::Pointer transfo = btkAffineSliceBySliceTransform::New();
                transfo->SetImage(m_ImagesLR[i]);
                transfo->Initialize(static_cast< itkAffineTransform* >(m_TransformsLR[i].GetPointer()));
                m_ResampleByInjectionFilter->SetTransform(i, transfo);
            }
            else if(m_TransformationType == EULER_3D)
            {
                btkEulerSliceBySliceTransform::Pointer transfo = btkEulerSliceBySliceTransform::New();
                transfo->SetImage(m_ImagesLR[i]);
                transfo->Initialize(static_cast< itkEulerTransform* >(m_TransformsLR[i].GetPointer()));
                m_ResampleByInjectionFilter->SetTransform(i, transfo);
            }
            else
            {
                m_ResampleByInjectionFilter->SetTransform(i,static_cast<btkSliceBySliceTransformBase*>(m_TransformsLR[i].GetPointer()));
            }

        }

        m_ResampleByInjectionFilter->UseReferenceImageOn();
        m_ResampleByInjectionFilter->SetReferenceImage(m_ImageHR);
        m_ResampleByInjectionFilter->SetImageMask(m_ImageMaskHR);
        m_ResampleByInjectionFilter->Update();

        if(it == 1)
        {
            HRold = HRini;
        }
        else
        {
            HRold = m_ImageHR;
        }

        m_ImageHR = m_ResampleByInjectionFilter->GetOutput();

        std::cout<<" Done. "<<std::endl;

        //compute error

        typedef itk::Euler3DTransform< double > EulerTransformType;
        EulerTransformType::Pointer identity = EulerTransformType::New();
        identity -> SetIdentity();

        typedef itk::LinearInterpolateImageFunction<itkImage,double>     InterpolatorType;
        typedef itk::NormalizedCorrelationImageToImageMetric< itkImage,itkImage > NCMetricType;

        InterpolatorType::Pointer interpolator = InterpolatorType::New();

        NCMetricType::Pointer nc = NCMetricType::New();
        nc -> SetFixedImage(  HRold );
        nc -> SetMovingImage( m_ImageHR );
        nc -> SetFixedImageRegion( HRold -> GetLargestPossibleRegion() );
        nc -> SetTransform( identity );
        nc -> SetInterpolator( interpolator );
        nc -> Initialize();

        previousMetric = currentMetric;
        currentMetric = - nc -> GetValue( identity -> GetParameters() );
        std::cout<<"previousMetric: "<<previousMetric<<", currentMetric: "<<currentMetric<<"\n";
        double delta = 0.0;

        if (it >= 2)
          delta = (currentMetric - previousMetric) / previousMetric;
        else
          delta = 1;

        if (delta < m_Epsilon) break;


    }










}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::ResampleByInjection()
{
    m_ResampleByInjectionFilter = ResamplerByInjectionType::New();



}

//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------



}
