
#include "btkHighResolutionSRFilter.h"

namespace btk
{
HighResolutionSRFilter::HighResolutionSRFilter()
{
    btkCoutMacro(HighResolutionSRFilter : Constructor );
    //resampler = Resampler::New();
    m_NlmTools = new btkNLMTool<float>;
}
//-----------------------------------------------------------------------------------------------------------
HighResolutionSRFilter::~HighResolutionSRFilter()
{
    btkCoutMacro(HighResolutionSRFilter : Destructor );

    if(m_NlmTools != NULL)
    {
        delete m_NlmTools;
        m_NlmTools = NULL;
    }
}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionSRFilter::Update()
{
    btkCoutMacro(HighResolutionSRFilter : Update Method );


    std::cout<<"Performing super resolution"<<std::endl;

    if(m_UseAffineFilter && !m_UseEulerFilter)
    {
        this->DoAffineReconstruction();
    }
    else if(!m_UseAffineFilter && m_UseEulerFilter)
    {
        this->DoRigidReconstruction();
    }
    else
    {
        //...
    }

}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionSRFilter::Initialize()
{
    std::string NameOfTransforms;

    if(!m_TransformsLR.empty())
    {
        NameOfTransforms = m_TransformsLR[0]->GetNameOfClass();

        if(NameOfTransforms == "SliceBySliceTransform" ||
           NameOfTransforms == "EulerSliceBySliceTransform")
        {
            m_UseEulerFilter = true;
            m_UseAffineFilter = false;
            m_UseSliceBySlice = true;
        }
        else if(NameOfTransforms == "AffineSliceBySliceTransform" )
        {
            m_UseAffineFilter = true;
            m_UseEulerFilter = false;
            m_UseSliceBySlice = true;
        }
        else if(NameOfTransforms == "AffineTransform")
        {
                m_UseAffineFilter = true;
                m_UseEulerFilter = false;
                m_UseSliceBySlice = false;
        }
        else if(NameOfTransforms == "EulerTransform")
        {
                m_UseAffineFilter = false;
                m_UseEulerFilter = true;
                m_UseSliceBySlice = false;
        }

        else
        {
            btkException("Wrong type of Transformation !");
        }
    }
    else
    {
        btkException("Transforms are not set !");
    }

}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionSRFilter::DoAffineReconstruction()
{
    std::cout<<"Performing super resolution"<<std::endl;
    itkImage::SizeType roiSize;
    itkImage::IndexType roiIndex;

    typedef btk::SuperResolutionAffineImageFilter< itkImage, itkImage >   Resampler;
    Resampler::Pointer resampler = Resampler::New();



    for(unsigned int i = 0; i < SuperClass::m_ImagesLR.size(); i++)
    {

        resampler->AddInput(SuperClass::m_ImagesLR[i]);
        if(SuperClass::m_MasksLR.size() > 0)
        {
            resampler->AddMask(SuperClass::m_MasksLR[i]);
            itkImage::RegionType roi = SuperClass::m_MasksLR[i]->GetAxisAlignedBoundingBoxRegion();
            roiIndex = roi.GetIndex();
            roiSize = roi.GetSize();
        }
        else
        {
            roiIndex = SuperClass::m_ImagesLR[i]->GetLargestPossibleRegion().GetIndex();
            roiSize = SuperClass::m_ImagesLR[i]->GetLargestPossibleRegion().GetSize();
        }
        itkImage::RegionType imageRegion;
        imageRegion.SetIndex(roiIndex);
        imageRegion.SetSize (roiSize);
        resampler -> AddRegion( imageRegion );

        itkImage::SizeType size = SuperClass::m_ImagesLR[i]->GetLargestPossibleRegion().GetSize();

        for(unsigned int j=0; j< size[2]; j++)//SuperClass::m_TransformsLR[i]->GetNumberOfSlices(); j++)
        {
            if(m_UseSliceBySlice)
            {
                resampler -> SetTransform(i, j, dynamic_cast<btkAffineSliceBySliceTransform*>(SuperClass::m_TransformsLR[i].GetPointer())-> GetSliceTransform(j) ) ;
            }
            else
            {
                resampler -> SetTransform(i, j, dynamic_cast<itkAffineTransform*>(SuperClass::m_TransformsLR[i].GetPointer()) ) ;
            }
        }
    }

    resampler -> UseReferenceImageOn();
    resampler -> SetReferenceImage( SuperClass::m_ImageHR );
    resampler -> SetIterations(m_Iter);
    resampler -> SetLambda( m_Lambda );

    if ( SuperClass::m_PsfType == Resampler::BOXCAR )
    {
        resampler -> SetPSF( Resampler::BOXCAR );

    }
    else
    {
        resampler -> SetPSF( Resampler::GAUSSIAN );
    }

    resampler -> Update();


    for (int i=0; i< SuperClass::m_Nloops; i++)
    {
      std::cout<<"Loop : "<<i+1<<std::endl;

      m_NlmTools->SetInput(resampler -> GetOutput());
      m_NlmTools->SetPaddingValue(0);
      m_NlmTools->SetDefaultParameters();
      m_NlmTools->ComputeOutput();

      itkImage::Pointer outputImage = itkImage::New();
      m_NlmTools->GetOutput(outputImage);

      resampler -> SetReferenceImage( outputImage );
      resampler -> Update();
    }
    //NLM denoising desired at the last step if number of loops > 0
    if(SuperClass::m_Nloops>0)
    {

      m_NlmTools->SetInput(resampler -> GetOutput());
      m_NlmTools->SetPaddingValue(0);
      m_NlmTools->SetDefaultParameters();
      m_NlmTools->ComputeOutput();

      itkImage::Pointer outputImage = itkImage::New();
      m_NlmTools->GetOutput(outputImage);

      resampler -> SetReferenceImage( outputImage );
    }

    SuperClass::m_OutputHRImage = resampler->GetOutput();
}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionSRFilter::DoRigidReconstruction()
{
    std::cout<<"Performing super resolution"<<std::endl;
    itkImage::SizeType roiSize;
    itkImage::IndexType roiIndex;
    typedef btk::SuperResolutionRigidImageFilter< itkImage, itkImage >   Resampler;
    Resampler::Pointer resampler = Resampler::New();

    for(unsigned int i = 0; i < SuperClass::m_ImagesLR.size(); i++)
    {

        resampler->AddInput(SuperClass::m_ImagesLR[i]);
        if(SuperClass::m_MasksLR.size() > 0)
        {
            resampler->AddMask(SuperClass::m_MasksLR[i]);
            itkImage::RegionType roi = SuperClass::m_MasksLR[i]->GetAxisAlignedBoundingBoxRegion();
            roiIndex = roi.GetIndex();
            roiSize = roi.GetSize();
        }
        else
        {
            roiIndex = SuperClass::m_ImagesLR[i]->GetLargestPossibleRegion().GetIndex();
            roiSize = SuperClass::m_ImagesLR[i]->GetLargestPossibleRegion().GetSize();
        }
        itkImage::RegionType imageRegion;
        imageRegion.SetIndex(roiIndex);
        imageRegion.SetSize (roiSize);
        resampler -> AddRegion( imageRegion );

        itkImage::SizeType size = SuperClass::m_ImagesLR[i]->GetLargestPossibleRegion().GetSize();

        for(unsigned int j=0; j< size[2]; j++)//SuperClass::m_TransformsLR[i]->GetNumberOfSlices(); j++)
        {
            if(m_UseSliceBySlice)
            {
                resampler -> SetTransform(i, j, dynamic_cast<btkEulerSliceBySliceTransform*>(SuperClass::m_TransformsLR[i].GetPointer())-> GetSliceTransform(j) ) ;
            }
            else
            {
                resampler -> SetTransform(i, j, dynamic_cast<itkEulerTransform*>(SuperClass::m_TransformsLR[i].GetPointer()) ) ;
            }


        }
    }

    resampler -> UseReferenceImageOn();
    resampler -> SetReferenceImage( SuperClass::m_ImageHR );
    resampler -> SetIterations(m_Iter);
    resampler -> SetLambda( m_Lambda );

    if ( SuperClass::m_PsfType == Resampler::BOXCAR )
    {
        resampler -> SetPSF( Resampler::BOXCAR );

    }
    else
    {
        resampler -> SetPSF( Resampler::GAUSSIAN );
    }

    resampler -> Update();


    for (int i=0; i< SuperClass::m_Nloops; i++)
    {
      std::cout<<"Loop : "<<i+1<<std::endl;

      m_NlmTools->SetInput(resampler -> GetOutput());
      m_NlmTools->SetPaddingValue(0);
      m_NlmTools->SetDefaultParameters();
      m_NlmTools->ComputeOutput();

      itkImage::Pointer outputImage = itkImage::New();
      m_NlmTools->GetOutput(outputImage);

      resampler -> SetReferenceImage( outputImage );
      resampler -> Update();
    }
    //NLM denoising desired at the last step if number of loops > 0
    if(SuperClass::m_Nloops>0)
    {

      m_NlmTools->SetInput(resampler -> GetOutput());
      m_NlmTools->SetPaddingValue(0);
      m_NlmTools->SetDefaultParameters();
      m_NlmTools->ComputeOutput();

      itkImage::Pointer outputImage = itkImage::New();
      m_NlmTools->GetOutput(outputImage);

      resampler -> SetReferenceImage( outputImage );
    }

    SuperClass::m_OutputHRImage = resampler->GetOutput();
}

}

