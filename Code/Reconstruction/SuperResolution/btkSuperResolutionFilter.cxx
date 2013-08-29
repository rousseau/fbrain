#include "btkSuperResolutionFilter.h"
#include "btkImageHelper.h"




namespace btk
{
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter():m_Lambda(0.01),m_ComputeSimulations(false)
{
    m_H =NULL;
    m_Y = NULL;
    m_H_Filter = NULL;
    m_H_Filter = H_Filter::New();

}
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::~SuperResolutionFilter()
{
    delete m_H;
    delete m_Y;

    m_H =NULL;
    m_Y = NULL;
}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::Initialize()
{

    m_H = new vnl_sparse_matrix< PrecisionType >();
    m_Y = new vnl_vector< PrecisionType >();


    m_NumberOfImages  = m_Images.size();

    m_ReferenceRegion = this->m_ReferenceImage->GetLargestPossibleRegion();
    Iterator RefIt(this->m_ReferenceImage, m_ReferenceRegion);

    m_X.set_size(m_ReferenceRegion.GetNumberOfPixels() );

    unsigned int linearSize = 0;

    for(RefIt.GoToBegin();  !RefIt.IsAtEnd(); ++RefIt, linearSize++)
    {
        m_X[linearSize] = RefIt.Get();
    }

    m_MasksObject.resize(m_NumberOfImages);
    m_Regions.resize(m_NumberOfImages);
    m_SimulatedImages.resize(m_NumberOfImages);
    for(unsigned int i = 0; i< m_NumberOfImages; i++)
    {
        m_MasksObject[i] = MaskType::New();
        m_MasksObject[i]->SetImage(m_Masks[i]);
       // m_Regions[i] = m_Images[i]->GetLargestPossibleRegion();
        m_Regions[i] = m_MasksObject[i]->GetAxisAlignedBoundingBoxRegion();

        m_SimulatedImages[i] = ImageType::New();

        m_SimulatedImages[i] = btk::ImageHelper< ImageType >::CreateNewImageFromPhysicalSpaceOf(m_Images[i]);

    }
}

//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::Update()
{

    // H computation
    m_H_Filter->SetImages(m_Images);

    m_H_Filter->SetH(m_H);

    m_H_Filter->SetY(m_Y);

    m_H_Filter->SetTransforms(m_Transforms);

    m_H_Filter->SetInverseTransforms(m_InverseTransforms);

    m_H_Filter->SetReferenceImage(m_ReferenceImage);

    m_H_Filter->SetMasks(m_MasksObject);

    m_H_Filter->Update();


    // such as Min(f(y - H*x) + lambda g(x))
    vnl_vector< PrecisionType > HtY;
    // Premult H with Y
    //Since m_H is a pointer to H matrix H_Filter don't need to return it !!
    m_H->pre_mult(*m_Y,HtY);


    // Cost Function
    VNLCostFunction CostFunction = VNLCostFunction(m_X.size());

    CostFunction.GetCostFunction()->SetH(*m_H);//Set H
    CostFunction.GetCostFunction()->SetLambda(m_Lambda);
    CostFunction.GetCostFunction()->SetY(*m_Y);
    CostFunction.GetCostFunction()->SetSRSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());

    CostFunction.GetCostFunction()->SetHtY(HtY); //Set the precomputed HtY

   vnl_conjugate_gradient optimizer(CostFunction);
   optimizer.set_max_function_evals(20);

    // Start minimization

    std::cout<<"Start minimization... "<<std::endl;
    optimizer.set_verbose(true);
    optimizer.minimize(m_X);

    //display optimizer result
    optimizer.diagnose_outcome();

    // convert X into float (maybe not needed)
    m_Xfloat = vnl_matops::d2f(m_X);

    if(m_ComputeSimulations)
    {
        this->SimulateLRImages();
    }


     m_X.clear();

    // Generate the output

    this->GenerateOutputData();

   // clear all
    m_H->clear();
    m_Y->clear();
    HtY.clear();



}

//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::SimulateLRImages()
{
    // H should previoulsy be computed
    vnl_vector< PrecisionType > simY;

    m_H->mult(m_Xfloat,simY);

    //Temporary variables

    ImageType::IndexType lrIndex;  //index of the current voxel in the LR image
    unsigned int lrLinearIndex = 0;


    unsigned int offset = 0;
    for (unsigned int im = 0; im < m_NumberOfImages; im++)
    {
        ImageType::SizeType lrSize = m_SimulatedImages[im]->GetLargestPossibleRegion().GetSize();
        m_SimulatedImages[im]->FillBuffer(0);
        //itkIteratorWithIndex itSim(m_SimulatedImages[im],m_Regions[im]);
        itkIteratorWithIndex itSim(m_SimulatedImages[im], m_SimulatedImages[im]->GetLargestPossibleRegion());
        for(itSim.GoToBegin(); !itSim.IsAtEnd(); ++itSim)
        {
            lrIndex = itSim.GetIndex();

            lrLinearIndex = offset + lrIndex[0] +lrIndex[1]*lrSize[0]+ lrIndex[2]*lrSize[0]*lrSize[1];

            itSim.Set(simY[lrLinearIndex]);
        }
        offset += m_SimulatedImages[im]->GetLargestPossibleRegion().GetNumberOfPixels();
    }


}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::GenerateOutputData()
{

    m_Output = ImageType::New();

    // Allocate data
    ImageType::IndexType outputStart;
    outputStart[0] = 0; outputStart[1] = 0; outputStart[2] = 0;

    const ImageType * referenceImage = this->GetReferenceImage();

    ImageType::SizeType outputSize = referenceImage -> GetLargestPossibleRegion().GetSize();

    ImageType::RegionType outputRegion;
    outputRegion.SetIndex(outputStart);
    outputRegion.SetSize(outputSize);

    m_Output -> SetRegions(outputRegion);
    m_Output -> Allocate();
    m_Output -> FillBuffer((PixelType) 0.);

    m_Output -> SetOrigin( referenceImage -> GetOrigin() );
    m_Output -> SetSpacing( referenceImage -> GetSpacing() );
    m_Output -> SetDirection( referenceImage -> GetDirection() );

    ImageType::IndexType hrIndex;
    ImageType::IndexType hrStart = outputRegion.GetIndex();
    ImageType::SizeType  hrSize  = outputRegion.GetSize();


    //TODO: If we iterate over output image and we check the value of the current index
    // in m_x(doing the inverse conversion), it may be faster.

    for (unsigned int i = 0; i<m_Xfloat.size(); i++)
    {
      ImageType::IndexType hrDiffIndex;
      hrDiffIndex[2] = i / (hrSize[0]*hrSize[1]);

      hrDiffIndex[1] = i - hrDiffIndex[2]*hrSize[0]*hrSize[1];
      hrDiffIndex[1] = hrDiffIndex[1] / hrSize[0];

      hrDiffIndex[0] = i - hrDiffIndex[2]*hrSize[0]*hrSize[1] - hrDiffIndex[1]*hrSize[0];


      hrIndex[0] = hrDiffIndex[0] + hrStart[0];
      hrIndex[1] = hrDiffIndex[1] + hrStart[1];
      hrIndex[2] = hrDiffIndex[2] + hrStart[2];

      m_Output -> SetPixel(hrIndex, m_Xfloat[i] );

    }


   //btk::ImageHelper< ImageType >::WriteImage(m_Output,"test.nii.gz");

    //m_X.clear();
}


}
