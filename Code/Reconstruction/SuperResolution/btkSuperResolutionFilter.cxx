#include "btkSuperResolutionFilter.h"
#include "btkImageHelper.h"




namespace btk
{
//-----------------------------------------------------------------------------------------------------------
SuperResolutionFilter::SuperResolutionFilter()
{

}
//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::Initialize()
{
    m_H_Filter = H_Filter::New();

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
        m_Regions[i] = m_Images[i]->GetLargestPossibleRegion();

        m_SimulatedImages[i] = ImageType::New();

        m_SimulatedImages[i] = btk::ImageHelper< ImageType >::CreateNewImageFromPhysicalSpaceOf(m_Images[i]);

    }
}

//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::Update()
{

    m_PSF = btk::GaussianPSF::New();

    //m_PSF = btk::BoxCarPSF::New();

    m_H_Filter->SetImages(m_Images);

    m_H_Filter->SetTransforms(m_Transforms);

    m_H_Filter->SetReferenceImage(m_ReferenceImage);

    m_H_Filter->SetMasks(m_MasksObject);

    m_H_Filter->SetPSF(m_PSF);

    m_H_Filter->Update();

    m_H = m_H_Filter->GetH();

    m_Y = m_H_Filter->GetY();

    std::cout<<"H computed !"<<std::endl;

    // TODO: Minimization of m_X
    // such as Min(f(y - H*x) + lambda g(x))



    VNLCostFunction CostFunction = VNLCostFunction(m_X.size());

    CostFunction.GetCostFunction()->SetH(m_H);
    CostFunction.GetCostFunction()->SetLambda(10);
    CostFunction.GetCostFunction()->SetY(m_Y);
    CostFunction.GetCostFunction()->SetSRSize(m_ReferenceImage->GetLargestPossibleRegion().GetSize());
    vnl_vector< PrecisionType > HtY;
    m_H.mult(m_Y,HtY);
    CostFunction.GetCostFunction()->SetHtY(HtY);


    vnl_conjugate_gradient optimizer(CostFunction);
    optimizer.set_max_function_evals(1000);

    // Start minimization

    std::cout<<"Start minimization... "<<std::endl;


    optimizer.set_verbose(true);
    optimizer.minimize(m_X);
    optimizer.diagnose_outcome();

    m_Xfloat = vnl_matops::d2f(m_X);

    this->GenerateOutputData();


    if(m_ComputeSimulations)
    {
        this->SimulateLRImages();
    }


    // Generate the output





}

//-----------------------------------------------------------------------------------------------------------
void SuperResolutionFilter::SimulateLRImages()
{
    // H should previoulsy be computed
    vnl_vector< PrecisionType > simY;

    std::cout<<m_H.cols()<<std::endl;
    std::cout<<m_H.rows()<<std::endl;

    std::cout<<m_X.size()<<std::endl;



    m_H.mult(m_Xfloat,simY);





    unsigned int offset = 0;

    for (unsigned int im = 0; im < m_NumberOfImages; im++)
    {
      ImageType::IndexType absIndex;

      ImageType::IndexType start = m_Regions[im].GetIndex();
      ImageType::SizeType  size  = m_Regions[im].GetSize();
      unsigned int nvoxels = m_Regions[im].GetNumberOfPixels();
      ImageType::IndexType diffIndex;

      for( unsigned int i=0; i<nvoxels; i++)
      {

        diffIndex[2] = i / (size[0]*size[1]);

        diffIndex[1] = i - diffIndex[2]*size[0]*size[1];
        diffIndex[1] = diffIndex[1] / size[0];

        diffIndex[0] = i - diffIndex[2]*size[0]*size[1] - diffIndex[1]*size[0];

        absIndex[0] = diffIndex[0] + start[0];
        absIndex[1] = diffIndex[1] + start[1];
        absIndex[2] = diffIndex[2] + start[2];

        m_SimulatedImages[im]->SetPixel(absIndex,simY[i + offset]);

      }

      offset = offset + nvoxels;

    }

    //m_X.clear();
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


    //ENH: If we iterate over output image and we check the value of the current index
    // in m_x(doing the inverse conversion), it may be faster.

    for (unsigned int i = 0; i<m_X.size(); i++)
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


    btk::ImageHelper< ImageType >::WriteImage(m_Output,"test.nii.gz");

    //m_X.clear();
}


}
