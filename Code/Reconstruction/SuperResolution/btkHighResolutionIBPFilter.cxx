#include "btkHighResolutionIBPFilter.h"

#include "itkBSplineInterpolationWeightFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "itkSubtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkStatisticsImageFilter.h"

#include "btkNLMTool.h"


namespace btk
{

HighResolutionIBPFilter::HighResolutionIBPFilter()
{
    btkCoutMacro(HighResolutionIBPFilter : Constructor );

    m_ImageHR = NULL;

    m_HRMaskFilter = NULL;
    m_SimuLRImagesFilter = NULL;
    m_PaddingValue = 0;
    m_HRMaskFilter = new btk::CreateHRMaskFilter();
    m_SimuLRImagesFilter = new btk::SimulateLRImageFilter();

}
//-----------------------------------------------------------------------------------------------------------
HighResolutionIBPFilter::~HighResolutionIBPFilter()
{
    btkCoutMacro(HighResolutionIBPFilter : Destructor );

    if(m_HRMaskFilter != NULL)
    {
        delete m_HRMaskFilter;
    }
    if(m_SimuLRImagesFilter != NULL)
    {
        delete m_SimuLRImagesFilter;
    }
}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionIBPFilter::Update()
{
    btkCoutMacro(HighResolutionIBPFilter : Update Method );

    this->Initialize();


    if(m_TransformType == AFFINE)
    {
        m_HRMaskFilter->SetTransformsAffine(m_TransformsLRAffine);
        m_HRMaskFilter->SetTransformType(m_TransformType);
    }
    else
    {
        m_HRMaskFilter->SetTransformsSbS(m_TransformsLRSbS);
        m_HRMaskFilter->SetTransformType(m_TransformType);
    }

    m_HRMaskFilter->SetHRImage(m_ImageHR);
    m_HRMaskFilter->SetInputLRImages(m_ImagesLR);
    m_HRMaskFilter->Update();

    m_ImageMaskHR = m_HRMaskFilter->GetOutput();




    this->InitializePSF();
    this->HComputation();






    double e = 0.1; //current change between two consecutive estimate.
    int i = 0;

    while( (e>0) && (i<m_Nloops) )
    {
      e = ComputeIterativeBackProjection(m_Nlm,m_Beta,m_MedianIBP);
      i++;
      std::cout<<"Loop "<<i+1<<", current e: "<<e<<std::endl;
    }

    //Copy current estimate to the output HR image
    itkDuplicator::Pointer duplicator = itkDuplicator::New();
    duplicator->SetInputImage( m_CurrentImageHR );
    duplicator->Update();
    m_OutputHRImage = duplicator->GetOutput();


}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionIBPFilter::HComputation()
{

      typedef itk::ContinuousIndex< double, 3 > itkContinuousIndex;
      typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;

      std::cout<<"Computing the matrix H (y=Hx) + fill y and x. \n";
      //Principle: for each voxel of the LR images, we compute the influence of each voxel of the PSF (centered at the current LR voxel) and add the corresponding influence value (PSF value * interpolation weight) in the matrix H

      // Set size of matrices
      unsigned int ncols = m_ImageHR->GetLargestPossibleRegion().GetNumberOfPixels();

      //m_offset is used to fill correctly the vector m_Y with the input LR image values (offset for the linear index)
      m_Offset.resize(m_ImagesLR.size());
      unsigned int nrows = 0;

      for(unsigned int im = 0; im < m_ImagesLR.size(); im++)
      {
        m_Offset[im] = nrows;
        nrows += m_ImagesLR[im]->GetLargestPossibleRegion().GetNumberOfPixels();
      }

      m_H.set_size(nrows, ncols);
      m_Y.set_size(nrows);
      m_Y.fill(0.0);
      m_X.set_size(ncols);
      m_X.fill(0.0);

      //linear index : an integer value corresponding to (x,y,z) triplet coordinates (ITK index)
      unsigned int lrLinearIndex = 0;
      unsigned int hrLinearIndex = 0;

      //Temporary variables
      itkImage::IndexType lrIndex;  //index of the current voxel in the LR image
      itkImage::PointType lrPoint;  //physical point location of lrIndex
      itkImage::IndexType psfIndex; //index of the current voxel in the PSF
      itkImage::PointType psfPoint; //physical point location of psfIndex
      itkImage::PointType transformedPoint; //Physical point location after applying affine transform
      itkContinuousIndex  hrContIndex;  //continuous index in HR image of psfPoint
      itkImage::IndexType hrIndex;  //index in HR image of interpolated psfPoint

      //We use linear interpolation for the estimation of point influence in matrix H
      typedef itk::BSplineInterpolationWeightFunction<double, 3, 1> itkBSplineFunction;

      itkBSplineFunction::Pointer bsplineFunction = itkBSplineFunction::New();
      itkBSplineFunction::WeightsType bsplineWeights;
      bsplineWeights.SetSize(8); // (bsplineOrder + 1)^3
      itkBSplineFunction::IndexType   bsplineStartIndex;
      itkBSplineFunction::IndexType   bsplineEndIndex;
      itkBSplineFunction::SizeType    bsplineSize;
      itkImage::RegionType            bsplineRegion;

      //Get the size of the HR image
      itkImage::SizeType  hrSize  = m_ImageHR->GetLargestPossibleRegion().GetSize();

      std::cout<<"loop over LR images\n";
      for(unsigned int i=0; i<m_ImagesLR.size(); i++)
      {
        std::cout<<"Adding image "<<i+1<<"\n";

        //Get the size of the current LR image
        itkImage::SizeType  lrSize  = m_ImagesLR[i]->GetLargestPossibleRegion().GetSize();

        //Instantiate an iterator over the current LR image
        itkIteratorWithIndex itLRImage(m_ImagesLR[i],m_ImagesLR[i]->GetLargestPossibleRegion());

        //Instantiate an iterator over the current PSF
        itkIteratorWithIndex itPSF(m_PSF[i],m_PSF[i]->GetLargestPossibleRegion());

        //Set the correct direction for the PSF of the current image
        m_PSF[i]->SetDirection(m_ImagesLR[i]->GetDirection());

        //Compute temporary variables for setting the correct origin of the PSF for every voxel of the LR image
        itkImage::SizeType psfSize = m_PSF[i]->GetLargestPossibleRegion().GetSize();
        itkContinuousIndex psfIndexCenter;
        psfIndexCenter[0] = (psfSize[0]-1)/2.0;
        psfIndexCenter[1] = (psfSize[1]-1)/2.0;
        psfIndexCenter[2] = (psfSize[2]-1)/2.0;
        itkImage::PointType psfPointCenter;
        m_PSF[i]->TransformContinuousIndexToPhysicalPoint(psfIndexCenter,psfPointCenter);
        itkImage::PointType psfOrigin;

        //Loop over the voxels of the current LR image
        for(itLRImage.GoToBegin(); !itLRImage.IsAtEnd(); ++itLRImage)
        {

          //Test on padding value (speed-up and keep H as sparse as possible)
          if(itLRImage.Get() > m_PaddingValue)
          {

            //Coordinate in the current LR image
            lrIndex = itLRImage.GetIndex();

            //std::cout<<"current LR index:"<<lrIndex[0]<<" "<<lrIndex[1]<<" "<<lrIndex[2]<<" -------------- \n";

            //Compute the corresponding linear index of lrIndex
            lrLinearIndex = m_Offset[i] + lrIndex[0] + lrIndex[1]*lrSize[0] + lrIndex[2]*lrSize[0]*lrSize[1];

            //Fill m_Y
            m_Y[lrLinearIndex] = itLRImage.Get();

            //Change the origin of the PSF so that itLRImage location corresponds to the center of the PSF
            m_ImagesLR[i]->TransformIndexToPhysicalPoint(lrIndex,lrPoint);
            psfOrigin[0] = lrPoint[0] - psfPointCenter[0];
            psfOrigin[1] = lrPoint[1] - psfPointCenter[1];
            psfOrigin[2] = lrPoint[2] - psfPointCenter[2];
            m_PSF[i]->SetOrigin(psfOrigin);

            //std::cout<<"current lr point:"<<lrPoint[0]<<" "<<lrPoint[1]<<" "<<lrPoint[2]<<" \n";
            //std::cout<<"psfSize:"<<psfSize[0]<<" "<<psfSize[1]<<" "<<psfSize[2]<<" \n";
            //std::cout<<"psfIndexCenter:"<<psfIndexCenter[0]<<" "<<psfIndexCenter[1]<<" "<<psfIndexCenter[2]<<"\n";
            //std::cout<<"psfPointCenter:"<<psfPointCenter[0]<<" "<<psfPointCenter[1]<<" "<<psfPointCenter[2]<<"\n";
            //std::cout<<"psf origin:"<<psfOrigin[0]<<" "<<psfOrigin[1]<<" "<<psfOrigin[2]<<" \n";

            //Loop over the m_PSF
            for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
            {

              if(itPSF.Get() > 0)
              {

                //Get coordinate in m_PSF
                psfIndex = itPSF.GetIndex();

                //Compute the physical point of psfIndex
                m_PSF[i]->TransformIndexToPhysicalPoint(psfIndex,psfPoint);

                //Apply estimated affine transform to psfPoint (need to apply the inverse since the transform goes from the HR image to the LR image)
                if(m_TransformType == AFFINE)
                {
                    transformedPoint = m_InverseTransformsLRAffine[i]->TransformPoint(psfPoint);
                }
                else
                {
                    transformedPoint = m_InverseTransformsLRSbS[i]->TransformPoint(psfPoint);
                }


                //Get back to the index in the HR image
                m_ImageHR->TransformPhysicalPointToContinuousIndex(transformedPoint,hrContIndex);
                //std::cout<<"hrContIndex no check: "<<hrContIndex[0]<<" "<<hrContIndex[1]<<" "<<hrContIndex[2]<<" \n ";

                //Check if the continuous index hrContIndex is inside the HR image
                if (m_ImageHR->GetLargestPossibleRegion().IsInside(hrContIndex))
                {

                  //std::cout<<"hrContIndex: "<<hrContIndex[0]<<" "<<hrContIndex[1]<<" "<<hrContIndex[2]<<" \n ";
                  //Get the interpolation weight using itkBSplineInterpolationWeightFunction
                  bsplineFunction->Evaluate(hrContIndex,bsplineWeights,bsplineStartIndex);
                  //std::cout<<"bsplineIndex: "<<bsplineStartIndex[0]<<" "<<bsplineStartIndex[1]<<" "<<bsplineStartIndex[2]<<"\n";

                  //Get the support size for interpolation
                  bsplineSize = bsplineFunction->GetSupportSize();

                  //Check if the bspline support region is inside the HR image
                  bsplineEndIndex[0] = bsplineStartIndex[0] + bsplineSize[0];
                  bsplineEndIndex[1] = bsplineStartIndex[1] + bsplineSize[1];
                  bsplineEndIndex[2] = bsplineStartIndex[2] + bsplineSize[2];

                  if( (m_ImageHR->GetLargestPossibleRegion().IsInside(bsplineStartIndex)) && (m_ImageHR->GetLargestPossibleRegion().IsInside(bsplineEndIndex)) )
                  {

                    //Set the support region
                    bsplineRegion.SetSize(bsplineSize);
                    bsplineRegion.SetIndex(bsplineStartIndex);

                    //Instantiate an iterator on HR image over the bspline region
                    itkIteratorWithIndex itHRImage(m_ImageHR,bsplineRegion);

                    //linear index of bspline weights
                    unsigned int weightLinearIndex = 0;

                    //Loop over the support region
                    for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage)
                    {

                      //Get coordinate in HR image
                      hrIndex = itHRImage.GetIndex();

                      //Compute the corresponding linear index
                      hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];

                      //Add weight*PSFValue to the corresponding element in H
                      m_H(lrLinearIndex,hrLinearIndex) += itPSF.Get() * bsplineWeights[weightLinearIndex];
                      weightLinearIndex += 1;
                      //std::cout<<"HR point:"<<hrIndex[0]<<" "<<hrIndex[1]<<" "<<hrIndex[2]<<", weight:"<<m_H(lrLinearIndex,hrLinearIndex)<<"\n";

                    } //end of loop over the support region

                  } //end of if support region inside HR image

                } //end of check of hrContIndex

              } //end of if itPSF > 0

            } //end of loop over PSF

          } //end of if on padding value

        } //end of loop over voxels of LR image

      } //end of loop over the set of LR images

      // Normalize m_H
      for(unsigned int i = 0; i < m_H.rows(); i++)
      {

        double sum = m_H.sum_row(i);

        vnl_sparse_matrix<float>::row & r = m_H.get_row(i);
        vnl_sparse_matrix<float>::row::iterator col_iter;

        for (col_iter = r.begin(); col_iter != r.end(); ++col_iter)
          (*col_iter).second = (*col_iter).second / sum;
      }

      //Fill m_X
      //Instantiate an iterator on HR image
      itkIteratorWithIndex itHRImage(m_ImageHR,m_ImageHR->GetLargestPossibleRegion());
      for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage){
        hrIndex = itHRImage.GetIndex();
        hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];
        m_X[hrLinearIndex] = itHRImage.Get();
      }
}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionIBPFilter::Initialize()
{
    m_SimulatedImagesLR.resize(m_ImagesLR.size());
    if(m_TransformType == AFFINE)
    {
        m_InverseTransformsLRAffine.resize(m_TransformsLRAffine.size());
    }
    else
    {
        m_InverseTransformsLRSbS.resize(m_TransformsLRSbS.size());
    }
   // m_InverseTransformsLR.resize(m_TransformsLR.size());

    typedef itk::ImageDuplicator< itkImage >   itkDuplicator;
    //duplicate the input image into the output images to keep all header information
    itkDuplicator::Pointer duplicator = itkDuplicator::New();
    duplicator->SetInputImage( m_ImageHR );
    duplicator->Update();
    m_OutputHRImage = duplicator->GetOutput();
    m_OutputHRImage->FillBuffer(0);

    itkDuplicator::Pointer duplicator2 = itkDuplicator::New();
    duplicator2->SetInputImage( m_ImageHR );
    duplicator2->Update();
    m_CurrentImageHR = duplicator2->GetOutput();

    itkDuplicator::Pointer duplicator3 = itkDuplicator::New();
    duplicator3->SetInputImage( m_ImageHR );
    duplicator3->Update();
    m_ImageMaskHR = duplicator3->GetOutput();

    for(unsigned int i=0;i<m_ImagesLR.size();i++)
    {

      //duplicate the input LR image into the output simulated image to keep all header information
      itkDuplicator::Pointer duplicator = itkDuplicator::New();
      duplicator->SetInputImage( m_ImagesLR[i] );
      duplicator->Update();
      m_SimulatedImagesLR[i] = duplicator->GetOutput();
      m_SimulatedImagesLR[i]->FillBuffer(0);
      //m_TransformsLR[i]->GetInverse(m_InverseTransformsLR[i]);

      if(m_TransformType == AFFINE)
      {
          m_InverseTransformsLRAffine[i] = itkAffineTransform::New();
          m_InverseTransformsLRAffine[i]->SetCenter( m_TransformsLRAffine[i]->GetCenter());
          m_TransformsLRAffine[i]->GetInverse(m_InverseTransformsLRAffine[i]);

      }
      else
      {
          m_InverseTransformsLRSbS[i] =  btkSliceBySliceTransform::New();
          m_InverseTransformsLRSbS[i]->SetFixedParameters( m_TransformsLRSbS[i]->GetFixedParameters());
          m_TransformsLRSbS[i]->GetInverse(m_InverseTransformsLRSbS[i]);
      }
    }


}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionIBPFilter::InitializePSF()
{
  //Principle: We build the PSF in LR image space (simple boxcar PSF = one anisotropic voxel) which is then interpolated or oversampled in SR space

  typedef itk::BSplineInterpolateImageFunction<itkImage, double, double>  itkBSplineInterpolator;
  typedef itk::ContinuousIndex< double, 3 > itkContinuousIndex;
  typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;

  std::cout<<"Initializing the PSF\n";

  //set the correct number of PSF (one PSF for one LR image -> this allows us to use images with different LR resolution)
  m_PSF.resize(m_ImagesLR.size());

  for(unsigned int i=0; i != m_PSF.size(); i++)
  {

    // 1- build the boxcar PSF in LR space (one anisotropic voxel)
    itkImage::Pointer LRPSF = itkImage::New();
    itkImage::IndexType lrIndex;
    itkImage::SizeType lrSize;
    //We enlarge by 1 voxel the LR PSF by null voxel for proper interpolated values.
    int border = 1; //to get proper interpolated value when close to the image boundary
    lrIndex[0] = 0;  lrIndex[1] = 0;  lrIndex[2] = 0;
    lrSize[0] = 1+2*border;   lrSize[1] = 1+2*border;   lrSize[2] = 1+2*border;
    itkImage::SpacingType lrSpacing = m_ImagesLR[i]->GetSpacing();

    //Allocate the LR PSF
    itkImage::RegionType lrRegion;
    lrRegion.SetSize(lrSize);
    lrRegion.SetIndex(lrIndex);
    LRPSF->SetRegions(lrRegion);
    LRPSF->SetSpacing(lrSpacing);
    LRPSF->Allocate();
    LRPSF->FillBuffer(0);
    std::cout<<"LRPSF size : "   <<lrSize[0]   <<" "<<lrSize[1]   <<" "<<lrSize[2]<<"\n";
    std::cout<<"LRPSF spacing : "<<lrSpacing[0]<<" "<<lrSpacing[1]<<" "<<lrSpacing[2]<<"\n";

    itkImage::IndexType lrIndexCenter;
    lrIndexCenter[0] = (lrSize[0]-1)/2.0;
    lrIndexCenter[1] = (lrSize[1]-1)/2.0;
    lrIndexCenter[2] = (lrSize[2]-1)/2.0;
    itkImage::PointType lrPointCenter;
    LRPSF->TransformIndexToPhysicalPoint(lrIndexCenter,lrPointCenter);
    std::cout<<"LRPSF index center : "<<lrIndexCenter[0]<<" "<<lrIndexCenter[1]<<" "<<lrIndexCenter[2]<<"\n";
    std::cout<<"LRPSF physical center : "<<lrPointCenter[0]<<" "<<lrPointCenter[1]<<" "<<lrPointCenter[2]<<"\n";
    LRPSF->SetPixel(lrIndexCenter,1.0);


    // 2- Allocate the corresponding HR PSF
    itkImage::IndexType hrIndex;
    itkImage::SizeType hrSize;
    hrSize[0] = (int)ceil(lrSpacing[0] / m_ImageHR->GetSpacing()[0]) + 2;
    hrSize[1] = (int)ceil(lrSpacing[1] / m_ImageHR->GetSpacing()[1]) + 2;
    hrSize[2] = (int)ceil(lrSpacing[2] / m_ImageHR->GetSpacing()[2]) + 2;
    hrIndex[0] = 0;  hrIndex[1] = 0;  hrIndex[2] = 0;

    std::cout<<"HR PSF spacing : "<<m_ImageHR->GetSpacing()[0]<<" "<<m_ImageHR->GetSpacing()[1]<<" "<<m_ImageHR->GetSpacing()[2]<<std::endl;
    std::cout<<"HR PSF size : "<<hrSize[0]<<" "<<hrSize[1]<<" "<<hrSize[2]<<std::endl;

    itkImage::RegionType hrRegion;
    hrRegion.SetSize(hrSize);
    hrRegion.SetIndex(hrIndex);
    m_PSF[i] = itkImage::New();
    m_PSF[i]->SetRegions(hrRegion);
    m_PSF[i]->SetSpacing(m_ImageHR->GetSpacing());
    m_PSF[i]->Allocate();
    m_PSF[i]->FillBuffer(0.0);
    itkContinuousIndex hrIndexCenter;
    hrIndexCenter[0] = (hrSize[0]-1)/2.0;
    hrIndexCenter[1] = (hrSize[1]-1)/2.0;
    hrIndexCenter[2] = (hrSize[2]-1)/2.0;
    itkImage::PointType hrPointCenter;
    m_PSF[i]->TransformContinuousIndexToPhysicalPoint(hrIndexCenter,hrPointCenter);

    std::cout<<"modification of the origin of m_PSF so that LRPSF and m_PSF are centered.\n";
    itkImage::PointType hrOrigin;
    hrOrigin[0] = lrPointCenter[0] - hrPointCenter[0];
    hrOrigin[1] = lrPointCenter[1] - hrPointCenter[1];
    hrOrigin[2] = lrPointCenter[2] - hrPointCenter[2];
    m_PSF[i]->SetOrigin(hrOrigin);
    m_PSF[i]->TransformContinuousIndexToPhysicalPoint(hrIndexCenter,hrPointCenter);
    std::cout<<"HRPSF origin : "<<hrOrigin[0]<<" "<<hrOrigin[1]<<" "<<hrOrigin[2]<<std::endl;
    std::cout<<"HRPSF physical center : "<<hrPointCenter[0]<<" "<<hrPointCenter[1]<<" "<<hrPointCenter[2]<<std::endl;

    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(m_InterpolationOrderPSF);
    bsInterpolator->SetInputImage(LRPSF);

    itkIteratorWithIndex itPSF(m_PSF[i],m_PSF[i]->GetLargestPossibleRegion());
    double nbSamples = 20; //parameter for oversampled HR PSF.
    std::vector<float> sigma(3); //parameter for 3D Gaussian PSF

    switch (m_PsfType)
    {
      case 0:
        std::cout<<"3D interpolated boxcar using "<<m_InterpolationOrderPSF<<" order B-Spline."<<std::endl;
        //Loop over voxels of HR PSF
        for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
        {

          //Coordinate in HR image
          hrIndex = itPSF.GetIndex();

          //Coordinate in physical space
          itkImage::PointType hrPoint;
          m_PSF[i]->TransformIndexToPhysicalPoint(hrIndex,hrPoint);

          //Continuous coordinate in LR image
          itkContinuousIndex lrContIndex;
          LRPSF->TransformPhysicalPointToContinuousIndex(hrPoint,lrContIndex);

          //Set interpolated value to m_PSF
          itPSF.Set(bsInterpolator->EvaluateAtContinuousIndex(lrContIndex));
        }
        break;

      case 1:
        std::cout<<"3D oversampled boxcar using nearest neighbour interpolation.\n";
        //set nearest neighbour interpolation mode
        bsInterpolator->SetSplineOrder(0);

        //Loop over voxels of HR PSF
        for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
        {

          double sum = 0;
          int x,y,z;
          //Coordinate in HR image
          hrIndex = itPSF.GetIndex();
          itkContinuousIndex hrContIndex;
          itkImage::PointType hrPoint;
          itkContinuousIndex lrContIndex;

          //Coordinate in physical space
          m_PSF[i]->TransformIndexToPhysicalPoint(hrIndex,hrPoint);
          //Continuous coordinate in LR image
          LRPSF->TransformPhysicalPointToContinuousIndex(hrPoint,lrContIndex);

          for(z=0; z<nbSamples; z++)
            for(y=0; y<nbSamples; y++)
              for(x=0; x<nbSamples; x++)
              {

                hrContIndex[0] = hrIndex[0] - 0.5 + 1.0*x/nbSamples;
                hrContIndex[1] = hrIndex[1] - 0.5 + 1.0*y/nbSamples;
                hrContIndex[2] = hrIndex[2] - 0.5 + 1.0*z/nbSamples;

                //Coordinate in physical space
                m_PSF[i]->TransformContinuousIndexToPhysicalPoint(hrContIndex,hrPoint);

                //Continuous coordinate in LR image
                LRPSF->TransformPhysicalPointToContinuousIndex(hrPoint,lrContIndex);

                sum += bsInterpolator->EvaluateAtContinuousIndex(lrContIndex);
              }

          //Set oversampled value to m_PSF
          itPSF.Set(sum);
        }
        break;

        case 2:
          std::cout<<"Compute a 3D Gaussian PSF.\n";
          //Note that we should also do oversampling to obtain a more accurate estimate of the PSF
          //Here, we use a piecewise constant approximation
          //Set the FWHM equal to the voxel size
          //FWHM = 2.3548 sigma
          sigma[0] = lrSpacing[0] / 2.3548;
          sigma[1] = lrSpacing[1] / 2.3548;
          sigma[2] = lrSpacing[2] / 2.3548;
          std::cout<<"Sigma :"<<sigma[0]<<" "<<sigma[1]<<" "<<sigma[2]<<"\n";
          for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
          {
            hrIndex = itPSF.GetIndex();
            float x = hrIndex[0] - hrIndexCenter[0];
            float y = hrIndex[1] - hrIndexCenter[1];
            float z = hrIndex[2] - hrIndexCenter[2];
            float value = (x*x)/(2*sigma[0]*sigma[0]) + (y*y)/(2*sigma[1]*sigma[1]) + (z*z)/(2*sigma[2]*sigma[2]);
            itPSF.Set(exp(-value));
          }
          break;

      default:
        std::cout<<"Invalid choice for the psf building"<<std::endl;
        break;
    }

    //for memory saving, we limit the number of non-null voxels
    for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
    {
      if(itPSF.Get() < 0.01)
        itPSF.Set(0);
    }

    //Normalization of the PSF
    double sum = 0;
    for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF) sum += itPSF.Get();
    if(sum>0)
      for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
        itPSF.Set( itPSF.Get() / sum );

    std::cout<<"Positive values of the PSF for the "<<i+1<<"th LR image : \n";
    for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
      if(itPSF.Get()>0)
      {
        hrIndex = itPSF.GetIndex();
        std::cout<<itPSF.Get()<<" ("<<hrIndex[0]<<", "<<hrIndex[1]<<", "<<hrIndex[2]<<") "<<std::endl;
      }

    std::cout<<"Set the HR origin to 0 since the PSF now should be centered with image voxels"<<std::endl<<std::endl;
    hrOrigin[0] = 0;
    hrOrigin[1] = 0;
    hrOrigin[2] = 0;
    m_PSF[i]->SetOrigin(hrOrigin);
  }
}
//-----------------------------------------------------------------------------------------------------------
void HighResolutionIBPFilter::UpdateX()
{
  typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;
  std::cout<<"Update x"<<std::endl;
  m_X.fill(0.0);
  itkImage::IndexType hrIndex;
  unsigned int hrLinearIndex = 0;

  itkImage::SizeType  hrSize  = m_CurrentImageHR->GetLargestPossibleRegion().GetSize();

  itkIteratorWithIndex itHRImage(m_CurrentImageHR,m_CurrentImageHR->GetLargestPossibleRegion());
  for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage)
  {
    hrIndex = itHRImage.GetIndex();
    hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];
    m_X[hrLinearIndex] = itHRImage.Get();
  }
}
//-----------------------------------------------------------------------------------------------------------
double HighResolutionIBPFilter::ComputeIterativeBackProjection( int & nlm, float & beta, int & medianIBP)
{
    typedef itk::BSplineInterpolateImageFunction< itkImage, double, double >          itkBSplineInterpolator;
    typedef itk::SubtractImageFilter < itkImage, itkImage >                          itkSubtractImageFilter;
    typedef itk::ResampleImageFilter< itkImage, itkImage >                            itkResampleFilter;
    typedef itk::ImageRegionIteratorWithIndex< itkImage >                           itkIteratorWithIndex;
    typedef itk::AddImageFilter < itkImage, itkImage >                               itkAddImageFilter;
    typedef itk::AbsoluteValueDifferenceImageFilter < itkImage, itkImage, itkImage >  itkAbsoluteValueDifferenceImageFilter;
    typedef itk::ImageRegionIterator< itkImage > itkIterator;
    typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
    typedef itk::StatisticsImageFilter< itkImage >                            itkStatisticsImageFilter;

    std::cout<<"Do iterated back projection"<<std::endl;
    std::cout<<"NLM Filtering type: "<<nlm<<std::endl;

    std::cout<<"Compute y-Hx"<<std::endl;
    this->UpdateX();
    m_SimuLRImagesFilter->SetH(m_H);
    m_SimuLRImagesFilter->SetLRImages(m_ImagesLR);
    m_SimuLRImagesFilter->SetX(m_X);
    m_SimuLRImagesFilter->SetOffset(m_Offset);
    m_SimuLRImagesFilter->Update();

    m_SimulatedImagesLR = m_SimuLRImagesFilter->GetOutput();

    std::cout<<"Update the current HR image"<<std::endl;
    //Initialize to 0 the output HR image
    m_OutputHRImage->FillBuffer(0);


    //parameters for interpolation (bspline interpolator)
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(m_InterpolationOrderIBP);

    std::string s = ""; // ? UseLess ?
    std::vector< itkImage::Pointer >   errorImages;            //difference image between simulated images and observed images.
    errorImages.resize(m_ImagesLR.size());

    for(unsigned int i=0; i< m_ImagesLR.size(); i++)
    {

      //compute the difference between LR input and the simulated LR image
      itkSubtractImageFilter::Pointer subtractFilter = itkSubtractImageFilter::New ();
      subtractFilter->SetInput1(m_ImagesLR[i]);
      subtractFilter->SetInput2(m_SimulatedImagesLR[i]); //TODO: Before perform this be sure that you compute the simulateLRImage
      subtractFilter->Update();

      std::ostringstream oss ;
      oss << i+1 ;
      //s = "ibp_"+oss.str()+"_substract.nii.gz";
      //data.WriteOneImage(subtractFilter->GetOutput(), s);

      //interpolate the LR difference
      itkResampleFilter::Pointer resample = itkResampleFilter::New();
      if(m_TransformType == AFFINE)
      {
          resample->SetTransform(m_TransformsLRAffine[i]);
      }
      else
      {
          resample->SetTransform(m_TransformsLRSbS[i]);
      }

      resample->SetInterpolator(bsInterpolator);
      resample->UseReferenceImageOn();
      resample->SetReferenceImage(m_CurrentImageHR);
      resample->SetInput(subtractFilter->GetOutput());
      resample->Update();

      errorImages[i] = resample->GetOutput();

      //s = "ibp_"+oss.str()+"_resample.nii.gz";
      //data.WriteOneImage(errorImages[i], s);



      //Add the interpolated differences
      itkAddImageFilter::Pointer addFilter = itkAddImageFilter::New ();
      addFilter->SetInput1(m_OutputHRImage);
      addFilter->SetInput2(resample->GetOutput());
      addFilter->Update();

     m_OutputHRImage = addFilter->GetOutput();
   /*    std::string s1 = "switch_out.nii.gz";
        data.WriteOneImage(data.m_outputHRImage, s1);
        */

      //s = "ibp_"+oss.str()+"_simulated.nii.gz";
      //data.WriteOneImage(data.m_simulatedInputLRImages[i], s);

    }

    //Iterator over the output image filled with backprojected error.
    itkIteratorWithIndex itImage(m_OutputHRImage,m_OutputHRImage->GetLargestPossibleRegion());

    /*
        //Normalize the resampled difference image
        for(itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage)
          itImage.Set( itImage.Get() / data.m_inputLRImages.size() );
    */

            //std::string s2 = "switch_in.nii.gz";

//FIXME : switch for filtering choice runs under MacOSX but not under Debian !!!
/*
    switch(medianIBP)
    {
      case 0:
      std::cout<<"Averaging the error maps (obtained for each LR image)\n";
        m_ImageHR->FillBuffer(0);

        for(unsigned int i=0; i<m_ImagesLR.size(); i++)
        {
          //Add the interpolated differences
          itkAddImageFilter::Pointer addFilter = itkAddImageFilter::New ();
          addFilter->SetInput1(m_OutputHRImage);
          addFilter->SetInput2(errorImages[i]);
          addFilter->Update();
          m_ImageHR = addFilter->GetOutput();
        }

        //data.WriteOneImage(m_OutputHRImage, s2);

        break;
      case 1:
          std::cout<<"Compute the median of the error map at each voxel\n";
          for(itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage)
          {
            itkImage::IndexType p = itImage.GetIndex();
            std::vector<float>  v;
            for(unsigned int i=0; i<m_ImagesLR.size(); i++)
              v.push_back(errorImages[i]->GetPixel(p));
            std::sort(v.begin(), v.end());

            itImage.Set(v.size() * v[(int)(v.size()/2)]);

          }

        break;
      default:
        std::cout<<"Invalid choice for the median IBP parameter (0 or 1).\n";
        break;
    }
*/
            //Normalize the resampled difference image
            for(itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage)
            {
              double value = itImage.Get() / m_ImagesLR.size();
              itImage.Set( value );
            }


    if(nlm==1)
    {
      std::cout<<"Smooth the error map using the current reconstructed image as reference for NLM filter --------------------------"<<std::endl;
      btkNLMTool<float> myTool;
      myTool.SetInput(m_OutputHRImage);
      myTool.SetMaskImage(m_ImageMaskHR);
      myTool.SetDefaultParameters();
      myTool.SetReferenceImage(m_CurrentImageHR);
      myTool.SetSmoothing(beta);
      myTool.SetBlockwiseStrategy(1); //0 pointwise, 1 block, 2 fast block


      myTool.ComputeOutput();
      myTool.GetOutput(m_OutputHRImage);
      //s = "ibp_nlm_error.nii.gz";
      //data.WriteOneImage(data.m_outputHRImage, s);
    }

    //Update the HR image correspondly
    itkAddImageFilter::Pointer addFilter2 = itkAddImageFilter::New ();
    addFilter2->SetInput1(m_OutputHRImage);
    addFilter2->SetInput2(m_CurrentImageHR);
    addFilter2->Update();

    m_OutputHRImage = addFilter2->GetOutput();
    //s = "ibp_updated.nii.gz";
    //data.WriteOneImage(data.m_currentHRImage, s);

    if(nlm==2)
    {
      std::cout<<"Smooth the current reconstructed image ------------------ \n";
      btkNLMTool<float> myTool;
      myTool.SetInput(m_OutputHRImage);
      myTool.SetMaskImage(m_ImageMaskHR);
      myTool.SetDefaultParameters();
      myTool.SetSmoothing(beta);
      myTool.SetBlockwiseStrategy(1); //0 pointwise, 1 block, 2 fast block
      myTool.ComputeOutput();
      myTool.GetOutput(m_OutputHRImage);
      //s = "ibp_nlm_smooth.nii.gz";
      //data.WriteOneImage(data.m_currentHRImage, s);
    }

    std::cout<<"Compute the changes between the two consecutive estimates\n";
    itkAbsoluteValueDifferenceImageFilter::Pointer absoluteValueDifferenceFilter = itkAbsoluteValueDifferenceImageFilter::New ();
    absoluteValueDifferenceFilter->SetInput1(m_OutputHRImage);
    absoluteValueDifferenceFilter->SetInput2(m_CurrentImageHR);
    absoluteValueDifferenceFilter->Update();

    double magnitude       = 0.0;
    double numberOfPoints  = 0.0;
    itkImage::Pointer  tmpImage = absoluteValueDifferenceFilter->GetOutput();
    itkIterator itTmpImage(tmpImage, tmpImage->GetLargestPossibleRegion());
    itkIterator itMaskHRImage(m_ImageMaskHR,m_ImageMaskHR->GetLargestPossibleRegion());

    for(itMaskHRImage.GoToBegin(),itTmpImage.GoToBegin(); !itMaskHRImage.IsAtEnd(); ++itMaskHRImage, ++itTmpImage)
    {
      if(itMaskHRImage.Get() > 0)
      {
        magnitude += itTmpImage.Get();
        numberOfPoints += 1;
      }
    }

    double meanMagnitude = magnitude/numberOfPoints;
    std::cout<<"Current mean change: "<<meanMagnitude<<std::endl;

    std::cout<<"Copying new estimate to current estimate HR image"<<std::endl;
    //Not efficient strategy but clearer to understand the code
    itkDuplicator::Pointer duplicator = itkDuplicator::New();
    duplicator->SetInputImage( m_OutputHRImage );
    duplicator->Update();
    m_CurrentImageHR = duplicator->GetOutput();

    itkStatisticsImageFilter::Pointer statisticsImageFilter = itkStatisticsImageFilter::New ();
    statisticsImageFilter->SetInput(m_CurrentImageHR);
    statisticsImageFilter->Update();
    std::cout << "Stat of the current HR estimate: "<<std::endl<<"Mean: " << statisticsImageFilter->GetMean() << std::endl;
    std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
    std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
    std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;

    double manjonCriterion = 0.002*statisticsImageFilter->GetSigma(); //As defined in Manjon et al. 2010
    std::cout<<"stopping criterion : "<<manjonCriterion<<std::endl;

    if(meanMagnitude < manjonCriterion)
      meanMagnitude = 0;

    return meanMagnitude;
}


}


