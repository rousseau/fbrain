/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 27/05/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.
  
==========================================================================*/

#include "btkSRHMatrixComputation.hxx"
#include "btkSliceBySliceTransform.h"
#include "itkEuler3DTransform.h"


namespace btk
{
//-------------------------------------------------------------------------------------------------

template < class TImage >
SRHMatrixComputation< TImage >::SRHMatrixComputation():m_PSFType(0)
{

}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void SRHMatrixComputation< TImage >::Initialize()
{
    m_IsHComputed = false;
    m_NumberOfLRImages = m_Images.size();

    m_Regions.resize(m_NumberOfLRImages);
    m_PSFs.resize(m_NumberOfLRImages);

    std::cout<<"PSF Type : "<<m_PSFType<<std::endl;
    for(unsigned int i = 0; i< m_NumberOfLRImages; i++)
    {
        //m_Regions[i] = m_Images[i]->GetLargestPossibleRegion();
        m_Regions[i] = m_Masks[i]->GetAxisAlignedBoundingBoxRegion();

        if(m_PSFType == 0)
        {
            m_PSFs[i] = btk::GaussianPSF::New();
            m_PSF = btk::GaussianPSF::New();
        }
        else if(m_PSFType == 1)
        {
            m_PSFs[i] = btk::BoxCarPSF::New();
            m_PSF = btk::BoxCarPSF::New();
        }
        else if(m_PSFType == 2)
        {
            m_PSFs[i] = btk::GaussianPSF::New();
            m_PSF = btk::GaussianPSF::New();
        }
        else if(m_PSFType == 3)
        {
            m_PSFs[i] = btk::SincPSF::New();
            m_PSF = btk::SincPSF::New();
        }
        else if(m_PSFType == 4)
        {
            m_PSFs[i] = btk::HybridPSF::New();
            m_PSF = btk::HybridPSF::New();
        }
        else
        {
            m_PSFs[i] = btk::GaussianPSF::New();
            m_PSF = btk::GaussianPSF::New();
        }

    }
}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void SRHMatrixComputation< TImage >::Update()
{


    this->Initialize();

    typename PsfImageType::Pointer PSF = PsfImageType::New();

    //We use linear interpolation for the estimation of point influence in matrix H
    typedef itk::BSplineInterpolationWeightFunction<double, 3, 1> itkBSplineFunction;
    itkBSplineFunction::Pointer bsplineFunction = itkBSplineFunction::New();
    itkBSplineFunction::WeightsType bsplineWeights;
    bsplineWeights.SetSize(8); // (bsplineOrder + 1)^3
    itkBSplineFunction::IndexType   bsplineStartIndex;
    itkBSplineFunction::IndexType   bsplineEndIndex;
    itkBSplineFunction::SizeType    bsplineSize;
    RegionType                      bsplineRegion;

    m_OutputImageRegion = m_ReferenceImage -> GetLargestPossibleRegion();
    IndexType start_hr  = m_OutputImageRegion.GetIndex();
    SizeType  size_hr   = m_OutputImageRegion.GetSize();

    //m_XSize : size of the SR image (used in other functions)
    m_XSize.width  = size_hr[0];
    m_XSize.height = size_hr[1];
    m_XSize.depth  = size_hr[2];

    IndexType end_hr;
    end_hr[0] = start_hr[0] + size_hr[0] - 1 ;
    end_hr[1] = start_hr[1] + size_hr[1] - 1 ;
    end_hr[2] = start_hr[2] + size_hr[2] - 1 ;

    // Differential continuous indexes to perform the neighborhood iteration
    SpacingType spacing_lr = m_Images[0] -> GetSpacing();
    SpacingType spacing_hr = m_ReferenceImage -> GetSpacing();

    //spacing_lr[2] is assumed to be the lowest resolution
    //compute the index of the PSF in the LR image resolution
    std::vector<ContinuousIndexType> deltaIndexes;
    int npoints =  spacing_lr[2] / (2.0 * spacing_hr[2]) ;
    ContinuousIndexType delta;
    delta[0] = 0.0; delta[1] = 0.0;


    for (int i = -npoints ; i <= npoints; i++ )
    {
        delta[2] = i * 0.5 / (double) npoints;
        deltaIndexes.push_back(delta);
    }

    // Set size of matrices
    unsigned int ncols = m_OutputImageRegion.GetNumberOfPixels();

    //m_X.set_size( ncols ); //FIXME: In this filter we don't need X !



    unsigned int nrows = 0;
    for(unsigned int im = 0; im < m_NumberOfLRImages; im++)
    {
        //nrows += m_Regions[im].GetNumberOfPixels();
        nrows += m_Images[im]->GetLargestPossibleRegion().GetNumberOfPixels();
    }

    m_H->set_size(nrows, ncols);

    m_Y->set_size(nrows);
    m_Y->fill(0.0);

    std::cout<<"size X : "<<ncols<<std::endl;
    std::cout<<"size Y : "<<nrows<<std::endl;


    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage(m_ReferenceImage);

    unsigned int lrLinearIndex;
    unsigned int hrLinearIndex;
    unsigned int offset = 0;
    double psfValue;



    unsigned int im;
    //#pragma omp parallel for private(im) schedule(dynamic)

    for(im = 0; im < m_Images.size(); im++)
    {

        std::cout<<"Processing image "<<im+1<<std::endl;

        SizeType lrSize = m_Images[im]->GetLargestPossibleRegion().GetSize();
        SpacingType lrSpacing = m_Images[im]->GetSpacing();
        ConstIteratorType lrIt( m_Images[im], m_Images[im]->GetLargestPossibleRegion() );

        // for all voxels
        for(lrIt.GoToBegin(); !lrIt.IsAtEnd(); ++lrIt)
        {
            IndexType lrIndex = lrIt.GetIndex();
            PointType lrPoint;
            m_Images[im]->TransformIndexToPhysicalPoint(lrIndex, lrPoint);
            PointType srPoint = m_Transforms[im]->TransformPoint(lrPoint);

            // if point is not in the sr image we skip it
            if(!interpolator->IsInsideBuffer(srPoint))
            {
                continue;
            }

            // if point is not in the mask we skip it
            if(!m_Masks[im]->GetImage()->GetPixel(lrIndex) > 0)
            {
                continue;
            }

            //compute the linear index corresponding to the index of lr image
            lrLinearIndex = offset + lrIndex[0] + lrIndex[1]*lrSize[0]
                            + lrIndex[2]*lrSize[0]*lrSize[1];

            //Fill Y
            m_Y->operator()(lrLinearIndex) = lrIt.Get();

            //Initialization of the PSF
            m_PSF->SetDirection(m_Images[im]->GetDirection());
            m_PSF->SetLrSpacing(lrSpacing);
            m_PSF->SetSpacing(m_ReferenceImage->GetSpacing());

            //TODO: Like Stan, but I don't understand why !
            SizeType psfSize;
            psfSize[0] = (int)ceil(m_Images[im]->GetSpacing()[0] / m_ReferenceImage->GetSpacing()[0]) +2;
            psfSize[1] = (int)ceil(m_Images[im]->GetSpacing()[1] / m_ReferenceImage->GetSpacing()[1]) +2;
            psfSize[2] = (int)ceil(m_Images[im]->GetSpacing()[2] / m_ReferenceImage->GetSpacing()[2]) +2;

            m_PSF->SetSize(psfSize);
            //We center the PSF on the physical LR point
           // std::cout<<lrPoint<<std::endl;

            m_PSF->ConstructImage();
            m_PSF->SetCenter(lrPoint);
            // Get PSF support (image)
            PSF = m_PSF->GetPsfImage();
            // Save the result (do that with a break point in debug mode,
            // otherwise PSF is saved at each voxel).
            //btk::ImageHelper< itk::Image< float, 3 > >::WriteImage(PSF,"PSFImage.nii.gz");

            ImageRegionConstIteratorWithIndex< itk::Image< float, 3 > > itPsf(PSF, PSF->GetLargestPossibleRegion());
            // Loop over PSF voxels
            for(itPsf.GoToBegin(); !itPsf.IsAtEnd(); ++itPsf)
            {
                IndexType psfIndex = itPsf.GetIndex();
                PointType psfInLrSpacePoint, transformedPoint;
                //Physical point of psf Index, in Lr Space
                m_Images[im]->TransformIndexToPhysicalPoint(psfIndex,psfInLrSpacePoint);

                psfValue = itPsf.Get();


                //std::cout<<"At Psf index : "<<psfIndex<<" value : "<<psfValue<<std::endl;

                if(!psfValue > 0)
                {
                    continue;
                }

                // psfPoint in sr space
                //m_Transforms[im]->GetInverse(inverse);
                //transformedPoint = m_Transforms[im]->TransformPoint(psfInLrSpacePoint);
                transformedPoint = m_InverseTransforms[im]->TransformPoint(psfInLrSpacePoint);

                // if point is not in the sr image we skip it
                if(interpolator->IsInsideBuffer(transformedPoint))
                {
                    ContinuousIndexType srContIndex;
//                    // continuous index in sr image
                    m_ReferenceImage->TransformPhysicalPointToContinuousIndex(transformedPoint, srContIndex);


                    // STAN VERSION :


//                    bool isInsideHR = true;

//                    // FIXME This checking should be done for all points first, and discard the point
//                    // if al least one point is out of the reference image

//                    if ( (srContIndex[0] < start_hr[0]) || (srContIndex[0] > end_hr[0]) ||
//                         (srContIndex[1] < start_hr[1]) || (srContIndex[1] > end_hr[1]) ||
//                         (srContIndex[2] < start_hr[2]) || (srContIndex[2] > end_hr[2]) )
//                       isInsideHR = false;

//                    if ( isInsideHR )
//                    {
//                      double hrValue = interpolator -> Evaluate( transformedPoint );

//                      for(unsigned int n=0; n<interpolator -> GetContributingNeighbors(); n++)
//                      {
//                        IndexType hrIndex = interpolator -> GetIndex(n);
//                        IndexType hrDiffIndex;

//                        hrDiffIndex[0] = hrIndex[0] - start_hr[0];
//                        hrDiffIndex[1] = hrIndex[1] - start_hr[1];
//                        hrDiffIndex[2] = hrIndex[2] - start_hr[2];

//                        hrLinearIndex = hrDiffIndex[0] + hrDiffIndex[1]*size_hr[0] + hrDiffIndex[2]*size_hr[0]*size_hr[1];
//                        //std::cout<<"old value in H("<<lrLinearIndex<<","<<hrLinearIndex<<") : "<<m_H->operator()(lrLinearIndex, hrLinearIndex)<<std::endl;

//                        m_H->operator()(lrLinearIndex, hrLinearIndex) += psfValue/* * interpolator->GetOverlap(n)*/;

//                        //std::cout<<"new value in H("<<lrLinearIndex<<","<<hrLinearIndex<<") : "<<m_H->operator()(lrLinearIndex, hrLinearIndex)<<std::endl;


//                      }
//                    }


                    // FRANCOIS VERSION :

                    //std::cout<<"hrContIndex: "<<srContIndex[0]<<" "<<srContIndex[1]<<" "<<srContIndex[2]<<std::endl;
                    //Get the interpolation weight using itkBSplineInterpolationWeightFunction
                    bsplineFunction->Evaluate(srContIndex,bsplineWeights,bsplineStartIndex);
                    //std::cout<<"bsplineIndex: "<<bsplineStartIndex[0]<<" "<<bsplineStartIndex[1]<<" "<<bsplineStartIndex[2]<<std::endl;

                    //Get the support size for interpolation
                    bsplineSize = bsplineFunction->GetSupportSize();

                    //Check if the bspline support region is inside the HR image
                    bsplineEndIndex[0] = bsplineStartIndex[0] + bsplineSize[0];
                    bsplineEndIndex[1] = bsplineStartIndex[1] + bsplineSize[1];
                    bsplineEndIndex[2] = bsplineStartIndex[2] + bsplineSize[2];

                    if(m_ReferenceImage->GetLargestPossibleRegion().IsInside(bsplineStartIndex)
                            && m_ReferenceImage->GetLargestPossibleRegion().IsInside(bsplineEndIndex))
                    {
                        //Set the support region
                        bsplineRegion.SetSize(bsplineSize);
                        bsplineRegion.SetIndex(bsplineStartIndex);

                        //Instantiate an iterator on HR image over the bspline region
                        ImageRegionConstIteratorWithIndex< ImageType > itHRImage(m_ReferenceImage,bsplineRegion);

                        //linear index of bspline weights
                        unsigned int weightLinearIndex = 0;

                        //Loop over the support region
                        for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage)
                        {

                            //Get coordinate in HR image
                            IndexType hrIndex = itHRImage.GetIndex();

                            //Compute the corresponding linear index
                            hrLinearIndex = hrIndex[0] + hrIndex[1]*size_hr[0] + hrIndex[2]*size_hr[0]*size_hr[1];

                            //Add weight*PSFValue to the corresponding element in H
                            //
                            //std::cout<<"old value in H("<<lrLinearIndex<<","<<hrLinearIndex<<") : "<<m_H->operator()(lrLinearIndex, hrLinearIndex)<<std::endl;
                            if(std::isnan(psfValue))
                            {
                                std::cout<<"NAN PSF"<<std::endl;
                                btkCoutVariable(im);
                                btkCoutVariable(lrLinearIndex);
                                btkCoutVariable(hrLinearIndex);
                                btkCoutVariable(weightLinearIndex);
                                btkCoutVariable(lrIndex);
                            }

                            if(std::isnan(bsplineWeights[weightLinearIndex]))
                            {
                                std::cout<<"BSPLINE IS NAN"<<std::endl;
                                btkCoutVariable(im);
                                btkCoutVariable(lrLinearIndex);
                                btkCoutVariable(hrLinearIndex);
                                btkCoutVariable(weightLinearIndex);
                                btkCoutVariable(lrIndex);
                            }

                            m_H->operator()(lrLinearIndex, hrLinearIndex) += psfValue * bsplineWeights[weightLinearIndex];

                            //std::cout<<"new value in H("<<lrLinearIndex<<","<<hrLinearIndex<<") : "<<m_H->operator()(lrLinearIndex, hrLinearIndex)<<std::endl;
                            weightLinearIndex++;

                        } //end of loop over the support region

                    }// end if bspline index inside sr image

                } // if psf point is inside sr image



            }// for PSF voxels


        }//for voxels

        offset += m_Images[im]->GetLargestPossibleRegion().GetNumberOfPixels();
        //std::cout<<"offset "<<offset<<"/"<<UINT_MAX<<std::endl;
        assert(offset < UINT_MAX);
    }//for im






    for (unsigned int i = 0; i < m_H->rows(); i++)
    {
        double sum = 0.0;
        sum = m_H->sum_row(i);

        vnl_sparse_matrix< PrecisionType >::row & r = m_H->get_row(i);
        vnl_sparse_matrix< PrecisionType >::row::iterator col_iter = r.begin();
        for ( ;col_iter != r.end(); ++col_iter)
        {

            (*col_iter).second = (*col_iter).second / sum;

            //std::cout<<(*col_iter).second

        }


    }

//    for(unsigned int i = 0; i< m_H->rows(); i++)
//    {
//        for(unsigned j = 0; j< m_H->cols(); j++)
//        {
//            std::cout<<m_H->operator()(i,j)<<" ";
//        }

//        std::cout<<std::endl;
//    }


    m_IsHComputed = true;


    std::cout<<"H computed !"<<std::endl;
    // DEBUG :
    std::cout<<"Testing Y, X and simulated Y..."<<std::endl;
    this->TestFillingOfY();
    this->TestFillingOfX();
    this->SimulateY();
}

//-------------------------------------------------------------------------------------------------
template< class TImage >
void
SRHMatrixComputation< TImage >::TestFillingOfY()
{
    if(!m_IsHComputed )
    {
        this->Update();
    }

    // H should previoulsy be computed
    m_ImagesFilledWithY.resize(m_NumberOfLRImages);
    unsigned int offset = 0;

    for(unsigned int i = 0; i< m_NumberOfLRImages; i++)
    {
        SizeType lrSize = m_Images[i]->GetLargestPossibleRegion().GetSize();
        m_ImagesFilledWithY[i] = ImageType::New();
        m_ImagesFilledWithY[i] = btk::ImageHelper< ImageType >::CreateNewImageFromPhysicalSpaceOf(m_Images[i]);
        m_ImagesFilledWithY[i]->FillBuffer(0.0);

        IteratorType it(m_ImagesFilledWithY[i], m_ImagesFilledWithY[i]->GetLargestPossibleRegion());
        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            IndexType lrIndex = it.GetIndex();

            unsigned int lrLinearIndex = offset + lrIndex[0] +lrIndex[1]*lrSize[0]+ lrIndex[2]*lrSize[0]*lrSize[1];

            it.Set(m_Y->operator()(lrLinearIndex));
        }
        offset += m_Images[i]->GetLargestPossibleRegion().GetNumberOfPixels();


        std::stringstream name;
        name<<"Y"<<i<<".nii.gz";

        btk::ImageHelper< ImageType >::WriteImage(m_ImagesFilledWithY[i], name.str());
    }



}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void
SRHMatrixComputation< TImage >::TestFillingOfX()
{
    typename ImageType::Pointer ImageX = btk::ImageHelper< ImageType>::CreateNewImageFromPhysicalSpaceOf(m_ReferenceImage);
    ImageX->FillBuffer(0.0);
    ConstIteratorType it(m_ReferenceImage, m_ReferenceImage->GetLargestPossibleRegion());
    unsigned int linearIndex = 0;
    m_X.set_size(m_ReferenceImage->GetLargestPossibleRegion().GetNumberOfPixels());
    // Fill X
    for(it.GoToBegin(); !it.IsAtEnd(); ++it, linearIndex++)
    {
        assert(linearIndex < UINT_MAX);
        m_X(linearIndex) = it.Get();
    }

    // Save X image
    IteratorType xit(ImageX, ImageX->GetLargestPossibleRegion());

    unsigned int xlinearIndex = 0;
    for(xit.GoToBegin(); !xit.IsAtEnd(); ++xit, xlinearIndex++)
    {
        assert(xlinearIndex< UINT_MAX);
        xit.Set(m_X(xlinearIndex));
    }

    btk::ImageHelper< ImageType >::WriteImage(ImageX, "XImage.nii.gz");

}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void
SRHMatrixComputation< TImage >::SimulateY()
{
    if(!m_IsHComputed )
    {
        this->Update();
    }

    vnl_vector< PrecisionType > SimY;
    SimY.set_size(m_Y->size());

    m_H->mult(m_X, SimY);

    // H should previoulsy be computed
    m_SimulatedImages.resize(m_NumberOfLRImages);
    unsigned int offset = 0;

    for(unsigned int i = 0; i< m_NumberOfLRImages; i++)
    {
        SizeType lrSize = m_Images[i]->GetLargestPossibleRegion().GetSize();
        m_SimulatedImages[i] = ImageType::New();
        m_SimulatedImages[i] = btk::ImageHelper< ImageType >::CreateNewImageFromPhysicalSpaceOf(m_Images[i]);
        m_SimulatedImages[i]->FillBuffer(0.0);

        IteratorType it(m_SimulatedImages[i], m_SimulatedImages[i]->GetLargestPossibleRegion());
        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            IndexType lrIndex = it.GetIndex();

            unsigned int lrLinearIndex = offset + lrIndex[0] +lrIndex[1]*lrSize[0]+ lrIndex[2]*lrSize[0]*lrSize[1];

            it.Set(SimY(lrLinearIndex));
        }
        offset += m_Images[i]->GetLargestPossibleRegion().GetNumberOfPixels();


        std::stringstream name;
        name<<"SimY"<<i<<".nii.gz";

        btk::ImageHelper< ImageType >::WriteImage(m_SimulatedImages[i], name.str());
    }
}
}//namespace
