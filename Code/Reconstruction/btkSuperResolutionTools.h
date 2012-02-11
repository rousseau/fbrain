/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 09/12/2011
  Author(s): François Rousseau (rousseau@unistra.fr)

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

#ifndef __btkSuperResolutionTools_h
#define __btkSuperResolutionTools_h

#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"
#include "itkBSplineInterpolationWeightFunction.h"
#include "itkSubtractImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkAbsoluteValueDifferenceImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "vnl/vnl_sparse_matrix.h"

#include "../Denoising/btkNLMTool.h"


#include <sstream>

class SuperResolutionDataManager;

class SuperResolutionTools
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 3>         itkImage;
  typedef itkImage::Pointer                 itkPointer;  
  typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
  typedef itk::ImageRegionIterator< itkImage > itkIterator;
  typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;
  typedef itk::ContinuousIndex<double,3>     itkContinuousIndex;


  typedef itk::ResampleImageFilter<itkImage, itkImage>                    itkResampleFilter;
  typedef itk::IdentityTransform<double, 3>                               itkIdentityTransform;
  typedef itk::LinearInterpolateImageFunction<itkImage, double>           itkLinearInterpolator;
  typedef itk::BSplineInterpolateImageFunction<itkImage, double, double>  itkBSplineInterpolator;
  typedef itk::SubtractImageFilter <itkImage, itkImage >                  itkSubtractImageFilter;
  typedef itk::AddImageFilter <itkImage, itkImage >                       itkAddImageFilter;
  typedef itk::BinaryThresholdImageFilter <itkImage, itkImage>            itkBinaryThresholdImageFilter;
  typedef itk::AbsoluteValueDifferenceImageFilter <itkImage, itkImage, itkImage>  itkAbsoluteValueDifferenceImageFilter;  
  typedef itk::StatisticsImageFilter<itkImage>                            itkStatisticsImageFilter;


  int                       m_psftype; // 0: 3D interpolated boxcar, 1: 3D oversampled boxcar
  std::vector<itkPointer>   m_PSF;
  int                       m_interpolationOrder;
  vnl_sparse_matrix<float>  m_H;
  vnl_vector<float>         m_Y;
  vnl_vector<float>         m_X;
  float                     m_paddingValue;
  std::vector<unsigned int> m_offset;
  
  
  SuperResolutionTools(){
    m_interpolationOrder = 1;   //linear interpolation for interpolated PSF
    m_psftype = 1;              //interpolated PSF by default (1: oversampled PSF)
    m_paddingValue = 0;         //0 is considered as background by default.
  };
  
  void SetPSFInterpolationOrder(int & order);
  void SetPSFComputation(int & type);
  void InitializePSF(SuperResolutionDataManager & data);
  void HComputation(SuperResolutionDataManager & data);
  void UpdateX(SuperResolutionDataManager & data);
  void SimulateLRImages(SuperResolutionDataManager & data);
  double IteratedBackProjection(SuperResolutionDataManager & data, int & nlm, float & beta);
  void CreateMaskHRImage(SuperResolutionDataManager & data);
};


void SuperResolutionTools::SetPSFInterpolationOrder(int & order)
{
  //order of the BSpline used for interpolated PSF
  m_interpolationOrder = order;
  std::cout<<"Order of the BSpline interpolator: "<<order<<". ";
  std::cout<<"This order is used for PSF computation if m_psftype=0.\n";
}

void SuperResolutionTools::SetPSFComputation(int & type)
{
  // type = 0 : interpolated PSF
  // type = 1 : oversampled PSF
  m_psftype = type;
}

void SuperResolutionTools::InitializePSF(SuperResolutionDataManager & data)
{
  //Principle: We build the PSF in LR image space (simple boxcar PSF = one anisotropic voxel) which is then interpolated or oversampled in SR space
  
  std::cout<<"Initializing the PSF\n";
  
  //set the correct number of PSF (one PSF for one LR image -> this allows us to use images with different LR resolution)
  m_PSF.resize(data.m_inputLRImages.size());
  
  for(uint i=0; i != m_PSF.size(); i++){
    
    // 1- build the boxcar PSF in LR space (one anisotropic voxel)
    itkPointer LRPSF = itkImage::New();
    itkImage::IndexType lrIndex;
    itkImage::SizeType lrSize;
    //We enlarge by 1 voxel the LR PSF by null voxel for proper interpolated values.
    int border = 1; //to get proper interpolated value when close to the image boundary
    lrIndex[0] = 0;  lrIndex[1] = 0;  lrIndex[2] = 0;
    lrSize[0] = 1+2*border;   lrSize[1] = 1+2*border;   lrSize[2] = 1+2*border;
    itkImage::SpacingType lrSpacing = data.m_inputLRImages[i]->GetSpacing();
    
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
    hrSize[0] = (int)ceil(lrSpacing[0] / data.m_spacing[0]) + 2;
    hrSize[1] = (int)ceil(lrSpacing[1] / data.m_spacing[1]) + 2;
    hrSize[2] = (int)ceil(lrSpacing[2] / data.m_spacing[2]) + 2;
    hrIndex[0] = 0;  hrIndex[1] = 0;  hrIndex[2] = 0;
    
    std::cout<<"HR PSF spacing : "<<data.m_spacing[0]<<" "<<data.m_spacing[1]<<" "<<data.m_spacing[2]<<"\n";
    std::cout<<"HR PSF size : "<<hrSize[0]<<" "<<hrSize[1]<<" "<<hrSize[2]<<"\n";
    
    itkImage::RegionType hrRegion;
    hrRegion.SetSize(hrSize);
    hrRegion.SetIndex(hrIndex);
    m_PSF[i] = itkImage::New(); 
    m_PSF[i]->SetRegions(hrRegion);
    m_PSF[i]->SetSpacing(data.m_spacing);
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
    std::cout<<"HRPSF origin : "<<hrOrigin[0]<<" "<<hrOrigin[1]<<" "<<hrOrigin[2]<<"\n";
    std::cout<<"HRPSF physical center : "<<hrPointCenter[0]<<" "<<hrPointCenter[1]<<" "<<hrPointCenter[2]<<"\n";
        
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(m_interpolationOrder);
    bsInterpolator->SetInputImage(LRPSF);
    
    itkIteratorWithIndex itPSF(m_PSF[i],m_PSF[i]->GetLargestPossibleRegion());
    double nbSamples = 20; //parameter for oversampled HR PSF.

    switch (m_psftype) {
      case 0:
        std::cout<<"3D interpolated boxcar using "<<m_interpolationOrder<<" order B-Spline.\n";
        //Loop over voxels of HR PSF
        for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF){
          
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
        for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF){
          
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
              for(x=0; x<nbSamples; x++){

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
      default:
        std::cout<<"Invalid choice for the psf building\n"; 
        break;
    } 
    
    //for memory saving, we limit the number of non-null voxels
    for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF){
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
      if(itPSF.Get()>0){
        hrIndex = itPSF.GetIndex();
        std::cout<<itPSF.Get()<<" ("<<hrIndex[0]<<", "<<hrIndex[1]<<", "<<hrIndex[2]<<") \n";
      }
    
    std::cout<<"Set the HR origin to 0 since the PSF now should be centered with image voxels\n\n";
    hrOrigin[0] = 0;
    hrOrigin[1] = 0;
    hrOrigin[2] = 0;
    m_PSF[i]->SetOrigin(hrOrigin);
  }
}

void SuperResolutionTools::HComputation(SuperResolutionDataManager & data)
{
  std::cout<<"Computing the matrix H (y=Hx) + fill y and x. \n";
  //Principle: for each voxel of the LR images, we compute the influence of each voxel of the PSF (centered at the current LR voxel) and add the corresponding influence value (PSF value * interpolation weight) in the matrix H
 
  // Set size of matrices
  unsigned int ncols = data.m_inputHRImage->GetLargestPossibleRegion().GetNumberOfPixels();
  
  //m_offset is used to fill correctly the vector m_Y with the input LR image values (offset for the linear index)
  m_offset.resize(data.m_inputLRImages.size());
  unsigned int nrows = 0;
  for(unsigned int im = 0; im < data.m_inputLRImages.size(); im++){
    m_offset[im] = nrows;
    nrows += data.m_inputLRImages[im]->GetLargestPossibleRegion().GetNumberOfPixels();
  }
  m_H.set_size(nrows, ncols);
  m_Y.set_size(nrows);
  m_Y.fill(0.0);
  m_X.set_size(ncols);
  m_X.fill(0.0);
  
  //linear index : an integer value corresponding to (x,y,z) triplet coordinates (ITK index)
  uint lrLinearIndex = 0;
  uint hrLinearIndex = 0;

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
  itkImage::SizeType  hrSize  = data.m_inputHRImage->GetLargestPossibleRegion().GetSize();
  
  std::cout<<"loop over LR images\n";
  for(uint i=0; i<data.m_inputLRImages.size(); i++){
    
    //Get the size of the current LR image
    itkImage::SizeType  lrSize  = data.m_inputLRImages[i]->GetLargestPossibleRegion().GetSize();

    //Instantiate an iterator over the current LR image
    itkIteratorWithIndex itLRImage(data.m_inputLRImages[i],data.m_inputLRImages[i]->GetLargestPossibleRegion());
    
    //Instantiate an iterator over the current PSF
    itkIteratorWithIndex itPSF(m_PSF[i],m_PSF[i]->GetLargestPossibleRegion());
    
    //Set the correct direction for the PSF of the current image
    m_PSF[i]->SetDirection(data.m_inputLRImages[i]->GetDirection());
    
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
    for(itLRImage.GoToBegin(); !itLRImage.IsAtEnd(); ++itLRImage){
 
      //Test on padding value (speed-up and keep H as sparse as possible)
      if(itLRImage.Get() > m_paddingValue){
        
        //Coordinate in the current LR image
        lrIndex = itLRImage.GetIndex();
        
        //std::cout<<"current LR index:"<<lrIndex[0]<<" "<<lrIndex[1]<<" "<<lrIndex[2]<<" -------------- \n";
        
        //Compute the corresponding linear index of lrIndex
        lrLinearIndex = m_offset[i] + lrIndex[0] + lrIndex[1]*lrSize[0] + lrIndex[2]*lrSize[0]*lrSize[1];
        
        //Fill m_Y
        m_Y[lrLinearIndex] = itLRImage.Get();
        
        //Change the origin of the PSF so that itLRImage location corresponds to the center of the PSF
        data.m_inputLRImages[i]->TransformIndexToPhysicalPoint(lrIndex,lrPoint);
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
        for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF){
                  
          if(itPSF.Get() > 0){
            
            //Get coordinate in m_PSF
            psfIndex = itPSF.GetIndex();
            
            //Compute the physical point of psfIndex
            m_PSF[i]->TransformIndexToPhysicalPoint(psfIndex,psfPoint);
            
            //Apply estimated affine transform to psfPoint (need to apply the inverse since the transform goes from the HR image to the LR image)           
            transformedPoint = data.m_inverseAffineTransform[i]->TransformPoint(psfPoint);
            
            //Get back to the index in the HR image
            data.m_inputHRImage->TransformPhysicalPointToContinuousIndex(transformedPoint,hrContIndex);
            //std::cout<<"hrContIndex no check: "<<hrContIndex[0]<<" "<<hrContIndex[1]<<" "<<hrContIndex[2]<<" \n ";

            //Check if the continuous index hrContIndex is inside the HR image
            if (data.m_inputHRImage->GetLargestPossibleRegion().IsInside(hrContIndex)) {
              
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
              
              if( (data.m_inputHRImage->GetLargestPossibleRegion().IsInside(bsplineStartIndex)) && (data.m_inputHRImage->GetLargestPossibleRegion().IsInside(bsplineEndIndex)) ){
                
                //Set the support region
                bsplineRegion.SetSize(bsplineSize);
                bsplineRegion.SetIndex(bsplineStartIndex);
                
                //Instantiate an iterator on HR image over the bspline region
                itkIteratorWithIndex itHRImage(data.m_inputHRImage,bsplineRegion);
                
                //linear index of bspline weights
                unsigned int weightLinearIndex = 0;
                
                //Loop over the support region
                for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage){
  
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
  for(uint i = 0; i < m_H.rows(); i++){
    
    double sum = m_H.sum_row(i);
    
    vnl_sparse_matrix<float>::row & r = m_H.get_row(i);
    vnl_sparse_matrix<float>::row::iterator col_iter;
    
    for (col_iter = r.begin(); col_iter != r.end(); ++col_iter)
      (*col_iter).second = (*col_iter).second / sum;
  }
  
  //Fill m_X
  //Instantiate an iterator on HR image
  itkIteratorWithIndex itHRImage(data.m_inputHRImage,data.m_inputHRImage->GetLargestPossibleRegion());
  for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage){
    hrIndex = itHRImage.GetIndex();
    hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];
    m_X[hrLinearIndex] = itHRImage.Get();
  }
}
void SuperResolutionTools::UpdateX(SuperResolutionDataManager & data)
{
  std::cout<<"Update x\n";
  m_X.fill(0.0);
  itkImage::IndexType hrIndex;
  uint hrLinearIndex = 0;
  itkImage::SizeType  hrSize  = data.m_currentHRImage->GetLargestPossibleRegion().GetSize();

  itkIteratorWithIndex itHRImage(data.m_currentHRImage,data.m_currentHRImage->GetLargestPossibleRegion());
  for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage){
    hrIndex = itHRImage.GetIndex();
    hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];
    m_X[hrLinearIndex] = itHRImage.Get();
  }
}
void SuperResolutionTools::SimulateLRImages(SuperResolutionDataManager & data)
{
  std::cout<<"Simulating LR images (= Hx)\n";
  
  //Update x
  UpdateX(data);
  
  //Compute H * x
  vnl_vector<float> Hx;
  m_H.mult(m_X,Hx);
  
  //resize the vector of simulated input LR images
  data.m_simulatedInputLRImages.resize(data.m_inputLRImages.size());
  
  //Temporary variables
  itkImage::IndexType lrIndex;  //index of the current voxel in the LR image
  uint lrLinearIndex = 0;

  for(uint i=0; i<data.m_inputLRImages.size(); i++){
    //duplicate the LR input image into the simulated LR images to keep all header information
    itkDuplicator::Pointer duplicator = itkDuplicator::New();
    duplicator->SetInputImage( data.m_inputLRImages[i] );
    duplicator->Update();
    data.m_simulatedInputLRImages[i] = duplicator->GetOutput();
    data.m_simulatedInputLRImages[i]->FillBuffer(0);
    
    //Get the size of the current LR image
    itkImage::SizeType  lrSize  = data.m_inputLRImages[i]->GetLargestPossibleRegion().GetSize();
    
    //Instantiate an iterator over the current LR image
    itkIteratorWithIndex itLRImage(data.m_simulatedInputLRImages[i],data.m_simulatedInputLRImages[i]->GetLargestPossibleRegion());

    
    //Loop over the voxels of the current LR image
    for(itLRImage.GoToBegin(); !itLRImage.IsAtEnd(); ++itLRImage){
 
      //Coordinate in the current LR image
      lrIndex = itLRImage.GetIndex();
      
      //Compute the corresponding linear index of lrIndex
      lrLinearIndex = m_offset[i] + lrIndex[0] + lrIndex[1]*lrSize[0] + lrIndex[2]*lrSize[0]*lrSize[1];
      
      //Fill the simulated input LR image
      itLRImage.Set(Hx[lrLinearIndex]);
    }
    
    itkStatisticsImageFilter::Pointer statisticsImageFilter = itkStatisticsImageFilter::New ();
    statisticsImageFilter->SetInput(data.m_simulatedInputLRImages[i]);
    statisticsImageFilter->Update();
    std::cout << "Stat of the LR image y=Hx. \nMean: " << statisticsImageFilter->GetMean() << std::endl;
    std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
    std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
    std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;        

  }
}

double SuperResolutionTools::IteratedBackProjection(SuperResolutionDataManager & data, int & nlm, float & beta)
{
  std::cout<<"Do iterated back projection\n";
  std::cout<<"NLM Filtering type: "<<nlm<<"\n";
  
  std::cout<<"Compute y-Hx\n";
  SimulateLRImages(data);
  
  std::cout<<"Update the current HR image\n";
  //Initialize to 0 the output HR image
  data.m_outputHRImage->FillBuffer(0);
    
  //parameters for interpolation (bspline interpolator)
  int interpolationOrder = 5;
  itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
  bsInterpolator->SetSplineOrder(interpolationOrder);
  
  std::string s = "";
  
  for(uint i=0; i<data.m_inputLRImages.size(); i++){

    //compute the difference between LR input and the simulated LR image
    itkSubtractImageFilter::Pointer subtractFilter = itkSubtractImageFilter::New ();
    subtractFilter->SetInput1(data.m_inputLRImages[i]);
    subtractFilter->SetInput2(data.m_simulatedInputLRImages[i]);
    subtractFilter->Update();

    //std::ostringstream oss ;
    //oss << i+1 ;
    //s = "ibp_"+oss.str()+"_substract.nii.gz";
    //data.WriteOneImage(subtractFilter->GetOutput(), s);
    
    //interpolate the LR difference 
    itkResampleFilter::Pointer resample = itkResampleFilter::New();
    resample->SetTransform(data.m_affineTransform[i]);
    resample->SetInterpolator(bsInterpolator);
    resample->UseReferenceImageOn();
    resample->SetReferenceImage(data.m_currentHRImage);
    resample->SetInput(subtractFilter->GetOutput());
    
    //s = "ibp_"+oss.str()+"_resample.nii.gz";
    //data.WriteOneImage(resample->GetOutput(), s);

    //Add the interpolated differences
    itkAddImageFilter::Pointer addFilter = itkAddImageFilter::New ();
    addFilter->SetInput1(data.m_outputHRImage);
    addFilter->SetInput2(resample->GetOutput());
    addFilter->Update();
    
    data.m_outputHRImage = addFilter->GetOutput();
    
    //s = "ibp_"+oss.str()+"_simulated.nii.gz";
    //data.WriteOneImage(data.m_simulatedInputLRImages[i], s);
    
  }
  
  //Normalize the resampled difference image
  itkIteratorWithIndex itImage(data.m_outputHRImage,data.m_outputHRImage->GetLargestPossibleRegion());
  for(itImage.GoToBegin(); !itImage.IsAtEnd(); ++itImage)
    itImage.Set( itImage.Get() / data.m_inputLRImages.size() );


  
  if(nlm==1){
    std::cout<<"Smooth the error map using the current reconstructed image as reference for NLM filter --------------------------\n";
    btkNLMTool<float> myTool;
    myTool.SetInput(data.m_outputHRImage);
    myTool.SetMaskImage(data.m_maskHRImage);
    myTool.SetDefaultParameters();
    myTool.SetReferenceImage(data.m_currentHRImage);  
    myTool.SetSmoothing(beta, myTool.m_refImage);  
    myTool.SetBlockwiseStrategy(1); //0 pointwise, 1 block, 2 fast block


    myTool.ComputeOutput();
    myTool.GetOutput(data.m_outputHRImage);    
    //s = "ibp_nlm_error.nii.gz";
    //data.WriteOneImage(data.m_outputHRImage, s);     
  }
 
  //Update the HR image correspondly
  itkAddImageFilter::Pointer addFilter2 = itkAddImageFilter::New ();
  addFilter2->SetInput1(data.m_outputHRImage);
  addFilter2->SetInput2(data.m_currentHRImage);
  addFilter2->Update();
  
  data.m_outputHRImage = addFilter2->GetOutput(); 
  //s = "ibp_updated.nii.gz";
  //data.WriteOneImage(data.m_currentHRImage, s);      
  
  if(nlm==2){
    std::cout<<"Smooth the current reconstructed image ------------------ \n";
    btkNLMTool<float> myTool;
    myTool.SetInput(data.m_outputHRImage);
    myTool.SetMaskImage(data.m_maskHRImage);
    myTool.SetDefaultParameters();
    myTool.SetSmoothing(beta, myTool.m_inputImage);
    myTool.SetBlockwiseStrategy(1); //0 pointwise, 1 block, 2 fast block
    myTool.ComputeOutput();
    myTool.GetOutput(data.m_outputHRImage);  
    //s = "ibp_nlm_smooth.nii.gz";
    //data.WriteOneImage(data.m_currentHRImage, s);                
  }  
  
  std::cout<<"Compute the changes between the two consecutive estimates\n";
  itkAbsoluteValueDifferenceImageFilter::Pointer absoluteValueDifferenceFilter = itkAbsoluteValueDifferenceImageFilter::New ();
  absoluteValueDifferenceFilter->SetInput1(data.m_outputHRImage);
  absoluteValueDifferenceFilter->SetInput2(data.m_currentHRImage);
  absoluteValueDifferenceFilter->Update();
  
  double magnitude       = 0.0;
  double numberOfPoints  = 0.0;
  itkPointer  tmpImage = absoluteValueDifferenceFilter->GetOutput();
  itkIterator itTmpImage(tmpImage, tmpImage->GetLargestPossibleRegion());
  itkIterator itMaskHRImage(data.m_maskHRImage,data.m_maskHRImage->GetLargestPossibleRegion());
  
  for(itMaskHRImage.GoToBegin(),itTmpImage.GoToBegin(); !itMaskHRImage.IsAtEnd(); ++itMaskHRImage, ++itTmpImage){
    if(itMaskHRImage.Get() > 0){
      magnitude += itTmpImage.Get();
      numberOfPoints += 1;
    }
  }
  
  double meanMagnitude = magnitude/numberOfPoints;
  std::cout<<"Current mean change: "<<meanMagnitude<<"\n";

  std::cout<<"Copying new estimate to current estimate HR image\n";
  //Not efficient strategy but clearer to understand the code
  itkDuplicator::Pointer duplicator = itkDuplicator::New();
  duplicator->SetInputImage( data.m_outputHRImage );
  duplicator->Update();
  data.m_currentHRImage = duplicator->GetOutput();
  
  itkStatisticsImageFilter::Pointer statisticsImageFilter = itkStatisticsImageFilter::New ();
  statisticsImageFilter->SetInput(data.m_currentHRImage);
  statisticsImageFilter->Update();
  std::cout << "Stat of the current HR estimate: \nMean: " << statisticsImageFilter->GetMean() << std::endl;
  std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
  std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
  std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;  
  
  double manjonCriterion = 0.002*statisticsImageFilter->GetSigma();
  std::cout<<"stopping criterion of manjon et al.: "<<manjonCriterion<<"\n";
  
  if(meanMagnitude < manjonCriterion)
    meanMagnitude = 0;
    
  return meanMagnitude;
}

void SuperResolutionTools::CreateMaskHRImage(SuperResolutionDataManager & data)
{
  std::cout<<"Create Mask HR Image by interpolating Masks of LR Images\n";
  
  //Initialize to 0 the output HR image
  data.m_maskHRImage->FillBuffer(0);
  
  //parameters for interpolation (bspline interpolator)
  int interpolationOrder = 0;
  itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
  bsInterpolator->SetSplineOrder(interpolationOrder);
  
  for(uint i=0; i<data.m_inputLRImages.size(); i++){
  
    //interpolate the LR mask 
    itkResampleFilter::Pointer resample = itkResampleFilter::New();
    resample->SetTransform(data.m_affineTransform[i]);
    resample->SetInterpolator(bsInterpolator);
    resample->UseReferenceImageOn();
    resample->SetReferenceImage(data.m_currentHRImage);
    resample->SetInput(data.m_inputLRImages[i]);
    
    //Add the interpolated mask
    itkAddImageFilter::Pointer addFilter = itkAddImageFilter::New ();
    addFilter->SetInput1(data.m_maskHRImage);
    addFilter->SetInput2(resample->GetOutput());
    addFilter->Update();
    
    data.m_maskHRImage = addFilter->GetOutput();
  }
  
  itkStatisticsImageFilter::Pointer statisticsImageFilter = itkStatisticsImageFilter::New ();
  statisticsImageFilter->SetInput(data.m_maskHRImage);
  statisticsImageFilter->Update();
  std::cout << "Stat of the HR mask: \nMean: " << statisticsImageFilter->GetMean() << std::endl;
  std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
  std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
  std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;        
  
  std::cout<<"Binarize the HR mask image\n";
  itkBinaryThresholdImageFilter::Pointer thresholdFilter = itkBinaryThresholdImageFilter::New();
  thresholdFilter->SetInput( data.m_maskHRImage );
  thresholdFilter->SetLowerThreshold(0.5);
  thresholdFilter->SetUpperThreshold( statisticsImageFilter->GetMaximum() + 1);
  thresholdFilter->SetInsideValue(1.0);
  thresholdFilter->SetOutsideValue(0.0); 
  thresholdFilter->Update(); 
  data.m_maskHRImage = thresholdFilter->GetOutput();

}

#endif
