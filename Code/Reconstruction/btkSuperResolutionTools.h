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

  typedef itk::ResampleImageFilter<itkImage, itkImage>  itkResampleFilter;
  typedef itk::IdentityTransform<double, 3> itkIdentityTransform;
  typedef itk::LinearInterpolateImageFunction<itkImage, double> itkLinearInterpolator;
  typedef itk::BSplineInterpolateImageFunction<itkImage, double, double>  itkBSplineInterpolator;
  
  int                       m_psftype; // 0: 3D interpolated boxcar, 1: 3D oversampled boxcar
  std::vector<itkPointer>   m_PSF;
  int                       m_interpolationOrder;
  
  SuperResolutionTools(){
    m_interpolationOrder = 1;   //linear interpolation for interpolated PSF
    m_psftype = 0;              //interpolated PSF by default (1: oversampled PSF)
  };
  
  void SetPSFInterpolationOrder(int & order);
  void SetPSFComputation(int & type);
  void InitializePSF(SuperResolutionDataManager & data);
  
};


void SuperResolutionTools::SetPSFInterpolationOrder(int & order)
{
  //order of the BSpline used for interpolated PSF
  m_interpolationOrder = order;
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
    //modification of the origin of m_PSF so that LRPSF and m_PSF are centered.
    itkImage::PointType hrOrigin;
    hrOrigin[0] = lrPointCenter[0] - hrPointCenter[0];
    hrOrigin[1] = lrPointCenter[1] - hrPointCenter[1];
    hrOrigin[2] = lrPointCenter[2] - hrPointCenter[2];
    m_PSF[i]->SetOrigin(hrOrigin);
    m_PSF[i]->TransformContinuousIndexToPhysicalPoint(hrIndexCenter,hrPointCenter);
    std::cout<<"HRPSF physical center : "<<hrPointCenter[0]<<" "<<hrPointCenter[1]<<" "<<hrPointCenter[2]<<"\n";
    
    
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(m_interpolationOrder);
    bsInterpolator->SetInputImage(LRPSF);
    
    itkIteratorWithIndex itPSF(m_PSF[i],m_PSF[i]->GetLargestPossibleRegion());
    double nbSamples = 10; //parameter for oversampled HR PSF.

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
    
    std::cout<<"Positive values of the PSF for the "<<i+1<<"th LR image : ";
    for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF) 
      if(itPSF.Get()>0)
        std::cout<<itPSF.Get()<<" ";
    std::cout<<"\n";
  }
  
  
  
  
  
  
}

#endif
