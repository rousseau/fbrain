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

class SuperResolutionDataManager;

class SuperResolutionTools
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 3>         itkImage;
  typedef itkImage::Pointer                 itkPointer;  
  typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
  
  int                       m_psftype; // 0: 3D interpolated boxcar, 1: 3D oversampled boxcar
  std::vector<itkPointer>   m_PSF;
  
  void InitializePSF(SuperResolutionDataManager & data);

  
};

void SuperResolutionTools::InitializePSF(SuperResolutionDataManager & data)
{
  //Principle: We build the PSF in LR image space (simple boxcar PSF = one anisotropic voxel) which is then interpolated or oversampled in SR space

  m_psftype = 0; //to be set elsewhere
  
  //we assume that the reconstruction HR image is isotropic
  float hrSpacing = data.m_spacing[0];
  
  //set the correct number of PSF (one PSF for one LR image -> this allows us to use images with different LR resolution)
  m_PSF.resize(data.m_inputLRImages.size());
  
  for(uint i=0; i != m_PSF.size(); i++){
    
    // 1- build the boxcar PSF in LR space (one anisotropic voxel)
    itkPointer LRPSF = itkImage::New();
    itkImage::IndexType lrIndex;
    itkImage::SizeType lrSize;
    lrIndex[0] = 0;  lrIndex[1] = 0;  lrIndex[2] = 0;
    lrSize[0] = 1;   lrSize[1] = 1;   lrSize[2] = 1;
    itkImage::SpacingType lrSpacing = data.m_inputLRImages[i]->GetSpacing();

    itkImage::RegionType lrRegion;
    lrRegion.SetSize(lrSize);
    lrRegion.SetIndex(lrIndex);
    LRPSF->SetRegions(lrRegion);
    LRPSF->SetSpacing(lrSpacing);
    LRPSF->Allocate();
    LRPSF->FillBuffer(1.0);
    std::cout<<"LRPSF spacing : "<<lrSpacing[0]<<" "<<lrSpacing[1]<<" "<<lrSpacing[2]<<"\n";
    
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

    
                       
  }
  
  
  switch (m_psftype) {
    case 0:
      std::cout<<"3D interpolated boxcar.\n";
      break;
    case 1:
      std::cout<<"3D oversampled boxcar.\n";
      break;
    default:
      std::cout<<"Invalid choice for the psf building\n"; 
      break;
  }
  
  
}

#endif
