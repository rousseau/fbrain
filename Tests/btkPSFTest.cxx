/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 
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

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "btkPSF.h"
#include "btkGaussianPSF.h"
#include "btkBoxCarPSF.h"
#include "btkSincPSF.h"
#include "btkHybridPSF.h"
#include "btkImageHelper.h"

#include "vnl/vnl_sparse_matrix.h"


#include "iostream"

int main (int, char* [])
{
    std::cout<<"PSF Test"<<std::endl;

    typedef itk::Image< float, 3 > ImageType;


    btk::GaussianPSF::Pointer psfG = btk::GaussianPSF::New();
    btk::BoxCarPSF::Pointer psfBC = btk::BoxCarPSF::New();
    btk::SincPSF::Pointer psfSinc = btk::SincPSF::New();
    btk::HybridPSF::Pointer psfHybrid = btk::HybridPSF::New();



    ImageType::Pointer image = ImageType::New();
    ImageType::Pointer outputImG = ImageType::New();
    ImageType::Pointer outputImBC = ImageType::New();
    ImageType::Pointer outputImSinC = ImageType::New();
    ImageType::Pointer outputImH = ImageType::New();



    ImageType::SizeType  size;
    ImageType::SpacingType spacing, lrSpacing;
    ImageType::IndexType index, centerI;
    index[0] = 0;
    index[1] = 0;
    index[2] = 0;

    spacing[0] = 0.74;
    spacing[1] = 0.74;
    spacing[2] = 0.74;

    lrSpacing[0] = 1.0;
    lrSpacing[1] = 1.0;
    lrSpacing[2] = 3.5;


    size[0] = 6;
    size[1] = 6;
    size[2] = 10;

    centerI[0] = size[0] /2;
    centerI[1] = size[1] /2;
    centerI[2] = size[2] /2;


    ImageType::RegionType region;
    region.SetIndex(index);
    region.SetSize(size);

    image->SetSpacing(spacing);
    image->SetRegions(region);

    outputImG->SetSpacing(spacing);
    outputImG->SetRegions(region);

    outputImBC->SetSpacing(spacing);
    outputImBC->SetRegions(region);

    image->Allocate();
    outputImG->Allocate();
    outputImBC->Allocate();

    itk::ImageRegionIteratorWithIndex<ImageType> it(image,image->GetLargestPossibleRegion());
    itk::ImageRegionIteratorWithIndex<ImageType> outIt(outputImG,outputImG->GetLargestPossibleRegion());
    itk::ImageRegionIteratorWithIndex<ImageType> outIt2(outputImBC,outputImBC->GetLargestPossibleRegion());

    psfG->SetDirection(image->GetDirection());
    psfG->SetSpacing(image->GetSpacing());

    psfSinc->SetDirection(image->GetDirection());
    psfSinc->SetSpacing(image->GetSpacing());
    psfSinc->SetSize(size);
    psfSinc->ConstructImage();
    btk::ImageHelper< ImageType >::WriteImage(psfSinc->GetPsfImage(),"PSFSinc.nii.gz");
    //return 0;

    //////////////////////////////////////////
    ImageType::SizeType size1;
    size1[0] = 10;
    size1[1] = 10;
    size1[2] = 10;

    psfG->SetSize(size);
    psfG->ConstructImage();
    btk::ImageHelper< ImageType >::WriteImage(psfG->GetPsfImage(),"PSFGaussian.nii.gz");



    /////////////////////////////////////////

    psfHybrid->SetDirection(image->GetDirection());
    psfHybrid->SetSpacing(image->GetSpacing());
    psfHybrid->SetSize(size);
    psfHybrid->SetXFunction(btk::HybridPSF::SINC);
    psfHybrid->SetYFunction(btk::HybridPSF::SINC);
    psfHybrid->SetZFunction(btk::HybridPSF::GAUSSIAN);

    psfHybrid->ConstructImage();
    btk::ImageHelper< ImageType >::WriteImage(psfHybrid->GetPsfImage(),"PSFHybrid.nii.gz");



    psfBC->SetDirection(image->GetDirection());
    psfBC->SetSpacing(image->GetSpacing());
    psfBC->SetSize(size);
    psfBC->SetLrSpacing(lrSpacing);
    psfBC->ConstructImage();


//    image->FillBuffer(0.0);
//    image->SetPixel(centerI, 1.0);
//    ImageType::PointType point;
//    image->TransformIndexToPhysicalPoint(centerI,point);
//    psfG->SetCenter(point);
//    psfBC->SetCenter(point);


//    for(it.GoToBegin(),outIt.GoToBegin(), outIt2.GoToBegin();
//        !it.IsAtEnd(),!outIt.IsAtEnd(), !outIt2.IsAtEnd();
//        ++it,++outIt,++outIt2)
//    {
//        ImageType::PointType currentPoint;
//        image->TransformIndexToPhysicalPoint(it.GetIndex(),currentPoint);
//        double valueG = psfG->Evaluate(currentPoint);
//        double valueBC = psfBC->Evaluate(currentPoint);
//        outIt.Set((float)valueG);
//        outIt2.Set((float)valueBC);
//    }




//    btk::ImageHelper< ImageType >::WriteImage(image,"Impulsional_Image.nii.gz");
//    btk::ImageHelper< ImageType >::WriteImage(outputImG,"Gaussian_PSF.nii.gz");
    btk::ImageHelper< ImageType >::WriteImage(psfBC->GetPsfImage(),"BoxCar_PSF.nii.gz");
    return 0;
}
