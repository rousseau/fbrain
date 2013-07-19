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
#include "btkImageHelper.h"


#include "iostream"

int main (int, char* [])
{
    std::cout<<"PSF Test"<<std::endl;

    typedef itk::Image< float, 3 > ImageType;

    btk::GaussianPSF::Pointer psfG = btk::GaussianPSF::New();
    btk::BoxCarPSF::Pointer psfBC = btk::BoxCarPSF::New();



    ImageType::Pointer image = ImageType::New();
    ImageType::Pointer outputImG = ImageType::New();
    ImageType::Pointer outputImBC = ImageType::New();



    ImageType::SizeType  size;
    ImageType::SpacingType spacing;
    ImageType::IndexType index, centerI;
    index[0] = 0;
    index[1] = 0;
    index[2] = 0;

    spacing[0] = 1;
    spacing[1] = 1;
    spacing[2] = 1;

    size[0] = 15;
    size[1] = 15;
    size[2] = 15;

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

    //////////////////////////////////////////
    ImageType::SizeType size1;
    size1[0] = 3;
    size1[1] = 3;
    size1[2] = 3;

    psfG->ConstructImage(size1);
    btk::ImageHelper< ImageType >::WriteImage(psfG->GetPsfImage(),"PSFTest.nii.gz");

    /////////////////////////////////////////

    psfBC->SetDirection(image->GetDirection());
    psfBC->SetSpacing(image->GetSpacing());

    image->FillBuffer(0.0);
    image->SetPixel(centerI, 1.0);
    ImageType::PointType point;
    image->TransformIndexToPhysicalPoint(centerI,point);
    psfG->SetCenter(point);
    psfBC->SetCenter(point);


    for(it.GoToBegin(),outIt.GoToBegin(), outIt2.GoToBegin();
        !it.IsAtEnd(),!outIt.IsAtEnd(), !outIt2.IsAtEnd();
        ++it,++outIt,++outIt2)
    {
        ImageType::PointType currentPoint;
        image->TransformIndexToPhysicalPoint(it.GetIndex(),currentPoint);
        double valueG = psfG->Evaluate(currentPoint);
        double valueBC = psfBC->Evaluate(currentPoint);
        outIt.Set((float)valueG);
        outIt2.Set((float)valueBC);
    }




    btk::ImageHelper< ImageType >::WriteImage(image,"Impulsional_Image.nii.gz");
    btk::ImageHelper< ImageType >::WriteImage(outputImG,"Gaussian_PSF.nii.gz");
    btk::ImageHelper< ImageType >::WriteImage(outputImBC,"BoxCar_PSF.nii.gz");
}
