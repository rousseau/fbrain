/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/08/2013
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


#include "btkRegistration.h"
#include "btkRigidRegistration.h"
#include "btkAffineRegistration.h"
#include "btkImageHelper.h"

#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"

const    unsigned int    Dimension = 3;
typedef  unsigned char           PixelType;

typedef itk::Image< PixelType, Dimension >  ImageType;
typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;

typedef itk::SubtractImageFilter<
    ImageType,
    ImageType,
    ImageType > DifferenceFilterType;
typedef itk::ImageRegionConstIteratorWithIndex< ImageType > IteratorType;

static void CreateEllipseImage(ImageType::Pointer image);
static void CreateSphereImage(ImageType::Pointer image, int* _translation);

int main (int, char* [])
{
    std::cout<<"Btk Registration test"<<std::endl;
    // Get the two images
    ImageType::Pointer  fixedImage  = ImageType::New();
    ImageType::Pointer  movingImageR  = ImageType::New();
    ImageType::Pointer movingImageA = ImageType::New();
    ImageType::Pointer resultR = ImageType::New();
    ImageType::Pointer resultA = ImageType::New();

    ImageType::Pointer mask = ImageType::New();
    bool testPassed = false;
    bool testRigid = false;
    bool testAffine = false;

    int t[3]= {25,25,25};
    CreateSphereImage(fixedImage,t);
    int t1[3]= {20,20,20};
    CreateSphereImage(movingImageR,t1);
    CreateEllipseImage(movingImageA);

    mask = btk::ImageHelper< ImageType >::CreateNewImageFromPhysicalSpaceOf(fixedImage);
    mask->FillBuffer(1);

    btk::Registration< ImageType >::Pointer Reg;


    std::cout<<"- Rigid Registration test"<<std::endl;
    Reg = btk::RigidRegistration< ImageType >::New();
    Reg->SetMovingImage(movingImageR);
    Reg->SetFixedImage(fixedImage);
    Reg->SetFixedImageMask(mask);
    Reg->StartRegistration();

    //  A resampling filter is created and the moving image is connected as  its input.

    ResampleFilterType::Pointer resamplerR = ResampleFilterType::New();
    resamplerR->SetInput( movingImageR);
    resamplerR->SetTransform( Reg->GetTransform() );
    resamplerR->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
    resamplerR->SetOutputOrigin(  fixedImage->GetOrigin() );
    resamplerR->SetOutputSpacing( fixedImage->GetSpacing() );
    resamplerR->SetOutputDirection( fixedImage->GetDirection() );
    resamplerR->SetDefaultPixelValue( 0 );
    resamplerR->Update();
    resultR = resamplerR->GetOutput();


    DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

    difference->SetInput1( fixedImage );
    difference->SetInput2( resultR );
    difference->Update();

    ImageType::Pointer diffR = difference->GetOutput();

    IteratorType itR(diffR, diffR->GetLargestPossibleRegion());
    unsigned int error = 0;
    unsigned int nbVoxel = 0;

    for(itR.GoToBegin(); !itR.IsAtEnd(); ++itR)
    {
        error+= itR.Get();
        nbVoxel++;
    }

    std::cout<<"  Difference between rigid result and fixed image : "<<error/nbVoxel<<std::endl;

    if(error/nbVoxel == 0)
    {
        testRigid = true;
        std::cout<<"Rigid Test Passed !"<<std::endl;
    }


    std::cout<<"- Affine Registration test"<<std::endl;
    Reg = btk::AffineRegistration< ImageType >::New();
    Reg->SetMovingImage(movingImageA);
    Reg->SetFixedImage(fixedImage);
    Reg->SetFixedImageMask(mask);
    Reg->StartRegistration();

    //  A resampling filter is created and the moving image is connected as  its input.

    ResampleFilterType::Pointer resamplerA = ResampleFilterType::New();
    resamplerA->SetInput( movingImageA);
    resamplerA->SetTransform( Reg->GetTransform() );
    resamplerA->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
    resamplerA->SetOutputOrigin(  fixedImage->GetOrigin() );
    resamplerA->SetOutputSpacing( fixedImage->GetSpacing() );
    resamplerA->SetOutputDirection( fixedImage->GetDirection() );
    resamplerA->SetDefaultPixelValue( 0 );
    resamplerA->Update();
    resultA = resamplerA->GetOutput();

    DifferenceFilterType::Pointer differenceA = DifferenceFilterType::New();

    differenceA->SetInput1( fixedImage );
    differenceA->SetInput2( resultR );
    differenceA->Update();

    ImageType::Pointer diffA = differenceA->GetOutput();

    IteratorType itA(diffA, diffA->GetLargestPossibleRegion());
    error = 0;
    nbVoxel = 0;

    for(itA.GoToBegin(); !itA.IsAtEnd(); ++itA)
    {
        error+= itA.Get();
        nbVoxel++;
    }

    std::cout<<"  Difference between affine result and fixed image : "<<error/nbVoxel<<std::endl;

    if(error/nbVoxel == 0)
    {
        testAffine = true;
        std::cout<<"Affine Test Passed !"<<std::endl;
    }

    if(!testAffine || !testRigid)
    {
        std::cout<<"Test failed !"<<std::endl;
        return EXIT_FAILURE;
    }


    std::cout<<"Test passed !"<<std::endl;

    return EXIT_SUCCESS;


}

void CreateEllipseImage(ImageType::Pointer image)
{
  typedef itk::EllipseSpatialObject< Dimension >   EllipseType;

  typedef itk::SpatialObjectToImageFilter<
    EllipseType, ImageType >   SpatialObjectToImageFilterType;

  SpatialObjectToImageFilterType::Pointer imageFilter =
    SpatialObjectToImageFilterType::New();

  ImageType::SizeType size;
  size[ 0 ] =  50;
  size[ 1 ] =  50;
  size[ 2 ] =  50;

  imageFilter->SetSize( size );

  ImageType::SpacingType spacing;
  spacing.Fill(1);
  imageFilter->SetSpacing(spacing);

  EllipseType::Pointer ellipse    = EllipseType::New();
  EllipseType::ArrayType radiusArray;
  radiusArray[0] = 5;
  radiusArray[1] = 10;
  radiusArray[2] = 10;
  ellipse->SetRadius(radiusArray);

  typedef EllipseType::TransformType                 TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  TransformType::OutputVectorType  translation;
  TransformType::CenterType        center;

  translation[ 0 ] =  20;
  translation[ 1 ] =  20;
  translation[ 2 ] =  20;
  transform->Translate( translation, false );

  ellipse->SetObjectToParentTransform( transform );

  imageFilter->SetInput(ellipse);

  ellipse->SetDefaultInsideValue(255);
  ellipse->SetDefaultOutsideValue(0);
  imageFilter->SetUseObjectValue( true );
  imageFilter->SetOutsideValue( 0 );

  imageFilter->Update();

  image->Graft(imageFilter->GetOutput());

}

void CreateSphereImage(ImageType::Pointer image, int* _translation)
{
 typedef itk::EllipseSpatialObject< Dimension >   EllipseType;

  typedef itk::SpatialObjectToImageFilter<
    EllipseType, ImageType >   SpatialObjectToImageFilterType;

  SpatialObjectToImageFilterType::Pointer imageFilter =
    SpatialObjectToImageFilterType::New();

  ImageType::SizeType size;
  size[ 0 ] =  50;
  size[ 1 ] =  50;
  size[ 2 ] =  50;

  imageFilter->SetSize( size );

  ImageType::SpacingType spacing;
  spacing.Fill(1);
  imageFilter->SetSpacing(spacing);

  EllipseType::Pointer ellipse    = EllipseType::New();
  EllipseType::ArrayType radiusArray;
  radiusArray[0] = 5;
  radiusArray[1] = 5;
  radiusArray[2] = 5;
  ellipse->SetRadius(radiusArray);

  typedef EllipseType::TransformType                 TransformType;
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();

  TransformType::OutputVectorType  translation;
  TransformType::CenterType        center;

  translation[ 0 ] =  _translation[0];
  translation[ 1 ] =  _translation[1];
  translation[ 2 ] =  _translation[2];
  transform->Translate( translation, false );

  ellipse->SetObjectToParentTransform( transform );

  imageFilter->SetInput(ellipse);

  ellipse->SetDefaultInsideValue(255);
  ellipse->SetDefaultOutsideValue(0);
  imageFilter->SetUseObjectValue( true );
  imageFilter->SetOutsideValue( 0 );

  imageFilter->Update();

  image->Graft(imageFilter->GetOutput());




}
