/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 02/04/2012
  Author(s): Marc Schweitzer (marc.schweitzer@unistra.fr)

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



/* ITK */
#include "itkImage.h"
#include "itkTransform.h"
#include "itkAffineTransform.h"

/* BTK */

#include "btkMacro.h"
#include "btkSliceBySliceTransform.h"
#include "btkImageHelper.h"
#include "btkFileHelper.h"
#include "btkSliceBySliceTransformBase.h"
#include "btkAffineSliceBySliceTransform.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkSuperResolutionFilter.h"
#include "btkIOTransformHelper.h"

#include "btkApplyTransformToImageFilter.h"
#include "btkNLMTool.h"


/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>

int main(int argc, char * argv[])
{

    /* Typedefs */

    const   unsigned int    Dimension = 3;
    typedef itk::Transform< double, Dimension > TransformType;


    typedef float  PixelType;
    //typedef short  PixelType;

    typedef itk::Image< PixelType, Dimension >  itkImage;
    typedef itk::Image< float, Dimension >       itkFloatImage;// For initialization of slice by slice transforms

    typedef itk::Image< unsigned char, Dimension >  ImageMaskType;
    typedef itk::ImageFileReader< ImageMaskType >   MaskReaderType;
    typedef itk::ImageMaskSpatialObject< Dimension > MaskType;

    typedef itkImage::RegionType               RegionType;
    typedef std::vector< RegionType >           RegionArrayType;

    typedef itk::ImageFileReader< itkImage >   ImageReaderType;
    typedef itk::ImageFileWriter< itkImage >   WriterType;

    typedef itk::TransformFileReader     TransformReaderType;
    typedef TransformReaderType::TransformListType* TransformListType;

    typedef itk::Transform<double, Dimension> itkTransformBase;
    typedef itk::MatrixOffsetTransformBase<double,Dimension,Dimension> MatrixTransformType;
    typedef itk::AffineTransform<double,Dimension>     itkAffineTransformation;
    typedef itk::Euler3DTransform<double>              itkEulerTransformation;
    typedef btk::SliceBySliceTransformBase< double, Dimension>  SliceBySliceTransfomType;
    typedef btk::AffineSliceBySliceTransform< double, Dimension>  btkAffineSliceBySliceTransform;
    typedef btk::EulerSliceBySliceTransform< double, Dimension>  btkEulerSliceBySliceTransform;
    typedef btk::SliceBySliceTransform<double, Dimension>   btkOldSliceBySliceTransform;
    typedef btk::ApplyTransformToImageFilter<itkImage, itkImage> Resampler;

    // TCLAP :
    TCLAP::CmdLine cmd("Apply the new Super-Resolution pipeline. Testing Version!", ' ', "Unversioned");

    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","low-resolution image mask file",false,"string",cmd);
    TCLAP::MultiArg<std::string> transArg("t","transform","transformations",false,"string",cmd);
    TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image for initialization. "
                                          "Typically the output of btkImageReconstruction is used." ,false,"","string",cmd);

    TCLAP::ValueArg<std::string> outArg("o","output","output image" ,false,"","string",cmd);


    //TODO: Add the others used arg

    std::vector< std::string > input;
    std::vector< std::string > mask;
    std::vector< std::string > transform;

    //std::vector< TransformReaderType > transforms;

    std::vector< itkImage::Pointer > inputsLRImages;
    std::vector< itkFloatImage::Pointer > inputsLRImagesFloat;
    std::vector< ImageMaskType::Pointer >          inputsLRMasks;
    std::vector< MatrixTransformType::Pointer > transforms;
     std::vector< MatrixTransformType::Pointer > inversetransforms;
    std::string refImage;
    std::string outImage;

    itkImage::Pointer referenceImage;



    // Parse the argv array.
    cmd.parse( argc, argv );
    input = inputArg.getValue();
    mask = maskArg.getValue();
    refImage = refArg.getValue();
    transform = transArg.getValue();
    outImage = outArg.getValue();





    int numberOfImages = input.size();



    inputsLRImages.resize(numberOfImages);
    inputsLRImagesFloat.resize(numberOfImages);
    inputsLRMasks.resize(numberOfImages);
    transforms.resize(numberOfImages);
    inversetransforms.resize(numberOfImages);


    inputsLRImages = btk::ImageHelper< itkImage > ::ReadImage(input);

    inputsLRImagesFloat = btk::ImageHelper< itkFloatImage >::ReadImage(input);

    inputsLRMasks = btk::ImageHelper< ImageMaskType >::ReadImage(mask);

    referenceImage = btk::ImageHelper< itkImage > ::ReadImage(refImage);

    btk::SuperResolutionFilter::Pointer SRFilter = btk::SuperResolutionFilter::New();

    SRFilter->SetImages(inputsLRImages);

    SRFilter->SetReferenceImage(referenceImage);
    SRFilter->ComputeSimulatedImages(true);

    if(transform.empty())
    {
        for(unsigned int i = 0; i< numberOfImages; i++)
        {
            transforms[i] = itkAffineTransformation::New();
            inversetransforms[i] = itkAffineTransformation::New();
            dynamic_cast<itkAffineTransformation*>(transforms[i].GetPointer())->SetIdentity();

            transforms[i]->GetInverse(inversetransforms[i]);

            SRFilter->AddTransform(transforms[i].GetPointer());
            SRFilter->AddInverseTransform(inversetransforms[i].GetPointer());
        }
    }
    else
    {
        //FIXME : When AffineSbs is register first, all sbs transform read are affine !!!!!

        //itk::TransformFactory<btkAffineSliceBySliceTransform>::RegisterTransform();
        //itk::TransformFactory<btkEulerSliceBySliceTransform>::RegisterTransform();
        itk::TransformFactory<btkOldSliceBySliceTransform>::RegisterTransform();
        std::vector< btkOldSliceBySliceTransform::Pointer > SbsTransforms;
        std::vector< btkOldSliceBySliceTransform::Pointer > InverseSbsTransforms;

        typedef itk::TransformFileReader     TransformReaderType;
        typedef TransformReaderType::TransformListType* TransformListType;

        SbsTransforms.resize(numberOfImages);
        InverseSbsTransforms.resize(numberOfImages);

        for(unsigned int i = 0; i< numberOfImages; i++)
        {
            SbsTransforms[i] = btkOldSliceBySliceTransform::New();
            InverseSbsTransforms[i] = btkOldSliceBySliceTransform::New();

            SbsTransforms[i] = btk::IOTransformHelper< btkOldSliceBySliceTransform >::ReadTransform(transform[i]);
            //dynamic_cast<btkEulerSliceBySliceTransform*>(transforms[i].GetPointer())->SetImage(inputsLRImagesFloat[i]);
            SbsTransforms[i]->SetImage(inputsLRImagesFloat[i]);

            SbsTransforms[i]->GetInverse(InverseSbsTransforms[i]);

            SRFilter->AddTransform(SbsTransforms[i].GetPointer());
            SRFilter->AddInverseTransform(InverseSbsTransforms[i].GetPointer());
        }

    }




    SRFilter->SetMasks(inputsLRMasks);
    SRFilter->Initialize();
    SRFilter->Update();

    referenceImage = SRFilter->GetOutput();

    btk::NLMTool<float>* myTool = new btk::NLMTool<float>();
    myTool->SetInput(referenceImage);
    myTool->SetPaddingValue(0);
    myTool->SetDefaultParameters();
    myTool->ComputeOutput();

    referenceImage = myTool->GetOutput();


    for(unsigned int i = 0; i< 0; i++)
    {
        SRFilter->SetReferenceImage(referenceImage);
        SRFilter->Initialize();
        SRFilter->Update();

        myTool->SetInput(referenceImage);
        myTool->SetPaddingValue(0);
        myTool->SetDefaultParameters();
        myTool->ComputeOutput();

        referenceImage = myTool->GetOutput();
    }



    std::vector< itkImage::Pointer > simImages = SRFilter->GetSimulatedImages();

    std::vector< std::string > simNames(3);

    simNames[0] = "Sim1.nii.gz";
    simNames[1] = "Sim2.nii.gz";
    simNames[2] = "Sim3.nii.gz";


    btk::ImageHelper< itkImage >::WriteImage(simImages,simNames);

    btk::ImageHelper< itkImage>::WriteImage(referenceImage,outImage);





    return EXIT_SUCCESS;

}
