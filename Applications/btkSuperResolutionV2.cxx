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
#include "btkPSF.h"
#include "btkGaussianPSF.h"
#include "btkBoxCarPSF.h"
#include "btkSincPSF.h"
#include "btkHybridPSF.h"

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

// Fonction for initialize PSF
btk::PSF::Pointer ChoosePSF(unsigned int _psf)
{
    btk::PSF::Pointer p = NULL;

    switch(_psf)
    {
        case 0:
            p = btk::BoxCarPSF::New();
            break;

        case 1:
            p = btk::GaussianPSF::New();
            break;

        case 2:
            p = btk::SincPSF::New();
            break;

        case 3:
            p = btk::HybridPSF::New();
            break;

        default:
            p = btk::GaussianPSF::New();
            break;
    }

    return p;

}

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
    TCLAP::CmdLine cmd("Apply the new Super-Resolution pipeline !", ' ', "v2");

    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","low-resolution image mask file",true,"string",cmd);
    TCLAP::MultiArg<std::string> transArg("t","transform","transformations",false,"string",cmd);
    TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image for initialization. "
                                          "Typically the output of btkImageReconstruction is used." ,false,"","string",cmd);
    TCLAP::MultiArg<std::string> simArg("s","sim","Simulated Low-Resolution Images",false,"string",cmd);

    TCLAP::ValueArg<std::string> outArg("o","output","output image" ,false,"","string",cmd);


    TCLAP::ValueArg<unsigned int> loopArg("","loop","Number of loop for super-resolution (default 1)" ,false,1,"uint",cmd);

    TCLAP::ValueArg<float> lambdaArg("","lambda","Regularisation term (default 0.01)" ,false,0.01,"float",cmd);

    TCLAP::ValueArg<unsigned int> psfArg("","psf","Psf type -> 0 : BoxCar, 1: Gaussian (default), 2: Sinc, 3 : hybrid (sinc on x & y, gaussian on z)" ,false,1,"uint",cmd);


    std::vector< std::string > input;
    std::vector< std::string > mask;
    std::vector< std::string > transform;
    std::vector< std::string > simulation;

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
    simulation = simArg.getValue();


    unsigned int loop = loopArg.getValue();

    int numberOfImages = input.size();

    unsigned int psf = psfArg.getValue();

    float lambda = lambdaArg.getValue();



    inputsLRImages.resize(numberOfImages);
    inputsLRImagesFloat.resize(numberOfImages);
    inputsLRMasks.resize(numberOfImages);
    transforms.resize(numberOfImages);
    inversetransforms.resize(numberOfImages);

    // Reading datas
    inputsLRImages = btk::ImageHelper< itkImage > ::ReadImage(input);

    inputsLRImagesFloat = btk::ImageHelper< itkFloatImage >::ReadImage(input);

    inputsLRMasks = btk::ImageHelper< ImageMaskType >::ReadImage(mask);

    referenceImage = btk::ImageHelper< itkImage > ::ReadImage(refImage);

    //Super Resolution Filter
    btk::SuperResolutionFilter::Pointer SRFilter = btk::SuperResolutionFilter::New();

    SRFilter->SetImages(inputsLRImages);

    SRFilter->SetReferenceImage(referenceImage);

    SRFilter->SetPSF(ChoosePSF(psf));

    SRFilter->SetLambda(lambda);

    //If  simulation
    if(!simulation.empty())
    {
        SRFilter->ComputeSimulatedImages(true);
    }
    else
    {
        SRFilter->ComputeSimulatedImages(false);
    }

    //if not transforms, we set identity
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
    //if transforms we read it
    else
    {
        //FIXME : When AffineSbs is register first, all sbs transform read are affine !!!!!

        //itk::TransformFactory<btkAffineSliceBySliceTransform>::RegisterTransform();
        //itk::TransformFactory<btkEulerSliceBySliceTransform>::RegisterTransform();

        //only SliceBySliceTransform, or EulerSLiceBySLiceTransform (the names are different, but it the same transforms)
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

            //Reading
            SbsTransforms[i] = btk::IOTransformHelper< btkOldSliceBySliceTransform >::ReadTransform(transform[i]);
            //dynamic_cast<btkEulerSliceBySliceTransform*>(transforms[i].GetPointer())->SetImage(inputsLRImagesFloat[i]);
            //Set the corresponding image in the transform (necessary for slice by slice...)
            SbsTransforms[i]->SetImage(inputsLRImagesFloat[i]);

            SbsTransforms[i]->GetInverse(InverseSbsTransforms[i]);

            //Add transform in SR filter
            SRFilter->AddTransform(SbsTransforms[i].GetPointer());
            //Add inverse too
            SRFilter->AddInverseTransform(InverseSbsTransforms[i].GetPointer());
        }

    }



    std::cout<<"Loop : "<<1<<std::endl;

    SRFilter->SetMasks(inputsLRMasks);
    SRFilter->Initialize();

    //Perform super resolution
    SRFilter->Update();

    //Get Output image
    referenceImage = SRFilter->GetOutput();

    //Denoising
    btk::NLMTool<float>* myTool = new btk::NLMTool<float>();
    myTool->SetInput(referenceImage);
    myTool->SetPaddingValue(0);
    myTool->SetDefaultParameters();
    myTool->ComputeOutput();

    referenceImage = myTool->GetOutput();

    //iterative process
    for(unsigned int i = 1; i< loop; i++)
    {
        std::cout<<"Loop : "<<i+1<<std::endl;

        SRFilter->SetReferenceImage(referenceImage);
        SRFilter->Initialize();
        SRFilter->Update();

        myTool->SetInput(referenceImage);
        myTool->SetPaddingValue(0);
        myTool->SetDefaultParameters();
        myTool->ComputeOutput();

        referenceImage = myTool->GetOutput();
    }


    //if simulation, we write it
    if(!simulation.empty())
    {
       std::cout<<"Simulated Low Resolution..."<<std::endl;
       std::vector< itkImage::Pointer > simImages = SRFilter->GetSimulatedImages();
       btk::ImageHelper< itkImage >::WriteImage(simImages,simulation);
    }

    //Write the final SR image
    btk::ImageHelper< itkImage>::WriteImage(referenceImage,outImage);

    //since btk::NLMTool has no smart pointer
    delete myTool;

    //end
    return EXIT_SUCCESS;

}
