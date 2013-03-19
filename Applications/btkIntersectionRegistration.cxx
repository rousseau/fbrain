/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 05/11/2012
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


/* ITK */
#include "itkImage.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkExtractImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegistrationMethod.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRigid3DTransform.h"
#include "itkRigid2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkHistogramMatchingImageFilter.h"

/* BTK */
#include "btkImageHelper.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkMotionCorrectionByIntersection.h"
#include "btkIOTransformHelper.h"
#include "btkResampleImageByInjectionFilter.h"
#include "btkSliceBySliceTransformBase.h"
#include "btkFileHelper.h"
#include "btkLowToHighResolutionFilter.hxx"
#include "btkCenteredEulerSliceBySliceTransform.h"
#include "btkSuperResolutionRigidImageFilter.h"
#include "btkNLMTool.h"
#include "btkMRIBiasCorrectionFilter.h"

/* OTHERS */
#include "iostream"
#include "sstream"
#include <tclap/CmdLine.h>



int main(int argc, char * argv[])
{

    const unsigned int Dimension = 3;
    typedef float PixelType;
    //typedef  double PixelType;
    typedef itk::Image< PixelType, Dimension > itkImage;
    typedef itk::Image< unsigned char, Dimension> itkMaskImage;
    typedef itk::ImageMaskSpatialObject<Dimension> Mask;
    typedef itk::MatrixOffsetTransformBase<double, Dimension > TransformBase;
    typedef btk::EulerSliceBySliceTransform< double, Dimension, PixelType > Transform;
    typedef btk::LowToHighResolutionFilter<itkImage> HighResFilter;
    typedef btk::ResampleImageByInjectionFilter< itkImage, itkImage>  ResampleFilter;
    typedef itk::ImageMaskSpatialObject< Dimension >  MaskType;
    typedef btk::SliceBySliceTransformBase< double, Dimension, PixelType > TransformBaseType;
    //typedef itk::MatrixOffsetTransformBase<double , Dimension> TransformBaseType;
    //typedef btk::SliceBySliceTransform<double,Dimension, PixelType> Transform;
    //typedef btk::CenteredEulerSliceBySliceTransform<double, Dimension, PixelType> Transform;
    typedef itk::ImageMaskSpatialObject< Dimension > MaskType;
    typedef itk::HistogramMatchingImageFilter<itkImage, itkImage, PixelType> HistoFilter;


    // TCLAP :
    TCLAP::CmdLine cmd("Apply the reconstruction with the intersection based method", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","mask image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> tranArg("t","transform","transforms to write",true,"string",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","output","reconstructed image",true,"","string",cmd);

    TCLAP::SwitchArg  VerboseArg("v","verbose","verbose Mode", cmd, false);
    TCLAP::SwitchArg  DenoisingArg("","denoising","denoising inputs", cmd, false);
    TCLAP::SwitchArg  BiasArg("","bias","bias correction", cmd, false);
    TCLAP::SwitchArg  InverseArg("","inverse","use inverse transformations", cmd, false);
    TCLAP::ValueArg<int> LoopArg("l","loop","Number of loops",false,1,"int",cmd);
    TCLAP::ValueArg<int> IterArg("","iter","Number of iteration for SR",false,1,"int",cmd);


    std::vector< std::string > input;
    std::vector< std::string > mask;
    std::vector<std::string> transfoNames;
    std::string output;

    std::vector< itkImage::Pointer > inputsImages;
    std::vector<MaskType::Pointer> masks;
    std::vector<itkMaskImage::Pointer> inputMasks;


    typedef btk::SuperResolutionRigidImageFilter< itkImage, itkImage >  SuperResolutionFilter;
    SuperResolutionFilter::Pointer SR_filter = SuperResolutionFilter::New();


    // Parse the argv array.
    cmd.parse( argc, argv );
    input = inputArg.getValue();
    mask = maskArg.getValue();
    transfoNames = tranArg.getValue();
    output = outputArg.getValue();
    bool verboseMode = VerboseArg.getValue();
    bool denoisingInput = DenoisingArg.getValue();
    bool biasInput = BiasArg.getValue();
    bool UseInverse = InverseArg.getValue();
    inputsImages = btk::ImageHelper<itkImage>::ReadImage(input);
    inputMasks = btk::ImageHelper<itkMaskImage>::ReadImage(mask);
    masks.resize(inputMasks.size());

    int loop = LoopArg.getValue();
    int numberOfLoops = IterArg.getValue();

    std::vector<TransformBase::Pointer> T;
    std::vector<Transform::Pointer> Identity;
    T.resize(inputMasks.size());
    Identity.resize(inputMasks.size());

    std::vector< Transform::Pointer > transforms;
    transforms.resize(inputsImages.size());

    bool computeRegistration = true;
    int ReferenceImageMatching = 0;
    bool histogramMatching = true;

    for(int i = 0; i< inputsImages.size(); i++)
    {
        T[i] = TransformBase::New();
        T[i]->SetIdentity();
        Identity[i] = Transform::New();
        Identity[i]->SetImage(inputsImages[i]);
        Identity[i]->SetIdentity();
        Identity[i]->Initialize();

        transforms[i] = Transform::New();
        transforms[i]->SetImage(inputsImages[i]);
        transforms[i]->SetIdentity();
        transforms[i]->Initialize();

        MaskType::Pointer mask = MaskType::New();
        mask -> SetImage( inputMasks[i] );
        SR_filter->AddInput(inputsImages[i]);
        SR_filter->AddMask(mask);

        itkImage::RegionType roi = mask -> GetAxisAlignedBoundingBoxRegion();
        SR_filter -> AddRegion( roi );

        if(denoisingInput)
        {
            btk::NLMTool<PixelType> myTool;
            myTool.SetInput(inputsImages[i]);
            myTool.SetPaddingValue(0);
            myTool.SetDefaultParameters();
            myTool.ComputeOutput();
            inputsImages[i] = myTool.GetOutput();
        }
        if(biasInput)
        {
            btk::MRIBiasCorrectionFilter<itkImage, itkImage>::Pointer bias = btk::MRIBiasCorrectionFilter<itkImage, itkImage>::New();
            bias->SetInput(inputsImages[i]);
            bias->SetMask(inputMasks[i]);
            bias->SetBiasFilterType(btk::MRIBiasCorrectionFilter<itkImage, itkImage>::N4_TYPE);
            bias->Update();

            inputsImages[i] = bias->GetOutput();
        }

    }

    if(histogramMatching)
    {
        for(int i = 0; i< inputsImages.size(); i++)
        {
            if(i != ReferenceImageMatching)
            {
                std::cout<<"Histogram matching..."<<std::endl;
                HistoFilter::Pointer filter = HistoFilter::New();
                filter->SetInput(inputsImages[i]);
                filter->SetReferenceImage(inputsImages[ReferenceImageMatching]);
                filter->SetNumberOfHistogramLevels(1024);
                filter->SetNumberOfMatchPoints(7);
                filter->SetThresholdAtMeanIntensity(true);
                filter->Update();
                inputsImages[i] = filter->GetOutput();
            }

        }
    }









    btk::MotionCorrectionByIntersection<itkImage>* IntersectionFilter = new btk::MotionCorrectionByIntersection<itkImage>();
    //---------------------------------------------------------------------
    if(computeRegistration)
    {

        IntersectionFilter->SetImages(inputsImages);
        IntersectionFilter->SetMasks(inputMasks);
        IntersectionFilter->SetVerboseMode(verboseMode);
        IntersectionFilter->SetUseSliceExclusion(false);
        IntersectionFilter->SetMaxLoop(loop);
        IntersectionFilter->Initialize();
        try
        {
            IntersectionFilter->Update();
        }
        catch(itk::ExceptionObject &exp)
        {
            std::cout<<exp<<std::endl;
            return EXIT_FAILURE;
        }

        if(UseInverse)
        {
           transforms = IntersectionFilter->GetInverseTransforms();
        }
        else
        {
           transforms = IntersectionFilter->GetTransforms();
        }



        //return 0;

    }
    else
    {
        //TODO: Test if files exists or not
       transforms =  btk::IOTransformHelper< Transform >::ReadTransform(transfoNames);
       //transforms.resize(inputsImages.size());
       for(int i = 0; i< transforms.size(); i ++)
       {
           //transforms[i] = Transform::New();
           transforms[i]->SetImage(inputsImages[i]);
           transforms[i]->Initialize();
       }

    }



    for(unsigned int i = 0; i< inputsImages.size(); i++)
    {
        for(unsigned int j=0; j< transforms[i] -> GetNumberOfSlices(); j++)
        {
            SR_filter-> SetTransform(i, j, transforms[i] -> GetSliceTransform(j) ) ;
            //std::cout<<"Transform image["<<i<<"] slice["<<j<<"] ->"<<transforms[i] -> GetSliceTransform(j)<<std::endl;
        }
    }

    // Construction of HighResolution image
    std::cout<<"Perform a High Resolution image ..."<<std::endl;
//    std::vector<itkImage::Pointer> images;
//    images.resize(inputsImages.size());
    HighResFilter::Pointer LowToHigh = HighResFilter::New();
    LowToHigh->SetImages(inputsImages);
    LowToHigh->SetMasks(inputMasks);
    LowToHigh->SetTransforms(T);
    try
    {
        LowToHigh -> Update();
    }
    catch(itk::ExceptionObject &exp)
    {
        std::cout<<exp<<std::endl;
        return EXIT_FAILURE;
    }

    //--------
//    std::vector<std::string> ResampledNames(inputsImages.size());
//    ResampledNames[0] = "tmp1.nii.gz";
//    ResampledNames[1] = "tmp2.nii.gz";
//    ResampledNames[2] = "tmp3.nii.gz";

//    typedef itk::ResampleImageFilter<itkImage,itkImage> ResampleType;

//    for(unsigned int i = 0; i< ResampledNames.size(); i++ )
//    {
//        ResampleType::Pointer resampler =  ResampleType::New();
//        resampler -> SetTransform( static_cast<TransformBase*>(transforms[i].GetPointer()) );
//        resampler -> SetInput( inputsImages[i] );
//        resampler -> SetReferenceImage( LowToHigh->GetOutput() );
//        resampler -> SetUseReferenceImage( true );
//        resampler -> SetDefaultPixelValue( 0 );
//        try
//        {
//        resampler -> Update();
//        }
//        catch(itk::ExceptionObject & exp)
//        {
//            throw (exp);
//        }

//        images[i] = resampler->GetOutput();
//    }

//    //-----

//    btk::ImageHelper<itkImage>::WriteImage(images, ResampledNames);





    //-----
    //
    // Injection :

    std::cout<<" Done !"<<std::endl;
    std::cout<<"Resample by Injection..."<<std::endl;
    ResampleFilter::Pointer resampler = ResampleFilter::New();
    for(unsigned int i = 0; i< inputsImages.size(); i++)
    {
        masks[i] = MaskType::New();
        masks[i] -> SetImage( inputMasks[i] );
        resampler->AddInput(inputsImages[i]);
        resampler->AddRegion(masks[i]->GetAxisAlignedBoundingBoxRegion());
        resampler->SetTransform(i,static_cast<TransformBaseType*>(transforms[i].GetPointer()) );

    }

    resampler -> UseReferenceImageOn();
    resampler -> SetReferenceImage( LowToHigh->GetOutput() );
    resampler -> SetImageMask(LowToHigh -> GetMaskCombination());

    try
    {
        //resampler -> Update();
    }
    catch(itk::ExceptionObject &exp)
    {
        std::cout<<exp<<std::endl;
        return EXIT_FAILURE;
    }


    //btk::ImageHelper<itkImage>::WriteImage(resampler->GetOutput(), "TMP_Reconstruction.nii.gz");
    //--------
    //
    //Super Resolution :


    std::cout<<"Performing super resolution"<<std::endl;
    SR_filter -> UseReferenceImageOn();
    SR_filter -> SetReferenceImage( LowToHigh->GetOutput() );
    SR_filter -> SetIterations(25);
    SR_filter -> SetLambda( 0.02 );
    SR_filter -> SetPSF( SuperResolutionFilter::GAUSSIAN );
    //SR_filter->SetOutliers(IntersectionFilter->GetOutliers());//NOT implemented well
    SR_filter -> Update();


    //btk::ImageHelper<itkImage>::WriteImage(Output, output);

    //return 0;
    for (int i=0; i<numberOfLoops; i++)
    {
      std::cout<<"Loop : "<<i+1<<std::endl;

      btk::NLMTool<PixelType> myTool;
      myTool.SetInput(SR_filter -> GetOutput());
      myTool.SetPaddingValue(0);
      myTool.SetDefaultParameters();
      myTool.ComputeOutput();

      SR_filter -> SetReferenceImage( myTool.GetOutput() );
      SR_filter -> Update();
    }
    //NLM denoising desired at the last step if number of loops > 0
    if(numberOfLoops>0)
    {

      btk::NLMTool<PixelType> myTool;
      myTool.SetInput(SR_filter -> GetOutput());
      myTool.SetPaddingValue(0);
      myTool.SetDefaultParameters();
      myTool.ComputeOutput();

      SR_filter -> SetReferenceImage( myTool.GetOutput() );
      SR_filter -> Update();
    }

    btk::NLMTool<PixelType> myTool;
    myTool.SetInput(SR_filter -> GetOutput());
    myTool.SetPaddingValue(0);
    myTool.SetDefaultParameters();
    myTool.ComputeOutput();

    std::cout<<"Done !"<<std::endl;

     //TODO : Cast into unsigned short at the end !
   itkImage::Pointer Output = itkImage::New();
   Output = myTool.GetOutput();

  btk::ImageHelper<itkImage>::WriteImage(Output, output);
  btk::IOTransformHelper< Transform >::WriteTransform(transforms,transfoNames);


    delete IntersectionFilter;

    return EXIT_SUCCESS;
}

