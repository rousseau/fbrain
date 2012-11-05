


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

/* BTK */
#include "btkImageHelper.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkMotionCorrectionByIntersection.h"
#include "btkIOTransformHelper.h"
#include "btkLowToHighResolutionFilter.hxx"
#include "btkResampleImageByInjectionFilter.h"
#include "btkSliceBySliceTransformBase.h"

/* OTHERS */
#include "iostream"
#include "sstream"
#include <tclap/CmdLine.h>



int main(int argc, char * argv[])
{

    const unsigned int Dimension = 3;
    //typedef float PixelType;
    typedef unsigned short PixelType;
    typedef itk::Image< PixelType, Dimension > itkImage;
    typedef itk::Image< unsigned char, Dimension> itkMaskImage;
    typedef itk::ImageMaskSpatialObject<Dimension> Mask;
    typedef itk::MatrixOffsetTransformBase<double, 3 > TransformBase;
    typedef btk::EulerSliceBySliceTransform< double, 3, PixelType > Transform;
    //typedef btk::LowToHighResolutionFilter<itkImage> ResampleFilter;
    typedef btk::ResampleImageByInjectionFilter< itkImage, itkImage>  ResampleFilter;
  typedef itk::ImageMaskSpatialObject< Dimension >  MaskType;
      typedef btk::SliceBySliceTransformBase< double, Dimension > TransformBaseType;



    // TCLAP :
    TCLAP::CmdLine cmd("Apply the reconstruction with the intersection based method", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","mask image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> tranArg("t","transform","transforms to write",true,"string",cmd);
    TCLAP::SwitchArg  VerboseArg("v","verbose","verbose Mode", cmd, false);

    std::vector< std::string > input;
    std::vector< std::string > mask;
    std::vector<std::string> transfoNames;

    std::vector< itkImage::Pointer > inputsImages;
    std::vector<MaskType::Pointer> masks;
    std::vector<itkMaskImage::Pointer> inputMasks;


    std::vector< Transform::Pointer > transforms;

    // Parse the argv array.
    cmd.parse( argc, argv );
    input = inputArg.getValue();
    mask = maskArg.getValue();
    transfoNames = tranArg.getValue();
    bool verboseMode = VerboseArg.getValue();
    inputsImages = btk::ImageHelper<itkImage>::ReadImage(input);
    inputMasks = btk::ImageHelper<itkMaskImage>::ReadImage(mask);
    masks.resize(inputMasks.size());

    //---------------------------------------------------------------------
    btk::MotionCorrectionByIntersection<itkImage> IntersectionFilter;


    IntersectionFilter.SetImages(inputsImages);
    IntersectionFilter.SetMasks(inputMasks);
    IntersectionFilter.SetVerboseMode(verboseMode);
    IntersectionFilter.SetUseSliceExclusion(true);
    IntersectionFilter.Initialize();
    try
    {
        IntersectionFilter.Update();
    }
    catch(itk::ExceptionObject &exp)
    {
        std::cout<<exp<<std::endl;
        return EXIT_FAILURE;
    }


    transforms = IntersectionFilter.GetTransforms();



    btk::IOTransformHelper< Transform >::WriteTransform(transforms,transfoNames);


    // Injection :


    ResampleFilter::Pointer resampler = ResampleFilter::New();
    for(unsigned int i = 0; i< inputsImages.size(); i++)
    {
        masks[i] = MaskType::New();
        masks[i] -> SetImage( inputMasks[i] );
        resampler->AddInput(inputsImages[i]);
        resampler->AddRegion(masks[i]->GetAxisAlignedBoundingBoxRegion());
        resampler->SetTransform(i,dynamic_cast<TransformBaseType*>(transforms[i].GetPointer()) );
    }

    resampler -> Update();

    btk::ImageHelper<itkImage>::WriteImage(resampler->GetOutput(), "TMP_Reconstruction.nii.gz");




    return EXIT_SUCCESS;
}

