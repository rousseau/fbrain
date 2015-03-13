/*
 Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 
 31 january 2014
 rousseau@unistra.fr
 
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
 */

/* Standard includes */
#include <tclap/CmdLine.h>
#include "vector"
#include "sstream"
#include "numeric"      // std::accumulate

/* Itk includes */
#include "itkImage.h"

#include "itkMatrixOffsetTransformBase.h"
#include "itkEuler3DTransform.h"
#include "btkIOTransformHelper.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"

#include "vnl/vnl_sparse_matrix.h"

#include <omp.h>

/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPandoraBoxImageFilters.h"
#include "btkPandoraBoxReconstructionFilters.h"
#include "btkPandoraBoxRegistrationFilters.h"


int main(int argc, char** argv)
{
  try {
    
    TCLAP::CmdLine cmd("Dummy", ' ', "0.1", true);
    
    TCLAP::MultiArg<std::string> inputImageArg        ("i","input","Low-resolution image files",true,"string",cmd);
    TCLAP::MultiArg<std::string> inputMaskArg         ("m","mask","Mask of low-resolution image files",false,"string",cmd);
    TCLAP::ValueArg<std::string> outputHRImageArg     ("o","output","High-resolution image file",true,"","string",cmd);
    TCLAP::ValueArg<std::string> inputHRImageArg      ("","init","Input high-resolution image file (used for initialization)",false,"","string",cmd);
    TCLAP::MultiArg<std::string> input3DtransformArg  ("","t3d","3D affine transforms",false,"string",cmd);
    TCLAP::MultiArg<std::string> inputSBStransformArg ("","sbs","Slice by slice affine transforms",false,"string",cmd);
    TCLAP::MultiArg<float>       spacingArg           ("s","spacing","resolution (spacing)) in mm",true,"float",cmd);
    TCLAP::ValueArg<int>         maxIterSRArg         ("","maxIterSR","maximum number of iterations for super-resolution",false,10,"int",cmd);
    
    
    // Parse the args.
    cmd.parse( argc, argv );
    
    // Get the value parsed by each arg.
    std::vector< std::string > input_image_filenames = inputImageArg.getValue();
    std::vector< std::string > input_mask_filenames  = inputMaskArg.getValue();
    std::string                output_HR_filename    = outputHRImageArg.getValue();
    std::string                input_HR_filename     = inputHRImageArg.getValue();
    std::vector< std::string > input_t3d_filenames   = input3DtransformArg.getValue();
    std::vector< std::string > input_sbs_filenames   = inputSBStransformArg.getValue();
    std::vector<float>         spacing               = spacingArg.getValue();
    
    int maxIterSR = maxIterSRArg.getValue();
    
    typedef itk::Image< float, 3 >                                   itkFloatImage;
    typedef itk::MatrixOffsetTransformBase<double,3,3>               itkTransformType;
    //itk::TransformFactory<itkTransformType>::RegisterTransform();
    typedef itk::ResampleImageFilter<itkFloatImage, itkFloatImage>   itkResampleFilter;
    
    //*******************************************************************************************************
    // READING IMAGE DATA
    //*******************************************************************************************************
    
    //Read LR images
    std::vector< itkFloatImage::Pointer > inputLRImages;
    inputLRImages = btk::ImageHelper< itkFloatImage > ::ReadImage(input_image_filenames);
    for(unsigned int i =0; i<inputLRImages.size(); i++)
      btk::PandoraBoxImageFilters::DisplayImageInfo(inputLRImages[i]);
    
    //Read LR masks
    std::vector< itkFloatImage::Pointer > inputLRMasks;
    if(input_mask_filenames.size() > 0)
      inputLRMasks = btk::ImageHelper< itkFloatImage > ::ReadImage(input_mask_filenames);
    else
      inputLRMasks = btk::ImageHelper< itkFloatImage > ::CreateNewImageFromPhysicalSpaceOf(inputLRImages,1.0);
    
    if(inputLRImages.size() != inputLRMasks.size())
    {
      std::cout<<"Not the same number of input LR images and input masks ! Program stops."<<std::endl;
      std::cout<<"Number of input LR images : "<<inputLRImages.size()<<std::endl;
      std::cout<<"Number of input LR masks  : "<<inputLRMasks.size()<<std::endl;
      exit(1);
    }
    
    //Convert input images into stacks
    std::vector< std::vector<itkFloatImage::Pointer> > inputLRStacks;
    inputLRStacks.resize( inputLRImages.size() );
    
    for(unsigned int i=0; i<inputLRImages.size() ; i++)
      btk::PandoraBoxReconstructionFilters::Convert3DImageToSliceStack(inputLRStacks[i], inputLRImages[i]);
    
    //Convert mask images into stacks
    std::vector< std::vector<itkFloatImage::Pointer> > inputMaskStacks;
    inputMaskStacks.resize( inputLRMasks.size() );
    
    for(unsigned int i=0; i<inputLRMasks.size() ; i++)
      btk::PandoraBoxReconstructionFilters::Convert3DImageToSliceStack(inputMaskStacks[i], inputLRMasks[i]);
    
    //*******************************************************************************************************
    // READING TRANSFORM DATA
    //*******************************************************************************************************
    
    std::cout<<"Reading 2D transforms if any"<<std::endl;
    typedef itk::TransformFileReader     TransformReaderType;
    typedef TransformReaderType::TransformListType * TransformListType;
    
    std::vector< std::vector<itkTransformType::Pointer> > affineSBSTransforms;
    affineSBSTransforms.resize( inputLRImages.size() );
    std::vector< std::vector<itkTransformType::Pointer> > inverseAffineSBSTransforms;
    inverseAffineSBSTransforms.resize( inputLRImages.size() );
    
    if(input_sbs_filenames.size() == inputLRImages.size() )
    {
      for(unsigned int i=0;i<input_sbs_filenames.size();i++)
      {
        std::cout<<"Reading slice by slice transform:"<<input_sbs_filenames[i]<<std::endl;
        TransformReaderType::Pointer transformReader = TransformReaderType::New();
        transformReader -> SetFileName( input_sbs_filenames[i] );
        transformReader -> Update();
        
        TransformListType transforms = transformReader->GetTransformList();
        std::cout<<"List size : "<<transforms->size()<<std::endl;
        
        TransformReaderType::TransformListType::const_iterator it = transforms->begin();
        
        for(it = transforms->begin(); it != transforms->end(); ++it)
        {
          affineSBSTransforms[i].push_back(dynamic_cast< itkTransformType * >( it->GetPointer() ) );
          std::cout<<"transfo :"<< affineSBSTransforms[i].back() <<std::endl;
        }
      }
    }
    else
    {
      std::cout<<"Setting slice by slice transforms with identity " << std::endl;
      for(unsigned int i=0;i< inputLRImages.size();i++)
      {
        unsigned int numberOfSlices = inputLRImages[i]->GetLargestPossibleRegion().GetSize()[2];
        affineSBSTransforms[i].resize(numberOfSlices);
        for(unsigned int j=0; j<numberOfSlices; j++)
        {
          //Set the slice by slice transform to identity
          affineSBSTransforms[i][j] = itkTransformType::New();
          affineSBSTransforms[i][j]->SetIdentity();
        }
      }
    }
    
    //Compute inverse slice-by-slice transform
    for(unsigned int i=0;i< inputLRImages.size();i++)
    {
      unsigned int numberOfSlices = inputLRImages[i]->GetLargestPossibleRegion().GetSize()[2];
      inverseAffineSBSTransforms[i].resize(numberOfSlices);
      for(unsigned int j=0; j<numberOfSlices; j++)
      {
        //Set the slice by slice transform to identity
        inverseAffineSBSTransforms[i][j] = itkTransformType::New();
        inverseAffineSBSTransforms[i][j]->SetCenter( affineSBSTransforms[i][j]->GetCenter() );
        affineSBSTransforms[i][j]->GetInverse( inverseAffineSBSTransforms[i][j] );
      }
    }
    
    //or read 3D transforms if any
    std::vector<itkTransformType::Pointer> affine3DTransforms;
    affine3DTransforms.resize(inputLRImages.size());
    
    if(input_t3d_filenames.size() == inputLRImages.size())
    {
      for(unsigned int i=0;i<input_t3d_filenames.size();i++)
      {
        std::cout<<"Reading 3D transform:"<<input_t3d_filenames[i]<<std::endl;
        affine3DTransforms[i] = btk::IOTransformHelper< itkTransformType >::ReadTransform( input_t3d_filenames[i] );
        std::cout<<"Transform :"  << affine3DTransforms[i]->GetNameOfClass()<<std::endl;
        std::cout<<"Center :"     << affine3DTransforms[i]->GetCenter()<<std::endl;
        std::cout<<"Offset :"     << affine3DTransforms[i]->GetOffset()<<std::endl;
        std::cout<<"Translation :"<< affine3DTransforms[i]->GetTranslation()<<std::endl;
        std::cout<<"Parameters :" << affine3DTransforms[i]->GetParameters()<<std::endl;
        std::cout<<"Matrix : "    << affine3DTransforms[i]->GetMatrix()<<std::endl;
        std::cout<<affine3DTransforms[i];
      }
    }
    else
    {
      for(unsigned int i=0;i< inputLRImages.size();i++)
      {
        //Set 3D transform to identity
        affine3DTransforms[i] = itkTransformType::New();
        affine3DTransforms[i]->SetIdentity();
      }
    }
    
    
    //*******************************************************************************************************
    // READ HR IMAGE IF ANY
    //*******************************************************************************************************
    itkFloatImage::Pointer inputHRImage;
    if(input_HR_filename != "")
      inputHRImage = btk::ImageHelper< itkFloatImage > ::ReadImage(input_HR_filename);
    //else, we should estimate one
    //need to set the size, the spacing, and the origin of the HR image
    
    
    itkFloatImage::Pointer outputHRImage;
    itkFloatImage::Pointer tmpImage;
    
    itkFloatImage::SpacingType outputSpacing;
    switch(spacing.size())
    {
      case 3:
        std::cout<<"Output spacing : "<<spacing[0]<<" "<<spacing[1]<<" "<<spacing[2]<<" mm"<<std::endl;
        for(unsigned int i=0; i<3; i++)
          outputSpacing[i] = spacing[i];
        break;
      case 1:
        std::cout<<"Only one spacing argument has been provided. Isotropic spacing is considered. Spacing : "<<spacing[0]<<" mm"<<std::endl;
        for(unsigned int i=0; i<3; i++)
          outputSpacing[i] = spacing[0];
        break;
      default:
        std::cout<<"Spacing has not been correctly setup. Please provide adequat output spacing."<<std::endl;
        exit(1);
    }
    
    //MULTIRESOLUTION (3 levels : typically 2mm, 1mm, 0.5mm)
    
    //AFFINE & RIGIDE : passer les paramètres
    
    
    //SLICE TO VOLUME REGISTRATION (cost function, optimizer (multi-start))
    //TEST (should we use ITK for registration?)
    
    itkFloatImage::Pointer fixedImage;
    itkFloatImage::SpacingType tmpSpacing;
    tmpSpacing[0] = 2;
    tmpSpacing[1] = 2;
    tmpSpacing[2] = 2;
    btk::PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(tmpImage, inputLRStacks[0]);
    btk::PandoraBoxImageFilters::DisplayImageInfo(tmpImage);
    btk::PandoraBoxImageFilters::ResampleImageUsingSpacing(tmpImage,fixedImage,tmpSpacing,0);
    
    /*
     const Doub ftol = 1e-10;//50;
     Amoeba am(ftol);
     VecDoub point(3);
     point[0] = 0 ; // Translation suivant x
     point[1] = 0 ; //Translation suivant y
     point[2] = 0; //Rotation de centre le milieu de l'image
     Doub del = 10;
     fdc myFDC;
     
     am.minimize(point, del, myFDC);
     */
    
    std::cout<<" TEST Simplex \n";
    
    itkFloatImage::Pointer slice;// = inputLRStacks[0][10];
    btk::PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(slice, inputLRStacks[0]);

    itkFloatImage::Pointer ref3DImage;
    btk::PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(ref3DImage, inputLRStacks[0]);
    std::vector<float> inputParam(6);
    std::vector<float> outputParam(6);
    //btk::PandoraBoxRegistrationFilters::Register3DImages(slice, slice, ref3DImage, ref3DImage, inputParam, outputParam);
    //std::cout<<outputParam[0]<<" "<<outputParam[1]<<" "<<outputParam[2]<<" "<<outputParam[3]<<" "<<outputParam[4]<<" "<<outputParam[5]<<std::endl;
    
    
    //2D-3D ITK : RecursiveGaussianImageFilter: The number of pixels along direction 2 is less than 4. This filter requires a minimum of four pixels along the dimension to be processed.
    
    
    //OUTLIER REJECTION
    
    //BIAS FIELD CORRECTION
    
    //LOOP
    
    
    //btk::PandoraBoxReconstructionFilters::IterativeBackProjection(outputHRImage, outputSpacing, inputLRStacks, inputMaskStacks, affineSBSTransforms, inverseAffineSBSTransforms, maxIterSR);
    //btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, output_HR_filename);
    
    exit(1);
    
    btk::PandoraBoxImageFilters::ResampleImageUsingSpacing(inputLRImages[0],outputHRImage,outputSpacing,3);
    
    //btk::PandoraBoxReconstructionFilters::ImageFusionByScatteredInterpolation(outputHRImage, inputLRStacks, inputMaskStacks, inverseAffineSBSTransforms);
    
    
    
    itkFloatImage::Pointer maskHRImage = btk::ImageHelper< itkFloatImage > ::CreateNewImageFromPhysicalSpaceOf(outputHRImage,1.0);
    btk::PandoraBoxReconstructionFilters::ImageFusionByInjection(outputHRImage, maskHRImage, inputLRStacks, affineSBSTransforms);
    btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, output_HR_filename);
    
    exit(1);
    
    //Estimate PSF
    itkFloatImage::Pointer psfImage;
    btk::PandoraBoxReconstructionFilters::ComputePSFImage(psfImage, outputSpacing, inputLRStacks[0][0]->GetSpacing() );
    
    //Compute the parameters (H,X,Y) of the observation model
    vnl_sparse_matrix<float> H;
    vnl_vector<float>        Y;
    vnl_vector<float>        X;
    btk::PandoraBoxReconstructionFilters::ComputerObservationModelParameters(H, Y, X, tmpImage, inputMaskStacks, inputLRStacks, inverseAffineSBSTransforms, psfImage);
    
    //Simulate observations
    std::vector< std::vector<itkFloatImage::Pointer> > simulatedLRStacks;
    btk::PandoraBoxReconstructionFilters::SimulateObservations(H, X, inputLRStacks, simulatedLRStacks);
    
    btk::PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(tmpImage, simulatedLRStacks[0]);
    btk::ImageHelper<itkFloatImage>::WriteImage(tmpImage, "lr_simu.nii.gz");
    
    btk::PandoraBoxImageFilters::ResampleImageUsingSpacing(tmpImage,outputHRImage,outputSpacing,1);
    btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, "simu.nii.gz");
    
    btk::PandoraBoxReconstructionFilters::ImageFusionByInjection(outputHRImage, maskHRImage, simulatedLRStacks, affineSBSTransforms);
    btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, "simu_inj.nii.gz");
    
    
    std::vector< std::vector<itkFloatImage::Pointer> > diffStacks;
    btk::PandoraBoxReconstructionFilters::ComputeModelError(inputLRStacks, simulatedLRStacks, diffStacks);
    
    btk::PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(tmpImage, diffStacks[0]);
    btk::PandoraBoxImageFilters::ResampleImageUsingSpacing(tmpImage,outputHRImage,outputSpacing,1);
    btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, "diff.nii.gz");
    
    btk::PandoraBoxReconstructionFilters::ImageFusionByInjection(outputHRImage, maskHRImage, diffStacks, affineSBSTransforms);
    btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, "diff_inj.nii.gz");
    
    
    
    
    //btk::PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(tmpImage, diffStacks[0]);
    //outputHRImage->FillBuffer(0.0);
    
    
    
    //btk::ImageHelper<itkFloatImage>::WriteImage(tmpImage, output_HR_filename);
    
    exit(1);
    
    itkTransformType::ParametersType params;
    itkTransformType::InputPointType center;
    itkTransformType::MatrixType     matrix;
    itkTransformType::OutputVectorType translation;
    
    params = affineSBSTransforms[0][0]->GetParameters();
    std::cout<<"params : "<<params<<std::endl;
    center = affineSBSTransforms[0][0]->GetCenter();
    std::cout<<"center : "<<center<<std::endl;
    
    std::cout<<params[0]<<std::endl;
    std::cout<<affineSBSTransforms[0][0]<<std::endl;
    
    itk::Euler3DTransform<double>::Pointer eulerTransform = itk::Euler3DTransform<double>::New();
    std::cout<<eulerTransform<<std::endl;
    
    // constant for converting degrees into radians
    const double dtr = ( vcl_atan(1.0) * 4.0 ) / 180.0;
    
    translation[1] = -5.4;
    
    std::cout<<"-------"<<std::endl;
    std::cout<<eulerTransform->GetParameters()<<std::endl;
    eulerTransform->SetRotation( dtr*90, dtr*0, dtr*0 );
    eulerTransform->SetTranslation( translation );
    std::cout<<eulerTransform->GetParameters()<<std::endl;
    std::cout<<"-------"<<std::endl;
    
    itkTransformType::Pointer totoTransform = itkTransformType::New();
    std::cout<<totoTransform<<std::endl;
    //    totoTransform->SetMatrix( eulerTransform->GetMatrix() );
    //    //totoTransform->SetOffset( eulerTransform->GetOffset() );
    //    totoTransform->SetTranslation( eulerTransform->GetTranslation() );
    //    std::cout<<totoTransform<<std::endl;
    
    std::vector<float> p(6);
    p[0] = 90;
    p[1] = 0;
    p[2] = 0;
    p[3] = 10;
    p[4] = 20;
    p[5] = 30;
    
    center[0] = 1;
    center[1] = 2;
    center[2] = 3;
    totoTransform->SetCenter( center );
    //btk::PandoraBoxTransform::ConvertParametersToRigidMatrix(totoTransform, p);
    std::cout<<totoTransform<<std::endl;
    
    std::vector<float> pp(6);
    //btk::PandoraBoxTransform::ConvertRigidMatrixToParameters(pp, totoTransform);
    std::cout<<pp[0]<<" "<<pp[1]<<" "<<pp[2]<<" "<<pp[3]<<" "<<pp[4]<<" "<<pp[5]<<std::endl;
    
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    
    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
}
