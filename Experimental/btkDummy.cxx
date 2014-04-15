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
#include "btkIOTransformHelper.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"


/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPandoraBoxImageFilters.h"
#include "btkPandoraBoxReconstructionFilters.h"


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

    typedef itk::Image< float, 3 >                                   itkFloatImage;
    typedef itk::MatrixOffsetTransformBase<double,3,3>               itkTransformType;
    itk::TransformFactory<itkTransformType>::RegisterTransform();
    typedef itk::ResampleImageFilter<itkFloatImage, itkFloatImage>   itkResampleFilter;

    //Read LR images
    std::vector< itkFloatImage::Pointer > inputLRImages;
    inputLRImages = btk::ImageHelper< itkFloatImage > ::ReadImage(input_image_filenames);

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

    std::cout<<"Reading 2D transforms if any"<<std::endl;
    typedef itk::TransformFileReader     TransformReaderType;
    typedef TransformReaderType::TransformListType * TransformListType;

    std::vector< std::vector<itkTransformType::Pointer> > affineSBSTransforms;
    affineSBSTransforms.resize( inputLRImages.size() );

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

    //Read HR image if any
    itkFloatImage::Pointer inputHRImage;
    if(input_HR_filename != "")
      inputHRImage = btk::ImageHelper< itkFloatImage > ::ReadImage(input_HR_filename);
    //else, we should estimate one
    //need to set the size, the spacing, and the origin of the HR image


    itkFloatImage::Pointer outputHRImage;

    itkFloatImage::SpacingType outputSpacing;
    for(unsigned int i=0; i<3; i++)
      outputSpacing[i] = spacing[i];

    btk::PandoraBoxImageFilters::ResampleImageUsingSpacing(inputLRImages[0],outputHRImage,outputSpacing,1);

    outputHRImage->FillBuffer(0.0);
    btk::PandoraBoxReconstructionFilters::ImageFusionByInjection(outputHRImage, inputLRStacks, affineSBSTransforms);


    btk::ImageHelper<itkFloatImage>::WriteImage(outputHRImage, output_HR_filename);
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------

    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
