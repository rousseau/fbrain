/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: v1 : 24/05/2011 v2 : 20/03/2013
  Author(s): Estanislao Oubel (oubel@unistra.fr)
             Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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
#include "itkImageMaskSpatialObject.h"

/* BTK */
#include "btkEulerSliceBySliceTransform.h"
#include "btkSimulateLowResolutionImagesFilter.hxx"
#include "btkImageHelper.h"
#include "btkIOTransformHelper.h"
#include "btkRandomSliceBySliceTransformGenerator.h"

/* TCLAP */
#include "tclap/CmdLine.h"



// typedefs
const unsigned int                                                      Dimension = 3;
typedef btk::EulerSliceBySliceTransform< double, Dimension >            TransformType;
typedef float                                                           PixelType;
typedef itk::Image< PixelType, Dimension >                              ImageType;
typedef itk::Image< unsigned char, Dimension >                          ImageMaskType;
typedef itk::ImageMaskSpatialObject< Dimension >                        MaskType;
typedef ImageType::RegionType                                           RegionType;
typedef std::vector< RegionType >                                       RegionArrayType;
typedef btk::SimulateLowResolutionImagesFilter< ImageType, ImageType >  ResamplerType;

int main( int argc, char *argv[] )
{

    try {

        std::vector<std::string> inputFile;
        std::vector<std::string> maskFile;
        std::vector<std::string> transformFile;
        std::vector<std::string> simFile;
        std::string refFile;

        // Parse arguments

        TCLAP::CmdLine cmd("Simulates a low resolution image from a high resolution "
                           "image (reconstructed, super-resolution, or acquired) and a "
                           "transformation between both images (transformation can be set or randomly computed)."
                           , ' ', "Unversioned");

        TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
        TCLAP::MultiArg<std::string> maskArg("m","mask","Low-resolution image mask file.",true ,"string",cmd);
        TCLAP::MultiArg<std::string> outArg  ("o","output","Simulated Images",true,"string",cmd);
        TCLAP::MultiArg<std::string> transArg("t","transform","transform file",false,"string",cmd);
        TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image.Typically the output of btkImageReconstruction is used." ,true,"","string",cmd);

        TCLAP::ValueArg<float> rotArg  ("","rotMax","Maximum angle in degree for random rotation (min = -max)" ,false,5,"float",cmd);
        TCLAP::ValueArg<float> traArg  ("","transMax","Maximum value in milimeter for random translation (min = -max)" ,false,5,"float",cmd);

        TCLAP::ValueArg<unsigned int> levelArg  ("","level","level of motion (1: 1/5 slices, 2: 3/5 slices, 3: 4/5 slices, 4 all slices)" ,false,2,"float",cmd);

        TCLAP::SwitchArg  VerboseArg("v","verbose","verbose Mode", cmd, false);
        TCLAP::SwitchArg  UseTransformArg("","UseTransform","Use existing transformations", cmd, false);


        std::cout<<"Simulate Low resolution image with transformation "<<std::endl;


        // Parse the argv array.
        cmd.parse( argc, argv );

        inputFile = inputArg.getValue();
        maskFile = maskArg.getValue();
        refFile = refArg.getValue();
        simFile = outArg.getValue();
        transformFile = transArg.getValue();

        float rotMax = rotArg.getValue();
        float traMax = traArg.getValue();

        unsigned int level = levelArg.getValue();

        bool UseTransforms = UseTransformArg.getValue();
        bool verbose = VerboseArg.getValue();


        ResamplerType::Pointer resampler = ResamplerType::New();
        // Filter setup

        ImageType::IndexType  roiIndex;
        ImageType::SizeType   roiSize;

        std::vector< ImageType::Pointer > inputImages;
        std::vector< ImageType::Pointer > outputImages;
        std::vector< ImageMaskType::Pointer > maskImages;
        std::vector< TransformType::Pointer> transforms;

        ImageType::Pointer refImage;




        // Reading inputs
        inputImages = btk::ImageHelper< ImageType >::ReadImage(inputFile);
        maskImages = btk::ImageHelper< ImageMaskType >::ReadImage(maskFile);
        refImage = btk::ImageHelper< ImageType >::ReadImage(refFile);

        if(UseTransforms)
        {
          // Read transformation
          transforms = btk::IOTransformHelper<TransformType>::ReadTransform(transformFile);
        }
        else
        {
            transforms.resize(inputImages.size());
        }

        for(unsigned int i =0; i< inputImages.size(); i++)
        {
            if(!UseTransforms)
            {
                //Random Slice by Slice transform Generator
                btk::RandomSliceBySliceTransformGenerator::Pointer  randomTransformGenerator = btk::RandomSliceBySliceTransformGenerator::New();
                randomTransformGenerator->SetImage(inputImages[i]);
                randomTransformGenerator->SetMaxTranslation(traMax);
                randomTransformGenerator->SetMaxRotation(rotMax);
                randomTransformGenerator->SetLevel(level);
                randomTransformGenerator->SetVerboseMode(verbose);
                randomTransformGenerator->Update();
                //transforms[i] = TransformType::New();
                transforms[i] = randomTransformGenerator->GetTransform();
            }



            resampler->AddInput(inputImages[i]);

            MaskType::Pointer mask = MaskType::New();
            mask -> SetImage( maskImages[i] );

            RegionType roi = mask -> GetAxisAlignedBoundingBoxRegion();
            roiIndex = roi.GetIndex();
            roiSize  = roi.GetSize();

            RegionType imageRegion;
            imageRegion.SetIndex(roiIndex);
            imageRegion.SetSize (roiSize);

            transforms[i]->SetImage(inputImages[i]);

            resampler->SetTransform(i,transforms[i]);
        }



        // Set reference image


        resampler -> SetReferenceImage( refImage );

        //Process
        resampler ->Update();

        outputImages.resize(inputImages.size());

        for(unsigned int i = 0; i< inputImages.size(); i++)
        {
            outputImages[i] = resampler->GetSimulatedImage(i);
        }


        btk::ImageHelper< ImageType >::WriteImage(outputImages,simFile);
        if(!UseTransforms)
        {
            btk::IOTransformHelper< TransformType >::WriteTransform(transforms,transformFile);
        }

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;

    }

    return EXIT_SUCCESS;
}

