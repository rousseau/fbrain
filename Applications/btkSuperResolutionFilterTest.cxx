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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkAffineTransform.h"

/* BTK */

#include "btkMacro.h"
#include "btkSliceBySliceTransform.h"
#include "btkSuperResolutionType.h"
#include "btkImageHelper.h"
#include "btkFileHelper.h"
#include "btkSliceBySliceTransformBase.h"
#include "btkAffineSliceBySliceTransform.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkSuperResolutionFilter.h"
#include "btkIOTransformHelper.h"

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


    // TCLAP :
    TCLAP::CmdLine cmd("Apply the new Super-Resolution pipeline. Testing Version!", ' ', "Unversioned");

    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","low-resolution image mask file",false,"string",cmd);
    TCLAP::MultiArg<std::string> transArg("t","transform","transform file",false,"string",cmd);
    TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image for initialization. "
                                          "Typically the output of btkImageReconstruction is used." ,false,"","string",cmd);
    TCLAP::ValueArg<int> loopArg  ("","loop","Number of loops (SR/denoising) (default = 5)",false, 5,"int",cmd);
    TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",true,"","string",cmd);

    TCLAP::SwitchArg    ComputeRegSwitchArg("","noReco","No Reconstruction is performed, the reference 3D image and the transforms files must be set", cmd, false);
    TCLAP::ValueArg<int>    TransfosTypeSwitchArg("","transfos","Type of Transformation to perform (0: Affine 3D, 1: Euler3D, 2: Affine SliceBySlice, 3: Euler SliceBySlice... ).",false,3,"int",cmd);
    TCLAP::ValueArg<int>    ReconTypeSwitchArg("","RecoType","Type of reconstruction to perform (0: SR 3D, 1: IBP, 2:... ).",false,0,"int",cmd);

    TCLAP::ValueArg<int> nlmArg           ("n","nlm","Type of filtering during IBP process (0: no filtering (default), 1: error map filtering, 2: current HR image filtering).", false,0,"int",cmd);
    TCLAP::ValueArg<float> betaArg        ("b","beta","Smoothing parameter for NLM filtering (default = 1).", false,1,"float",cmd);
    TCLAP::ValueArg<int> simArg           ("s","sim","Simulation of LR images based on the input HR image and the input LR images (0: no simulation, 1: simulation).", false,0,"int",cmd);
    TCLAP::ValueArg<int> ibpOrderArg      ("","ibpOrder","Order for the B-spline interpolation during image backprojections (0: nearest neighbor, 1: trilinear etc.)", false,5,"int",cmd);
    TCLAP::ValueArg<int> psfArg           ("p","psftype","Type of the PSF (0: interpolated boxcar, 1: oversampled boxcar (default), 2: Gaussian)", false,1,"int",cmd);
    TCLAP::ValueArg<int> medArg           ("","medianIBP","Type of filtering on the error map (0: mean of error maps (default), 1: median)", false,0,"int",cmd);


    //TODO: Add the others used arg

    std::vector< std::string > input;
    std::vector< std::string > mask;
    std::vector< std::string > transform;



    btk::TRANSFORMATION_TYPE transfosType;
    btk::RECONSTRUCTION_TYPE recoType;




    //std::vector< TransformReaderType > transforms;

    std::vector< itkImage::Pointer > inputsLRImages;
    std::vector< ImageMaskType::Pointer >          inputsLRMasks;
    std::vector< itkTransformBase::Pointer >     inputsLRTransfos;
    std::vector< btkAffineSliceBySliceTransform::Pointer >    inputsLRAffineTransfos;
    std::vector< btkEulerSliceBySliceTransform::Pointer >     inputsLREulerTransfos;
    std::string refImage;
    std::string outImage;

    itkImage::Pointer SuperResolutionImage;



    // Parse the argv array.
    cmd.parse( argc, argv );
    input = inputArg.getValue();
    mask = maskArg.getValue();
    transform = transArg.getValue();
    refImage = refArg.getValue().c_str();
    outImage = outArg.getValue().c_str();

    int loops                    = loopArg.getValue();
    int nlm                      = nlmArg.getValue();
    int simulation               = simArg.getValue();
    float beta                   = betaArg.getValue();
    int ibpOrder                 = ibpOrderArg.getValue();
    int psftype                  = psfArg.getValue();
    int medianIBP                = medArg.getValue();
    transfosType = (btk::TRANSFORMATION_TYPE)TransfosTypeSwitchArg.getValue();
    recoType = (btk::RECONSTRUCTION_TYPE)ReconTypeSwitchArg.getValue();


    int numberOfImages = input.size();

    bool computeTransfo = ComputeRegSwitchArg.getValue();

    inputsLRImages.resize(numberOfImages);
    inputsLRMasks.resize(numberOfImages);
    inputsLRAffineTransfos.resize(numberOfImages);
    inputsLREulerTransfos.resize(numberOfImages);
    inputsLRTransfos.resize(numberOfImages);



    itk::TransformFactory< MatrixTransformType >::RegisterTransform();
    itk::TransformFactory< btkAffineSliceBySliceTransform >::RegisterTransform();
    itk::TransformFactory< btkEulerSliceBySliceTransform >::RegisterTransform();
    itk::TransformFactory< btkOldSliceBySliceTransform >::RegisterTransform();

    std::cout<<"Testing SuperResolution Pipeline :"<<std::endl;

    btk::SuperResolutionFilter * SuperResolutionFilter = NULL;
    SuperResolutionFilter =  new btk::SuperResolutionFilter();
    if(SuperResolutionFilter != NULL)
    {

       inputsLRImages = btk::ImageHelper< itkImage >::ReadImage(input);
       inputsLRMasks = btk::ImageHelper< ImageMaskType >::ReadImage(mask);


        //std::vector<btkOldSliceBySliceTransform::Pointer> tr = btk::IOTransformHelper< btkOldSliceBySliceTransform >::ReadTransformArray(transform);
        for(int i=0; i<numberOfImages; i++)
        {
            if(transform.size() > 0)
            {
               bool TExist = btk::FileHelper::FileExist(transform[i]);

               if(!TExist)
               {
                   // If transfo doesn't exist we must compute them and write them
                   computeTransfo = true;
//                   itkAffineTransformation::Pointer affine = itkAffineTransformation::New();
//                   affine->SetIdentity();
//                   itkEulerTransformation::Pointer euler = itkEulerTransformation::New();
//                   euler->SetIdentity();

                   switch(transfosType)
                   {
                   case btk::AFFINE:
                       inputsLRTransfos[i] = itkAffineTransformation::New();

                       break;

                   case btk::EULER_3D:
                       inputsLRTransfos[i] = itkEulerTransformation::New();

                       break;

                   case btk::SLICE_BY_SLICE_AFFINE:
                       inputsLRTransfos[i] = btkAffineSliceBySliceTransform::New();
                       reinterpret_cast< btkAffineSliceBySliceTransform*>(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                       reinterpret_cast<btkAffineSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->Initialize();
                       break;

                   case btk::SLICE_BY_SLICE_EULER:
                       inputsLRTransfos[i] = btkEulerSliceBySliceTransform::New();
                       reinterpret_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                       reinterpret_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->Initialize();
                       break;

                   default:
                       inputsLRTransfos[i] = btkEulerSliceBySliceTransform::New();
                       reinterpret_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                       reinterpret_cast<btkEulerSliceBySliceTransform*>(inputsLRTransfos[i].GetPointer())->Initialize();
                      break;


                   }
               }
               else
               {
                   // Transfos are set we should read them
                   computeTransfo = false;
                   // test transfoHelper:
                   //inputsLRTransfos[i]=btk::IOTransformHelper< btkOldSliceBySliceTransform >::ReadTransform(transform[i]);
                   std::cout<<"Reading transform:"<<transform[i]<<std::endl;
                   TransformReaderType::Pointer transformReader = TransformReaderType::New();
                   transformReader -> SetFileName( transform[i] );
                   transformReader -> Update();
                   TransformListType transforms = transformReader->GetTransformList();
                   TransformReaderType::TransformListType::const_iterator titr = transforms->begin();

                   const char * className = titr->GetPointer()->GetNameOfClass();

                   if(strcmp(className,"MatrixOffsetTransformBase") == 0 || strcmp(className,"AffineTransform") == 0)
                   {
                       inputsLRTransfos[i] = static_cast<itkAffineTransformation*>(titr->GetPointer());
                       //inputsLRTransfos[i]->SetImage(inputsLRImages[i]);
                       //inputsLRTransfos[i]->Initialize(static_cast<itkAffineTransformation*>(titr->GetPointer()));
                   }

                   if(strcmp(className,"Euler3DTransform") == 0)
                   {
                       inputsLRTransfos[i] = static_cast<itkEulerTransformation*>(titr->GetPointer());
                       //inputsLRTransfos[i]->SetImage(inputsLRImages[i]);
                       //inputsLRTransfos[i]->Initialize(static_cast<itkEulerTransformation*>(titr->GetPointer()));

                   }

                   if(strcmp(className,"AffineSliceBySliceTransform") == 0)
                   {
                       inputsLRTransfos[i] = btkAffineSliceBySliceTransform::New();
                       inputsLRTransfos[i] = dynamic_cast< btkAffineSliceBySliceTransform* >(titr->GetPointer());
                       reinterpret_cast<btkAffineSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                   }

                   if(strcmp(className,"EulerSliceBySliceTransform") == 0 || strcmp(className,"SliceBySliceTransform") == 0)
                   {
                       inputsLRTransfos[i] = btkEulerSliceBySliceTransform::New();
                       inputsLRTransfos[i] = static_cast< btkEulerSliceBySliceTransform* >(titr->GetPointer());
                       static_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                       std::cout<<inputsLRTransfos[i]<<std::endl;

                   }
               }
            }
            else
            {
                // If transfo doesn't exist we must compute them and write them
                computeTransfo = true;
//                   itkAffineTransformation::Pointer affine = itkAffineTransformation::New();
//                   affine->SetIdentity();
//                   itkEulerTransformation::Pointer euler = itkEulerTransformation::New();
//                   euler->SetIdentity();

                switch(transfosType)
                {
                case btk::AFFINE:
                    inputsLRTransfos[i] = itkAffineTransformation::New();

                    break;

                case btk::EULER_3D:
                    inputsLRTransfos[i] = itkEulerTransformation::New();

                    break;

                case btk::SLICE_BY_SLICE_AFFINE:
                    inputsLRTransfos[i] = btkAffineSliceBySliceTransform::New();
                    reinterpret_cast< btkAffineSliceBySliceTransform*>(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                    reinterpret_cast<btkAffineSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->Initialize();
                    break;

                case btk::SLICE_BY_SLICE_EULER:
                    inputsLRTransfos[i] = btkEulerSliceBySliceTransform::New();
                    reinterpret_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                    reinterpret_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->Initialize();
                    break;

                default:
                    inputsLRTransfos[i] = btkEulerSliceBySliceTransform::New();
                    reinterpret_cast<btkEulerSliceBySliceTransform* >(inputsLRTransfos[i].GetPointer())->SetImage(inputsLRImages[i]);
                    reinterpret_cast<btkEulerSliceBySliceTransform*>(inputsLRTransfos[i].GetPointer())->Initialize();
                   break;


                }
            }
        }



        bool srImExist = btk::FileHelper::FileExist(refImage);
        if(srImExist)
        {
            SuperResolutionImage = btk::ImageHelper< itkImage >::ReadImage(refImage);
            SuperResolutionFilter->SetImageHR(SuperResolutionImage);
            SuperResolutionFilter->SetComputeRegistration(false);
            //TODO: Verify that transforms are set !
        }
        else
        {
            SuperResolutionFilter->SetComputeRegistration(true);
        }



        if(!computeTransfo)
        {
           SuperResolutionFilter->SetUseMotionCorrection(false);
        }

        SuperResolutionFilter->SetImagesLR(inputsLRImages);
        SuperResolutionFilter->SetImagesMaskLR(inputsLRMasks);
        SuperResolutionFilter->SetTransformationType(transfosType);
        std::cout<<inputsLRTransfos[0]<<std::endl;
        SuperResolutionFilter->SetTransformsLR(inputsLRTransfos);
        //SuperResolutionFilter->SetTransformsLRSbS(inputsLRTransfos);


        if(transfosType == btk::SLICE_BY_SLICE)
        {
            // Since IBP only works with affine transform :
            SuperResolutionFilter->SetReconstructionType(btk::SR);
        }
        else
        {
            // Since IBP only works with affine transform :
            SuperResolutionFilter->SetReconstructionType(btk::SR);//(btk::IBP);
        }

        SuperResolutionFilter->SetParameters(nlm,beta,loops,medianIBP,psftype, 0.1,5,1,1);
        try
        {
            SuperResolutionFilter->Update();
        }
        catch(std::string &exception)
        {
            std::cerr << "Error : "<< exception <<std::endl;
            return EXIT_FAILURE;


        }


        SuperResolutionImage = SuperResolutionFilter->GetOutput();

        std::vector<TransformType::Pointer> transfos = SuperResolutionFilter->GetTransformsLR();

        std::cout<<transfos[0]<<std::endl;

        if(computeTransfo)
        {
            btk::IOTransformHelper< TransformType >::WriteTransform(transfos,transform);
        }


    }

    // Write HR image
    btk::ImageHelper<itkImage>::WriteImage(SuperResolutionImage, outImage);


    delete SuperResolutionFilter;

    return EXIT_SUCCESS;

}
