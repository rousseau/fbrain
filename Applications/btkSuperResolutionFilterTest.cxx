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

/* BTK */
#include "btkSuperResolutionFilter.h"
#include "btkMacro.h"
#include "btkSliceBySliceTransform.h"
#include "btkSuperResolutionType.h"

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

/* OTHERS */
#include "iostream"
#include "tclap/CmdLine.h"

int main(int argc, char * argv[])
{
    /* Typedefs */

    const   unsigned int    Dimension = 3;
    typedef itk::Transform< double, Dimension > TransformType;


    typedef float  PixelType;

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


    typedef itk::MatrixOffsetTransformBase<double,Dimension,Dimension> MatrixTransformType;
    typedef itk::AffineTransform<double,Dimension>     itkAffineTransformation;
    typedef btk::SliceBySliceTransform< double, Dimension>  SliceBySliceTransfomType;

    // TCLAP :
    TCLAP::CmdLine cmd("Apply the new Super-Resolution pipeline. Testing Version!", ' ', "Unversioned");

    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","low-resolution image mask file",false,"string",cmd);
    TCLAP::MultiArg<std::string> transArg("t","transform","transform file",false,"string",cmd);
    TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image for initialization. "
                                          "Typically the output of btkImageReconstruction is used." ,true,"","string",cmd);
    TCLAP::ValueArg<int> loopArg  ("","loop","Number of loops (SR/denoising) (default = 5)",false, 5,"int",cmd);
    TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",true,"","string",cmd);

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



    btk::TRANSFORMATION_TYPE transfoType;
    btk::RECONSTRUCTION_TYPE reconstructionType= btk::IBP;

    bool computeTransfo = false;

    //std::vector< TransformReaderType > transforms;

    std::vector< itkImage::Pointer > inputsLRImages;
    std::vector< ImageMaskType::Pointer >          inputsLRMasks;
    std::vector< itkAffineTransformation::Pointer >     inputsLRAffineTransfos;
    std::vector< SliceBySliceTransfomType::Pointer >     inputsLRSbSTransfos;

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

    int numberOfImages = input.size();

    inputsLRImages.resize(numberOfImages);
    inputsLRMasks.resize(numberOfImages);
    inputsLRAffineTransfos.resize(numberOfImages);
    inputsLRSbSTransfos.resize(numberOfImages);
    //transforms.resize(numberOfImages);



    itk::TransformFactory< MatrixTransformType >::RegisterTransform();
    itk::TransformFactory< SliceBySliceTransfomType >::RegisterTransform();

    std::cout<<"Testing SuperResolution Pipeline :"<<std::endl;

    btk::SuperResolutionFilter * SuperResolutionFilter = NULL;
    SuperResolutionFilter =  new btk::SuperResolutionFilter();
    if(SuperResolutionFilter != NULL)
    {
        for(int i=0; i<numberOfImages; i++)
        {
            // Reading input images :
            std::cout<<"Reading image : "<<input[i].c_str()<<std::endl;
            ImageReaderType::Pointer imageReader = ImageReaderType::New();
            imageReader -> SetFileName( input[i].c_str() );
            imageReader -> Update();
            inputsLRImages[i] = imageReader->GetOutput();

            // Reading input Masks :
            // add region
            if ( mask.size() > 0 )
            {
              std::cout<<"Reading mask image : "<<mask[i].c_str()<<std::endl;
              MaskReaderType::Pointer maskReader = MaskReaderType::New();
              maskReader -> SetFileName( mask[i].c_str() );
              maskReader -> Update();
              inputsLRMasks[i] = maskReader->GetOutput();


            }
            if(transform.size() > 0)
            {



                std::cout<<"Reading Affine Transform : "<<transform[i]<<std::endl;
                TransformReaderType::Pointer reader = TransformReaderType::New();
                reader->SetFileName( transform[i].c_str() );
                reader->Update();
                TransformListType transformList = reader->GetTransformList();
                TransformReaderType::TransformListType::const_iterator titr = transformList->begin();

                const char * className = titr -> GetPointer() -> GetNameOfClass();

                btkCoutVariable(className);
                //TODO: Add two vector of transfo (one for Affine, and the other for SbS)
                if(strcmp(className,"MatrixOffsetTransformBase") == 0)
                {
                    SuperResolutionFilter->SetTransformationType(btk::AFFINE);
                    inputsLRAffineTransfos[i] = static_cast< itkAffineTransformation* >( titr->GetPointer() );
                    transfoType = btk::AFFINE ;



                }
                else
                {
                    if(strcmp(className,"SliceBySliceTransform") == 0)
                    {

                        SuperResolutionFilter->SetTransformationType(btk::SLICE_BY_SLICE);
                        inputsLRSbSTransfos[i] = static_cast< SliceBySliceTransfomType* >( titr->GetPointer() );
                        transfoType = btk::SLICE_BY_SLICE;

                    }
                    else
                    {
                        std::cout << className << " is not a valid transform. SuperResolution only take MatrixOffsetTransformBase and SliceBySlice transforms" << std::endl;
                        return EXIT_FAILURE;
                    }
                }



            }
            else
            {
                // If there are no transform : do the identity and let the algorithm compute them :
                inputsLRAffineTransfos[i] = itkAffineTransformation::New();
                inputsLRSbSTransfos[i] = SliceBySliceTransfomType::New();
                computeTransfo = true;
            }




        }

        ImageReaderType::Pointer imageReader = ImageReaderType::New();
        imageReader -> SetFileName( refImage.c_str() );
        imageReader -> Update();


        SuperResolutionFilter->SetImagesLR(inputsLRImages);
        SuperResolutionFilter->SetImagesMaskLR(inputsLRMasks);
        SuperResolutionFilter->SetTransformationType(transfoType);
        //SuperResolutionFilter->SetReconstructionType(reconstructionType);

        if(transfoType == btk::AFFINE)
        {
            SuperResolutionFilter->SetTransformsLRAffine(inputsLRAffineTransfos);

        }
        else
        {
            SuperResolutionFilter->SetTransformsLRSbS(inputsLRSbSTransfos);
        }

        SuperResolutionFilter->SetImageHR(imageReader->GetOutput());
        SuperResolutionFilter->SetParameters(nlm,beta,loops,medianIBP,psftype,1,5);
        SuperResolutionFilter->Update();

        SuperResolutionImage = SuperResolutionFilter->GetOutput();


    }

    // Write HR image
    btkCoutMacro(Writing Output image);
    typedef itk::ImageFileWriter< itkImage >  WriterType;

    WriterType::Pointer writer =  WriterType::New();
    writer-> SetFileName( outImage );
    writer-> SetInput( SuperResolutionImage );
    writer-> Update();

    delete SuperResolutionFilter;

    return 0;

}
