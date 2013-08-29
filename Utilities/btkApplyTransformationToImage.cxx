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
#include "itkResampleImageFilter.h"

#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkTransformFactoryBase.h"
#include "itkTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkAffineTransform.h"
#include "itkRigid2DTransform.h"
#include "itkRigid3DTransform.h"
#include "itkEuler2DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkMatrixOffsetTransformBase.h"

/* BTK */

#include "btkMacro.h"
#include "btkSliceBySliceTransform.h"
#include "btkImageHelper.h"
#include "btkFileHelper.h"
#include "btkSliceBySliceTransformBase.h"
#include "btkAffineSliceBySliceTransform.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkCenteredEulerSliceBySliceTransform.h"
#include "btkIOTransformHelper.h"
#include "btkApplyTransformToImageFilter.h"



/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>

int main (int argc, char* argv[])
{

    /* Typedefs */

    const   unsigned int    Dimension = 3;
    typedef itk::Transform< double, Dimension > TransformType;
    typedef itk::TransformFileReader     TransformReaderType;
    typedef TransformReaderType::TransformListType* TransformListType;

    typedef float PixelType;
    typedef itk::Image<PixelType, Dimension> itkImage;
    typedef itk::Image<PixelType,2> itkImage2D;

    typedef itk::MatrixOffsetTransformBase<double,Dimension,Dimension> MatrixTransformType;
    typedef itk::AffineTransform<double,Dimension>     itkAffineTransformation;
    typedef itk::Euler3DTransform<double>              itkEulerTransformation;
    typedef btk::SliceBySliceTransformBase< double, Dimension>  SliceBySliceTransfomType;
    typedef btk::AffineSliceBySliceTransform< double, Dimension>  btkAffineSliceBySliceTransform;
    typedef btk::EulerSliceBySliceTransform< double, Dimension>  btkEulerSliceBySliceTransform;
    typedef btk::SliceBySliceTransform<double, Dimension>   btkOldSliceBySliceTransform;
    typedef btk::CenteredEulerSliceBySliceTransform<double , Dimension> btkCenteredEulerSliceBySliceTransform;

    typedef btk::ApplyTransformToImageFilter<itkImage, itkImage> Resampler;





    // TCLAP :
    TCLAP::CmdLine cmd("Apply the new Super-Resolution pipeline. Testing Version!", ' ', "Unversioned");
    TCLAP::ValueArg<int> dimArg("d","dimension","Dimension of input image",false,3,"int",cmd);
    TCLAP::ValueArg<std::string> inputArg("i","input","Input Image file",true,"","string",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","Reference Image file",true,"","string",cmd);
    TCLAP::ValueArg<std::string> transfoArg("t","transfo","Transform file",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",true,"","string",cmd);
    TCLAP::SwitchArg invArg  ("","inv","inverse the transformation", cmd, false);

    // Parse the argv array.
    cmd.parse( argc, argv );

    std::string inputFileImage = inputArg.getValue();
    std::string inputFileTransform = transfoArg.getValue();
    std::string outputFileImage = outArg.getValue();
    std::string referenceFileImage = refArg.getValue();

    bool inverseTheTransform = invArg.getValue();

    itk::TransformFactoryBase::RegisterDefaultTransforms();




    itk::TransformFactory< MatrixTransformType >::RegisterTransform();
    itk::TransformFactory< btkEulerSliceBySliceTransform >::RegisterTransform();
    itk::TransformFactory< btkOldSliceBySliceTransform >::RegisterTransform();
    itk::TransformFactory< btkAffineSliceBySliceTransform >::RegisterTransform();
    itk::TransformFactory< btkCenteredEulerSliceBySliceTransform >::RegisterTransform();


    itkImage::Pointer inputImage = btk::ImageHelper<itkImage>::ReadImage(inputFileImage);
    itkImage::Pointer referenceImage = btk::ImageHelper<itkImage>::ReadImage(referenceFileImage);

    MatrixTransformType::Pointer transform;
    MatrixTransformType::Pointer inverse ;


    std::cout<<"Reading transform:"<<inputFileTransform<<std::endl;
    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader -> SetFileName( inputFileTransform );
    transformReader -> Update();
    TransformListType transforms = transformReader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator titr = transforms->begin();


    const char * className = titr->GetPointer()->GetNameOfClass();

    //We should set the correct image if transform are sliceBySlice transforms

    if(strcmp(className,"AffineSliceBySliceTransform") == 0)
    {

        if(inverseTheTransform)
        {
            inverse = btkAffineSliceBySliceTransform::New();
            btkAffineSliceBySliceTransform::Pointer inverseTransform = btkAffineSliceBySliceTransform::New();

            inverseTransform = static_cast<btkAffineSliceBySliceTransform*>(titr->GetPointer());
            inverseTransform->SetImage(inputImage);
            inverseTransform->GetInverse(static_cast<btkAffineSliceBySliceTransform* >(inverse.GetPointer()));
        }
        else
        {
            transform = btkAffineSliceBySliceTransform::New();
            transform = static_cast<btkAffineSliceBySliceTransform*>(titr->GetPointer());
            static_cast<btkAffineSliceBySliceTransform*>(transform.GetPointer())->SetImage(inputImage);
        }

    }
    else if(strcmp(className,"EulerSliceBySliceTransform") == 0 || strcmp(className,"SliceBySliceTransform") == 0)
    {

        if(inverseTheTransform)
        {
            inverse = btkEulerSliceBySliceTransform::New();
            btkEulerSliceBySliceTransform::Pointer inverseTransform = btkEulerSliceBySliceTransform::New();

            inverseTransform = static_cast<btkEulerSliceBySliceTransform*>(titr->GetPointer());
            inverseTransform->SetImage(inputImage);
            inverseTransform->GetInverse(static_cast<btkEulerSliceBySliceTransform* >(inverse.GetPointer()));
        }
        else
        {
            transform = btkEulerSliceBySliceTransform::New();
            transform = static_cast<btkEulerSliceBySliceTransform*>(titr->GetPointer());
            static_cast<btkEulerSliceBySliceTransform*>(transform.GetPointer())->SetImage(inputImage);

        }

    }
    else if(strcmp(className,"CenteredEulerSliceBySliceTransform") == 0 )
    {

        if(inverseTheTransform)
        {
            inverse = btkCenteredEulerSliceBySliceTransform::New();
            btkCenteredEulerSliceBySliceTransform::Pointer inverseTransform = btkCenteredEulerSliceBySliceTransform::New();

            inverseTransform = static_cast<btkCenteredEulerSliceBySliceTransform*>(titr->GetPointer());
            inverseTransform->SetImage(inputImage);
            inverseTransform->GetInverse(static_cast<btkCenteredEulerSliceBySliceTransform* >(inverse.GetPointer()));
        }
        else
        {
            transform = btkCenteredEulerSliceBySliceTransform::New();
            transform = static_cast<btkCenteredEulerSliceBySliceTransform*>(titr->GetPointer());
            static_cast<btkCenteredEulerSliceBySliceTransform*>(transform.GetPointer())->SetImage(inputImage);

        }

    }
    else
    {
        if(inverseTheTransform)
        {
            inverse =  MatrixTransformType::New();
            transform = static_cast<MatrixTransformType*>(titr->GetPointer());
            transform->GetInverse(inverse);
        }
        else
        {
            transform = static_cast<MatrixTransformType*>(titr->GetPointer());

        }

    }


    Resampler::Pointer resampler = Resampler::New();

    resampler->SetInputImage(inputImage);
    resampler->SetReferenceImage(referenceImage);
    if(inverseTheTransform)
    {
        resampler->SetTransform(inverse);
    }
    else
    {
       resampler->SetTransform(transform);
    }

    resampler->Initialize();

    try
    {
      resampler->Update();
    }
    catch(itk::ExceptionObject & excp)
    {
        std::cerr << "Error while saving the resampled image" << std::endl;
        std::cerr << excp << std::endl;
        std::cout << "[FAILED]" << std::endl;
        return EXIT_FAILURE;
    }

    btk::ImageHelper<itkImage>::WriteImage(resampler->GetOutputImage(),outputFileImage);

    return EXIT_SUCCESS;







}
