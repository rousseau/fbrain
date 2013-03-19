/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 13/06/2012
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
#include "btkAffineSliceBySliceTransform.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkImageHelper.h"
#include "btkIOTransformHelper.h"
#include "btkApplyTransformToImageFilter.h"
#include "btkMathFunctions.h"

/* ITK */
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"


/* Others */
#include "iostream"
#include <tclap/CmdLine.h>
#include "stdio.h"


const unsigned int Dimension = 3;
typedef float PixelType;
typedef itk::Image< PixelType, Dimension >  itkImage;

int main(int argc, char* argv[])
{
    // TCLAP :
    TCLAP::CmdLine cmd("SimulateMotionSliceBySlice, Apply transformations on slices of input image.", ' ', "Unversioned");
    TCLAP::ValueArg<std::string> inputArg("i","input","input image.",true,"","string",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","reference image.",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","output","output image.",true,"","string",cmd);
    TCLAP::ValueArg<std::string> transArg("t","transform","output transformation.",true,"","string",cmd);
    TCLAP::ValueArg<int> motionLevelArg("l","level","level of motion. (1,2,3 or 4 for all slices)",false,2,"int",cmd);
    TCLAP::ValueArg<double> rotationLevelArg("","rotMax","Rotation max in degree.",false,5,"int",cmd);
    TCLAP::ValueArg<double> translationLevelArg("","traMax","Translation max in mm.",false,5,"int",cmd);


    typedef btk::EulerSliceBySliceTransform<double,3>      RigidSliceBySliceTransformation;
    typedef itk::Euler3DTransform<double>                  RigidTransformation;
    typedef itk::ImageRegionIteratorWithIndex<itkImage> Iterator;
    typedef itk::LinearInterpolateImageFunction<itkImage> Interpolator;
    typedef itk::ContinuousIndex<itkImage> ContinuousIndex;
    typedef btk::ApplyTransformToImageFilter<itkImage, itkImage > Resampler;

    cmd.parse(argc,argv);

    std::string inputFileName;
    std::string outputFileName;
    std::string transformFileName;
    std::string refFileName;

    inputFileName = inputArg.getValue();
    outputFileName = outputArg.getValue();
    transformFileName = transArg.getValue();
    refFileName = refArg.getValue();
    int level = motionLevelArg.getValue();
    double rotMax = rotationLevelArg.getValue();
    double traMax = translationLevelArg.getValue();

     rotMax = btk::MathFunctions::DegreesToRadians(rotMax);

     double rotMin =  -rotMax;
     double traMin =  -traMax;

    int threshold = 0;
    switch(level)
    {
    case 1:
        threshold = 8;
        break;

    case 2:
        threshold = 5;
        break;
    case 3:
        threshold = 3;
        break;

    case 4:
        threshold = 0;//all slice
        break;

    default:
        threshold = 8;
        break;

    }

    Resampler::Pointer resampleFilter = Resampler::New();

    itkImage::Pointer inputImage = itkImage::New();
    itkImage::Pointer outputImage = itkImage::New();
    itkImage::Pointer referenceImage = itkImage::New();
    inputImage = btk::ImageHelper<itkImage>::ReadImage(inputFileName);
    referenceImage = btk::ImageHelper<itkImage>::ReadImage(refFileName);


    std::cout<<"The Program will compute a motion on several slices of input Image"<<std::endl;



    // Initialize the Global Transformation
    RigidSliceBySliceTransformation::Pointer globalTransform = RigidSliceBySliceTransformation::New();


    globalTransform->SetImage(inputImage);
    globalTransform->Initialize();

    //Slice N°i random transformation :
    RigidTransformation::Pointer RandomTransform = RigidTransformation::New();
    srand(time(NULL));

    int NSlices = globalTransform->GetNumberOfSlices();

    RigidTransformation::ParametersType randomParams(RandomTransform->GetNumberOfParameters());


    for(int i = 0; i<NSlices; i++)
    {
        randomParams = globalTransform->GetSliceParameters(i);
        globalTransform->SetSliceParameters(i, randomParams);


        std::cout<<"Slice N° : "<<i<<std::endl;
        int rndm = rand()%10;
        RandomTransform->SetIdentity();
        RigidTransformation::ParametersType randomParams(RandomTransform->GetNumberOfParameters());

        if(rndm >= threshold)
        {
            std::cout<<"Apply a Motion"<<std::endl;

            int rot_trans = rand()%3;

            //Rotation
            if(rot_trans == 0)
            {
                int paramR = rand()%3 ; // x y and z
                randomParams[paramR] = rand()/((double)RAND_MAX/std::abs(rotMax - rotMin )) - std::abs(rotMin);
            }
            //Translation
            else if(rot_trans == 1)
            {
                //if translation on z wanted change %2 in %3
                int paramT = (rand()%2) + 3; // only x and y          
                randomParams[paramT] = rand()/((double)RAND_MAX/ std::abs(traMax - traMin)) - std::abs(traMin);
            }
            //Rotation and translation
            else
            {
                int paramR = rand()%3;
                //if translation on z wanted change %2 in %3
                int paramT = (rand()%2) + 3;
                randomParams[paramR] = rand()/((double)RAND_MAX/std::abs(rotMax - rotMin )) - std::abs(rotMin);
                randomParams[paramT] = rand()/((double)RAND_MAX/ std::abs(traMax - traMin)) - std::abs(traMin);
            }


            RandomTransform->SetParameters(randomParams);
        }
        else
        {
            std::cout<<"No Motion apply "<<std::endl;
        }

        globalTransform->SetSliceParameters(i, RandomTransform->GetParameters());
        for(int x=0;x<6;x++)
        {
            std::cout<<"Parameters "<<x<<" = "<<randomParams[x]<<std::endl;
            //Erase old parameters

        }
        randomParams = globalTransform->GetSliceParameters(i);

    }


    Resampler::Pointer resampler = Resampler::New();

    resampler->SetInputImage(inputImage);
    resampler->SetReferenceImage(referenceImage);
    resampler->SetTransform(globalTransform);

    try
    {
        resampler->Update();
    }
    catch(itk::ExceptionObject &e)
    {
        std::cout<<e<<std::endl;
    }


    btk::ImageHelper< itkImage >::WriteImage(resampler->GetOutputImage(),outputFileName);
    btk::IOTransformHelper<RigidSliceBySliceTransformation>::WriteTransform(globalTransform,transformFileName);

    return EXIT_SUCCESS;





}
