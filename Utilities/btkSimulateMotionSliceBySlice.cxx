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

/* ITK */
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

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
    TCLAP::ValueArg<std::string> outputArg("o","output","output image.",true,"","string",cmd);
    TCLAP::ValueArg<int> motionLevelArg("l","level","level of motion. (1,2 or 3 the greater)",false,1,"int",cmd);
    TCLAP::ValueArg<double> rotationLevelArg("","rotMax","Rotation max in radian.",false,0.2,"int",cmd);
    TCLAP::ValueArg<double> translationLevelArg("","traMax","Translation max in mm.",false,1.5,"int",cmd);


    typedef itk::ResampleImageFilter< itkImage, itkImage > Resampler;
    typedef btk::EulerSliceBySliceTransform<double,3>      RigidSliceBySliceTransformation;
    typedef itk::Euler3DTransform<double>                  RigidTransformation;

    cmd.parse(argc,argv);

    std::string inputFileName;
    std::string outputFileName;

    inputFileName = inputArg.getValue();
    outputFileName = outputArg.getValue();
    int level = motionLevelArg.getValue();
    double rotMax = rotationLevelArg.getValue();
    double traMax = translationLevelArg.getValue();

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

    default:
        threshold = 8;
        break;

    }

    Resampler::Pointer resampleFilter = Resampler::New();

    itkImage::Pointer inputImage = itkImage::New();
    itkImage::Pointer outputImage = itkImage::New();
    inputImage = btk::ImageHelper<itkImage>::ReadImage(inputFileName);

    std::cout<<"The Program will compute a motion on several slices of input Image"<<std::endl;



    // Initialize the Global Transformation
    RigidSliceBySliceTransformation::Pointer globalTransform = RigidSliceBySliceTransformation::New();


    globalTransform->SetImage(inputImage);
    globalTransform->Initialize();

    //Slice N°i random transformation :
    RigidTransformation::Pointer RandomTransform = RigidTransformation::New();
    srand(time(NULL));

    int NSlices = globalTransform->GetNumberOfSlices();

    for(int i = 0; i<NSlices; i++)
    {

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
                int paramR = rand()%2 ;
                randomParams[paramR] = rand()/(double)RAND_MAX * rotMax;
            }
            //Translation
            else if(rot_trans == 1)
            {
                int paramT = rand()%(2) + 3;
                randomParams[paramT] = rand()/(double)RAND_MAX * traMax;
            }
            //Rotation and translation
            else
            {
                int paramR = rand()%2;
                int paramT = rand()%2 + 3;
                randomParams[paramR] = rand()/(double)RAND_MAX * rotMax;
                randomParams[paramT] = rand()/(double)RAND_MAX * traMax;
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
            randomParams[x] = 0.0;
        }

    }

    resampleFilter->SetOutputOrigin(inputImage->GetOrigin());
    resampleFilter->SetOutputSpacing(inputImage->GetSpacing());
    resampleFilter->SetOutputDirection(inputImage->GetDirection());
    resampleFilter->SetTransform(globalTransform);

    resampleFilter->SetInput(inputImage);
    resampleFilter->SetSize(inputImage->GetLargestPossibleRegion().GetSize());
    resampleFilter->Update();


    outputImage =  resampleFilter->GetOutput();
    btk::ImageHelper< itkImage >::WriteImage(outputImage,outputFileName);

    return EXIT_SUCCESS;





}
