#ifndef BTKMOTIONCORRECTIONBYINTERSECTION_TXX
#define BTKMOTIONCORRECTIONBYINTERSECTION_TXX

#include "btkMotionCorrectionByIntersection.h"

namespace btk
{
//-------------------------------------------------------------------------------------------------
template<typename TImage>
MotionCorrectionByIntersection<TImage>::MotionCorrectionByIntersection():m_VerboseMode(true),m_MaxLoop(3),m_VerboseDbg(false)
  ,m_CurrentError(0.0),m_UseSliceExclusion(true),m_NumberOfParameters(6)
{
#ifndef NDEBUG
        m_VerboseDbg = false; //true
#else
        m_VerboseDbg = false;
#endif
}
//-------------------------------------------------------------------------------------------------
template <typename TImage>
void MotionCorrectionByIntersection<TImage>::Initialize()
{
    if(m_Images.empty() || m_Masks.empty() || (m_Images.size() != m_Masks.size()))
    {
        btkException("Missing Images or Masks");
    }
    else
    {
        m_NumberOfImages = m_Images.size();
        m_Transforms.resize(m_NumberOfImages);
        m_InverseTransforms.resize(m_NumberOfImages);
        m_BestError.resize(m_NumberOfImages);
        m_BlurredImages.resize(m_NumberOfImages);


        for(int i = 0; i< m_NumberOfImages; i++)
        {
            m_Transforms[i] = TransformType::New();
            m_InverseTransforms[i] = TransformType::New();


            typename ImageType::SizeType size = m_Images[i]->GetLargestPossibleRegion().GetSize();

            // For the center of the transform !
            itk::ContinuousIndex<double ,3> centerIndex;
            centerIndex[0] = (size[0]-1)/2.0 ;
            centerIndex[1] = (size[1]-1)/2.0 ;
            centerIndex[2] = (size[2]-1)/2.0;

            typename ImageType::PointType center;

            m_Images[i]->TransformContinuousIndexToPhysicalPoint(centerIndex,center);

            m_Transforms[i]->SetImage(m_Images[i]);
            m_Transforms[i]->Initialize();

            m_InverseTransforms[i]->SetImage(m_Images[i]);
            m_InverseTransforms[i]->Initialize();
            m_Transforms[i]->GetInverse(m_InverseTransforms[i]);
            m_BestError[i].resize(m_Images[i]->GetLargestPossibleRegion().GetSize()[2]);

            m_BlurredImages[i] = ImageType::New();
            m_BlurredImages[i] = m_Images[i];


        }

        m_X.set_size(m_NumberOfParameters);
        m_X.fill(0.0);
        m_ScaleX.set_size(m_NumberOfParameters);
        m_ScaleX.fill(1.0);

        //Like all itk Example

        //        const double rotationScale = 0.05;
        //        const double translationScale = 0.001;

        const double rotationScale = 0.5;
        //const double translationScale = 0.5/200.0;

        const double translationScale = 1.0/1000.0;


//        //Euler angles
        m_ScaleX[0] = 1.0;///rotationScale;
        m_ScaleX[1] = 1.0;///rotationScale;
        m_ScaleX[2] = 1.0;///rotationScale;
//        //translation
        m_ScaleX[3] = translationScale;
        m_ScaleX[4] = translationScale;
        m_ScaleX[5] = translationScale;



        // reference :

        m_ReferenceStack = 1;
        m_ReferenceSlice = m_Images[m_ReferenceStack]->GetLargestPossibleRegion().GetSize()[2]/2;

        std::cout<<"Reference image : "<<m_ReferenceStack<<std::endl;
        std::cout<<"Reference slice : "<<m_ReferenceSlice<<std::endl;


    }


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::Update()
{
    srand(time(NULL));

    typename TransformType::ParametersType params;
    params.set_size(m_NumberOfParameters);


    typename TransformType::ParametersType initialParams, simulationParams, zeroParams, MinBounds, MaxBounds;
    typename TransformType::ParametersType ActiveParameters;
    initialParams.set_size(m_NumberOfParameters);
    ActiveParameters.set_size(m_NumberOfParameters);
    zeroParams.set_size(m_NumberOfParameters);
    zeroParams.Fill(0.0);
    initialParams.Fill(0.0);
    simulationParams.Fill(0.0);
    ActiveParameters.Fill(1.0);

   ActiveParameters[5] = 0;



    MinBounds.set_size(m_NumberOfParameters);
    MaxBounds.set_size(m_NumberOfParameters);

    MinBounds[0] = MinBounds[1] = MinBounds[2] = -45; //degrees
    MinBounds[3] = MinBounds[4] = MinBounds[5] = -20;//mm

    MaxBounds[0] = MaxBounds[1] = MaxBounds[2] = 45; //degrees
    MaxBounds[3] = MaxBounds[4] = MaxBounds[5] = 20;//mm



    if(m_VerboseMode)
    {
        std::cout<<"Scales : "<<m_ScaleX<<std::endl;
    }

    std::string Exclusion;

    (m_UseSliceExclusion) ? Exclusion = "ON" : Exclusion = "OFF";
    std::cout<<" * ---------------------------------------- * "<<std::endl;
    std::cout<<" * Registration by intersection of slices * "<<std::endl;
    std::cout<<" * ---------------------------------------- * "<<std::endl<<std::endl;

    std::cout<<" * ----------------! INFOS !----------------- * "<<std::endl;
    std::cout<<"    Number of loops : "<<m_MaxLoop<<std::endl;
    std::cout<<"    Processing Slices Exclusion : "<<Exclusion<<std::endl;
    std::cout<<" * ---------------------------------------- * "<<std::endl<<std::endl;
    std::cout<<"Processing ..."<<std::endl;

    //TODO : Clean this part

    for(unsigned int i = 0; i< m_NumberOfImages; i++)
    {
        m_Transforms[i]->Initialize();
        m_Transforms[i]->GetInverse(m_InverseTransforms[i]);
    }


    for(unsigned loop = 0; loop < m_MaxLoop; loop++) //TODO: add a while loop, and check m_MaxLoop, AND epsilon between two loop (if a parameters change, of sum of errors).
    {
        //****** Multi-Resolution ****
//        this->BlurImages(m_MaxLoop - loop);
//        //************************************
//        // last loop we use original images
//        if(loop == m_MaxLoop -1)
//        {
            for(unsigned int i = 0; i< m_Images.size(); i++)
            {
                m_BlurredImages[i] = m_Images[i];
            }

//        }

//        std::ofstream Rx, Ry, Rz, Tx, Ty, Tz;
//        Rx.open("RxError.txt");
//        Ry.open("RyError.txt");
//        Rz.open("RzError.txt");
//        Tx.open("TxError.txt");
//        Ty.open("TyError.txt");
//        Tz.open("TzError.txt");
        //************************************

        std::cout<<"\rLoop N°: "<<loop+1<<" "<<std::flush;

        for(unsigned int i = 0; i < m_NumberOfImages; i++)
        {
            if(m_VerboseMode)
            {
                std::cout<<"**************** "<<std::endl;
                std::cout<<"Moving Image n° : "<<i<<std::flush;
            }



            typename ImageType::SizeType sizeMov = m_Images[i]->GetLargestPossibleRegion().GetSize();

            unsigned int smov = 0;

            for(smov = 0; smov < sizeMov[2]; smov++ )
            {
                //---------
                typename TransformType::ParametersType p = m_Transforms[i]->GetSliceParameters(smov);
                initialParams[0] = p[0];
                initialParams[1] = p[1];
                initialParams[2] = p[2];
                // Euler 3D
                initialParams[3] = p[3];
                initialParams[4] = p[4];
                initialParams[5] = p[5];


                //
                if(smov == 15 && i == 0)
                {
                    std::cout<<"Initial params : "<<initialParams<<std::endl;
                }
                //---------
                m_CurrentError = 0.0;
                if(m_VerboseMode)
                {
                    std::cout<<"Moving slice n° : "<<smov<<std::endl;
                }


                if(!(smov == m_ReferenceSlice && i == m_ReferenceStack))
                {

                    typename btk::SlicesIntersectionITKCostFunction<ImageType>::Pointer f = btk::SlicesIntersectionITKCostFunction<ImageType>::New();
                    f->SetNumberOfParameters(m_NumberOfParameters);

                    //itk::PowellOptimizer::Pointer optimizer = itk::PowellOptimizer::New();
                    //itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New();// getValue return a array of double
//                    itk::OnePlusOneEvolutionaryOptimizer::Pointer optimizer = itk::OnePlusOneEvolutionaryOptimizer::New();
//                    typedef itk::Statistics::NormalVariateGenerator GeneratorType;
//                    GeneratorType::Pointer generator = GeneratorType::New();
//                    generator->Initialize(12345);
//                    optimizer->MaximizeOff();
//                    optimizer->SetNormalVariateGenerator( generator );
//                    optimizer->Initialize( 10.0 );
//                    optimizer->SetEpsilon( 1.0 );

//                    optimizer->SetMaximumIteration( 200 );

                    //itk::AmoebaOptimizer::Pointer optimizer = itk::AmoebaOptimizer::New();
                    //itk::RegularStepGradientDescentOptimizer::Pointer optimizer = itk::RegularStepGradientDescentOptimizer::New();
                    //itk::ExhaustiveOptimizer::Pointer optimizer =  itk::ExhaustiveOptimizer::New();
                    //itk::ConjugateGradientOptimizer::Pointer optimizer = itk::ConjugateGradientOptimizer::New();
                    //itk::GradientDescentOptimizer::Pointer optimizer = itk::GradientDescentOptimizer::New();
                    //btk::SimulatedAnnealingOptimizer::Pointer optimizer = btk::SimulatedAnnealingOptimizer::New();

                    f->SetVerboseMode(m_VerboseDbg);
                    //f->SetImages(m_Images);
                    f->SetImages(m_BlurredImages);
                    f->SetMasks(m_Masks);
                    f->SetTransforms(m_Transforms);
                    f->SetInverseTransforms(m_InverseTransforms);
                    f->SetMovingImageNum(i);
                    f->SetMovingSliceNum(smov);
                    f->Initialize();


                    //optimizer.minimize(m_X);
                    btk::SmartStepGradientDescentOptimizer::Pointer optimizer = btk::SmartStepGradientDescentOptimizer::New();
                    optimizer->SetCostFunction(f.GetPointer());
                    optimizer->SetNumberOfIterations(300);
                    optimizer->SetMaxStep(1.0);
                    optimizer->SetMinStep(0.01);
                    optimizer->SetMinBounds(MinBounds);
                    optimizer->SetMaxBounds(MaxBounds);
                    optimizer->SetUseBounds(false);
                    optimizer->SetOptimizedParameters(ActiveParameters);


//                    if(smov==15)
//                    {
//                        optimizer->SetVerboseMode(true);
//                    }
//                    else
//                    {
//                        optimizer->SetVerboseMode(false);
//                    }

                    //optimizer->MinimizeOn();
                    //m_Optimizer->MaximizeOn();
//                    optimizer->SetMaximumStepLength( 0.1 );
//                    optimizer->SetMinimumStepLength( 0.001 );
//                    optimizer->SetNumberOfIterations( 200 );
//                    optimizer->SetRelaxationFactor( 0.8 );
                    //optimizer->SetOptimizeWithRestarts(true);
                    //optimizer->SetMaximize(false);
                    //optimizer->SetMinimize(true) ;
                    //optimizer->SetMaximumIteration( 1000 );
                    //optimizer->SetMaximumLineIteration(1000);
                    //optimizer->SetMetricWorstPossibleValue(250.0);

                    // Exchousitve optimiser
//                    itk::ExhaustiveOptimizer::StepsType steps( 6 );
//                    steps[1] = 0.1;
//                    steps[2] = 0.1;
//                    steps[3] = 0.1;
//                    steps[4] = 20;
//                    steps[5] = 20;

//                    steps[6] = 20;
//                    optimizer->SetNumberOfSteps( steps );

                    //Regular step gradient descent
                    //optimizer->MinimizeOn();
//                    optimizer->SetMaximumStepLength( 0.1 );
//                    optimizer->SetMinimumStepLength( 0.01 );
//                    optimizer->SetNumberOfIterations( 1000 );
//                    optimizer->SetRelaxationFactor( 0.8 ); //0.8

                    //Conjugate gradient:
//                    optimizer->SetNumberOfIterations(5000);
//                    optimizer->SetLearningRate( 1.0 );
//                    optimizer->MinimizeOn();
//                    optimizer->MaximizeOff();
                    //optimizer->SetOptimizedParameters(ActiveParameters);

                    optimizer->SetInitialPosition( initialParams );


//                    if(i == 0 && loop == 0)
//                    {
//                        double value = 20.0 ;
//                        double Rmin = btk::MathFunctions::DegreesToRadians(-value);
//                        double Rmax = btk::MathFunctions::DegreesToRadians(value);
//                        simulationParams = this->SimulateMotion(Rmin,Rmax,-value,value);

//                        //optimizer->SetInitialPosition( simulationParams );
//                        optimizer->SetInitialPosition( initialParams );
//                    }
//                    else
//                    {
//                        optimizer->SetInitialPosition( initialParams );
//                        //std::cout<<"Simulation :"<<simulationParams<<std::endl;
//                    }


                    optimizer->SetScales(m_ScaleX);

                    double initialError = f->GetValue(zeroParams);



                    if(f->GetIntersection())
                    {
                        try
                        {
                            optimizer->StartOptimization();
                        }
                        catch(itk::ExceptionObject &obj)
                        {
                            btkCoutMacro("Error : "<<obj);
                        }

                        // if error = max(double) there are no intersection and we do nothing with this slice
                        //m_CurrentError = optimizer->GetCurrentCost();
                        m_CurrentError = optimizer->GetValue();
                        //m_CurrentError = optimizer->GetCurrentValue();

                        if(m_VerboseMode)
                        {
                            std::cout<<"Best error : "<<m_CurrentError<<std::endl;
                        }

                        m_X = optimizer->GetCurrentPosition();
                        params = m_Transforms[i]->GetSliceParameters(smov);

                        if(m_CurrentError == initialError) //  or same as initial
                        {
                            //m_X.fill(0.0);
                            for(unsigned int x = 0; x< params.size(); x++)
                            {
                                m_X[x] = params[x];
                            }
                        }

                        if(m_VerboseMode)
                        {
                            std::cout<<"Final Parameters : "<<m_X<<std::endl;
                        }
                        //***************
//                        if(i==0 && (smov == 16 || smov == 15))
//                        {
//                            std::cout<<"Error : "<<initialError<<std::endl;
//                            std::cout<<"Best error : "<<m_CurrentError<<std::endl;
//                            std::cout<<"Final Parameters : "<<m_X<<std::endl;
//                        }
                        //**************


                        for(unsigned int x = 0; x< params.size(); x++)
                        {

                            params[x] = m_X[x];
                            if(x < 3)
                            {
                                params[x] = MathFunctions::DegreesToRadians(m_X[x]);
                                //params[x] = m_X[x];
                            }
////                            else
////                            {
////                                params[x] = m_X[x];
////                            }
//                            //Centered Euler 3D
//                            else if(x > 5)
//                            {
//                                params[x] = m_X[x];
//                            }

                        }

                        m_Transforms[i]->SetSliceParameters(smov,params);
                        // DEBUG PART
                        //***********************************************
//                        if(i == 0)
//                        {
//                            if(Rx.is_open())
//                            {
//                                Rx<<smov<<" "<<params[0]  <<std::endl;
//                            }
//                            if(Ry.is_open())
//                            {
//                                Ry<<smov<<" "<< params[1] <<std::endl;
//                            }
//                            if(Rz.is_open())
//                            {
//                                Rz<<smov<<" "<<params[2]<<std::endl;
//                            }
//                            // Euler 3D
//                            if(Tx.is_open())
//                            {
//                                Tx<<smov<<" "<<params[3] <<std::endl;
//                            }
//                            if(Ty.is_open())
//                            {
//                                Ty<<smov<<" "<<params[4] <<std::endl;
//                            }
//                            if(Tz.is_open())
//                            {
//                                Tz<<smov<<" "<<params[5]  <<std::endl;
//                            }

//                        }

                        //***********************************************

                    }
                    this->UpdateInfos();
                    m_X.fill(0.0);

                }
                m_BestError[i][smov] = m_CurrentError;

            }
        }
        //DEBUG
//        Rx.close();
//        Ry.close();
//        Rz.close();
//        Tx.close();
//        Ty.close();
//        Tz.close();

        // perform the sum of errors :
        double error = 0.0;
        for(unsigned int i = 0; i< m_BestError.size(); i++)
        {
            for(unsigned j = 0; j< m_BestError[i].size(); j++)
            {
                error +=m_BestError[i][j];
            }
        }

        std::cout<<"**** ERROR : "<<error<<" ****"<<std::endl;

    }

    if(m_UseSliceExclusion)
    {
        //this->SlicesExclusion();
    }



    std::cout<<" Done !"<<std::endl;


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::UpdateInfos()
{
    for(unsigned int i = 0; i<m_NumberOfImages; i++)
    {
        m_Transforms[i]->GetInverse(m_InverseTransforms[i]);
    }
}

//-------------------------------------------------------------------------------------------------
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::SlicesExclusion()
{
    std::vector<double> Mean,Variance, StdDeviation, MedianValues;
    std::vector< std::vector<double> > Median;
    Median.resize(m_BestError.size());


    Mean.resize(m_BestError.size());
    Variance.resize(m_BestError.size());
    StdDeviation.resize(m_BestError.size());
    Median.resize(m_BestError.size());
    MedianValues.resize(m_BestError.size());
    m_Outliers.resize(m_NumberOfImages);
    std::fill(Mean.begin(),Mean.end(), 0.0);
    std::fill(Variance.begin(),Variance.end(), 0.0);
    std::fill(StdDeviation.begin(),StdDeviation.end(), 0.0);
    std::fill(MedianValues.begin(),MedianValues.end(), 0.0);

    typename TransformType::ParametersType initialParams;
    initialParams.set_size(6);
    initialParams.Fill(0);

    std::cout<<"Processing slice exclusion..."<<std::flush;

    for(unsigned int i =0; i< m_BestError.size(); i++)
    {
        Median[i].resize(m_BestError[i].size());
        for(unsigned int j = 0; j< m_BestError[i].size(); j++)
        {
            Mean[i] += m_BestError[i][j];
            Median[i][j] = m_BestError[i][j];
        }

        Mean[i] = Mean[i]/m_BestError[i].size();

    }


    for(unsigned int i =0; i< m_BestError.size(); i++)
    {
        std::sort(Median[i].begin(),Median[i].end());
        MedianValues[i] = Median[i][Median[i].size()/2];
        for(unsigned int j = 0; j< m_BestError[i].size(); j++)
        {
            Variance[i] += (m_BestError[i][j] - Mean[i]) * (m_BestError[i][j] - Mean[i]);
        }

        Variance[i] = Variance[i]/m_BestError[i].size() - 1;
        StdDeviation[i] = std::sqrt(Variance[i]);
    }

    if(m_VerboseMode)
    {
        for(unsigned int im = 0; im< m_Images.size(); im++)
        {
            std::cout<<"Mean Error for Image ["<<im<<"] : "<<Mean[im]<<std::endl;
            std::cout<<"Median Value for Image ["<<im<<"] : "<<MedianValues[im]<<std::endl;
            std::cout<<"Variance for Error of Image ["<<im<<"] : "<<Variance[im]<<std::endl;
            std::cout<<"Standard deviation for Error of Image ["<<im<<"] : "<<StdDeviation[im]<<std::endl;
        }
    }

    for(unsigned int im = 0; im< m_NumberOfImages; im++)
    {
        m_Outliers[im].resize(m_Images[im]->GetLargestPossibleRegion().GetSize()[2]);
        std::fill(m_Outliers[im].begin(),m_Outliers[im].end(),false);

        for(unsigned slice = 0; slice< m_Images[im]->GetLargestPossibleRegion().GetSize()[2]; slice++)
        {
            //if(m_BestError[im][slice] > StdDeviation[im])
            if(m_BestError[im][slice] > 1.25 * MedianValues[im])
            {
                typename TransformType::ParametersType tmpP;
                tmpP = m_Transforms[im]->GetSliceParameters(slice);
                initialParams[3] = tmpP[3];
                initialParams[4] = tmpP[4];
                initialParams[5] = tmpP[5];

                if(m_BestError[im][slice] == DBL_MAX)
                {
                    m_Transforms[im]->SetSliceParameters(slice,initialParams);
                }
                else
                {
                    m_Outliers[im][slice] = true;// Use that for check if we use the slice or not for the super resolution.
                    m_Transforms[im]->SetSliceParameters(slice,initialParams);
                    //std::cout<<" Outlier : ["<<im<<"]["<<slice<<"] "<<std::endl;
                    //std::cout<<"Error : "<<m_BestError[im][slice]<<"Median Error : "<<MedianValues[im]<<std::endl;
                }



            }

        }


    }

    std::cout<<" Done !"<<std::endl;


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
typename MotionCorrectionByIntersection<TImage>::
TransformType::ParametersType MotionCorrectionByIntersection<TImage>::SimulateMotion(double _Rmin, double _Rmax, double _Tmin, double _Tmax)
{
    typename TransformType::ParametersType params;
    params.set_size(6);
    params.Fill(0.0);

    for(unsigned int i = 0; i<6 ;i++)
    {
        if(i < 3)
        {
            params[i] = rand()/((double)RAND_MAX/std::abs(_Rmax - _Rmin)) - std::abs(_Rmin);
        }
        else
        {
            params[i] = rand()/((double)RAND_MAX/std::abs(_Tmax - _Tmin)) - std::abs(_Tmin);
        }

    }
    return params;
}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::
BlurImages(double level)
{
    typedef itk::DiscreteGaussianImageFilter<ImageType, ImageType> gaussianFilter;
    typename gaussianFilter::ArrayType sigma;
    sigma[0] = level;
    sigma[1] = level;
    sigma[2] = 0.0001;
    for(unsigned int i =0; i< m_Images.size(); i++)
    {
        typename gaussianFilter::Pointer blur = gaussianFilter::New();
        blur->SetInput(m_Images[i]);
        blur->SetVariance(sigma);
        blur->Update();
        m_BlurredImages[i] = blur->GetOutput();


    }


}

//-------------------------------------------------------------------------------------------------

}//namespace



#endif // BTKMOTIONCORRECTIONBYINTERSECTION_TXX
