#ifndef BTKMOTIONCORRECTIONBYINTERSECTION_TXX
#define BTKMOTIONCORRECTIONBYINTERSECTION_TXX

#include "btkMotionCorrectionByIntersection.h"

namespace btk
{
//-------------------------------------------------------------------------------------------------
template<typename TImage>
MotionCorrectionByIntersection<TImage>::MotionCorrectionByIntersection():m_VerboseMode(false),m_MaxLoop(3),m_VerboseDbg(false)
  ,m_CurrentError(0.0),m_UseSliceExclusion(false)
{
#ifndef NDEBUG
        m_VerboseDbg = true;
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
            m_Transforms[i] = Transform::New();
            m_InverseTransforms[i] = Transform::New();
            m_Transforms[i]->SetImage(m_Images[i]);
            m_Transforms[i]->Initialize();
            typename ImageType::SizeType size = m_Images[i]->GetLargestPossibleRegion().GetSize();
            m_Transforms[i]->GetInverse(m_InverseTransforms[i]);
            m_BestError[i].resize(m_Images[i]->GetLargestPossibleRegion().GetSize()[2]);

            m_BlurredImages[i] = ImageType::New();


        }

        m_X.set_size(6);
        m_X.fill(0.0);
        m_ScaleX.set_size(6);
        m_ScaleX.fill(1.0);

        //Like all itk Example

        //        const double rotationScale = 0.05;
        //        const double translationScale = 0.001;

        const double rotationScale = 0.5;
        const double translationScale = 0.5/200.0;


//        //Euler angles
//        m_ScaleX[0] = 1.0/rotationScale;
//        m_ScaleX[1] = 1.0/rotationScale;
//        m_ScaleX[2] = 1.0/rotationScale;
//        //translation
//        m_ScaleX[3] = 1.0/translationScale;
//        m_ScaleX[4] = 1.0/translationScale;
//        m_ScaleX[5] = 1.0/translationScale;

        //translation
        m_ScaleX[3] = translationScale;
        m_ScaleX[4] = translationScale;
        m_ScaleX[5] = translationScale;


        // reference :

        m_ReferenceStack = 1;
        m_ReferenceSlice = m_Images[m_ReferenceStack]->GetLargestPossibleRegion().GetSize()[2]/2;



    this->BlurImages(1);
    }


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::Update()
{
    srand(time(NULL));

    typename Transform::ParametersType params;
    params.set_size(9);

    typename Transform::ParametersType initialParams, simulationParams, zeroParams;
    typename Transform::ParametersType DeltaSimplex;
    initialParams.set_size(6);
    DeltaSimplex.set_size(6);
    zeroParams.set_size(6);
    zeroParams.Fill(0.0);
    initialParams.Fill(0.0);
    simulationParams.Fill(0.0);

    DeltaSimplex[0] = 0.2;
    DeltaSimplex[1] = 0.2;
    DeltaSimplex[2] = 0.5;

    DeltaSimplex[3] = 0.1;
    DeltaSimplex[4] = 0.1;
    DeltaSimplex[5] = 0.1;

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

    for(unsigned loop = 0; loop < m_MaxLoop; loop++)
    {


        //****** Multi-Resolution ****
        this->BlurImages(m_MaxLoop - loop);
        //************************************
        std::ofstream Rx, Ry, Rz, Tx, Ty, Tz;
        Rx.open("RxError.txt");
        Ry.open("RyError.txt");
        Rz.open("RzError.txt");
        Tx.open("TxError.txt");
        Ty.open("TyError.txt");
        Tz.open("TzError.txt");
        //************************************

        std::cout<<"Loop N°: "<<loop+1<<std::flush;

        for(unsigned int i = 0; i < m_NumberOfImages; i++)
        {
            if(m_VerboseMode)
            {

                std::cout<<"**************** "<<std::endl;
                std::cout<<"Loop : "<<loop<<"/r"<<std::flush;
                std::cout<<"Moving Image n° : "<<i<<std::flush;
            }



            typename ImageType::SizeType sizeMov = m_Images[i]->GetLargestPossibleRegion().GetSize();
            //TODO : Parallelize this
            // NOTE :  Parallelization with a call of UpdateInfos may cause segmentation fault
            unsigned int smov = 0;
            //#pragma omp parallel for private(smov) schedule(dynamic)
            for(smov = 0; smov < sizeMov[2]; smov++ )
            {
                //---------
                typename Transform::ParametersType p = m_Transforms[i]->GetSliceParameters(smov);
                initialParams[0] = p[0];
                initialParams[1] = p[1];
                initialParams[2] = p[2];

                initialParams[3] = p[6];
                initialParams[4] = p[7];
                initialParams[5] = p[8];

                //std::cout<<"Initial params : "<<initialParams<<std::endl;
                //---------
                m_CurrentError = 0.0;
                if(m_VerboseMode)
                {
                    std::cout<<"Moving slice n° : "<<smov<<std::endl;
                }


                if(!(smov == m_ReferenceSlice && i == m_ReferenceStack))
                {

                    typename btk::SlicesIntersectionITKCostFunction<ImageType>::Pointer f = btk::SlicesIntersectionITKCostFunction<ImageType>::New();
                    f->SetNumberOfParameters(6);

                    //itk::PowellOptimizer::Pointer optimizer = itk::PowellOptimizer::New();
                    //itk::AmoebaOptimizer::Pointer optimizer = itk::AmoebaOptimizer::New();
                    itk::RegularStepGradientDescentOptimizer::Pointer optimizer = itk::RegularStepGradientDescentOptimizer::New();
                    //itk::ExhaustiveOptimizer::Pointer optimizer =  itk::ExhaustiveOptimizer::New();

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
                    optimizer->SetCostFunction(f.GetPointer());
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

                    optimizer->SetMaximumStepLength( 0.1 );
                    optimizer->SetMinimumStepLength( 0.001 );
                    optimizer->SetNumberOfIterations( 200 );
                    optimizer->SetRelaxationFactor( 0.8 );

                    if(i == 0 && loop == 0)
                    {
                        double value = 5.0 ;
                        double Rmin = btk::MathFunctions::DegreesToRadians(-value);
                        double Rmax = btk::MathFunctions::DegreesToRadians(value);
                        simulationParams = this->SimulateMotion(Rmin,Rmax,-value,value);

                        optimizer->SetInitialPosition( simulationParams );
                    }
                    else
                    {
                        optimizer->SetInitialPosition( initialParams );
                        //std::cout<<"Simulation :"<<simulationParams<<std::endl;
                    }


                    //optimizer->SetMaximumNumberOfIterations(10000);

                    //optimizer->AutomaticInitialSimplexOff();
                    //optimizer->SetInitialSimplexDelta(DeltaSimplex);

                    //optimizer->SetStepLength(1.0);
                    // Taking from ITK's examples
                    //optimizer->SetParametersConvergenceTolerance( 1e-6 ); // about 0.005 degrees
                    //optimizer->SetFunctionConvergenceTolerance( 1e-4 );  // variation in metric value

                    //optimizer->SetMaximumNumberOfIterations( 200 );

                    optimizer->SetScales(m_ScaleX);

                    double initialError = f->GetValue(zeroParams);

                    if(initialError != DBL_MAX)
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
                        //m_CurrentError = optimizer->GetValue();
                        //m_CurrentError = optimizer->GetCurrentValue();

                        if(m_VerboseMode)
                        {
                            std::cout<<"Best error : "<<m_CurrentError<<std::endl;
                        }

                        m_X = optimizer->GetCurrentPosition();

                        if(m_CurrentError == DBL_MAX || m_CurrentError == initialError)
                        {
                            m_X.fill(0.0);
                        }


                        if(m_VerboseMode)
                        {
                            std::cout<<"Final Parameters : "<<m_X<<std::endl;
                        }

                        params = m_Transforms[i]->GetSliceParameters(smov);


                        for(unsigned int x = 0; x< params.size(); x++)
                        {
                            if(x < 3)
                            {
                                //params[x] = MathFunctions::DegreesToRadians(m_X[x]);
                                params[x] = m_X[x];
                            }
                            else if(x > 5)
                            {
                                params[x] = m_X[x-3];
                            }

                        }

                        m_Transforms[i]->SetSliceParameters(smov,params);
                        //***********************************************
                        if(i == 0)
                        {
                            if(Rx.is_open())
                            {
                                Rx<<smov<<" "<<params[0] <<std::endl;
                            }
                            if(Ry.is_open())
                            {
                                Ry<<smov<<" "<< params[1]<<std::endl;
                            }
                            if(Rz.is_open())
                            {
                                Rz<<smov<<" "<<params[2]<<std::endl;
                            }
                            if(Tx.is_open())
                            {
                                Tx<<smov<<" "<<params[6]<<std::endl;
                            }
                            if(Ty.is_open())
                            {
                                Ty<<smov<<" "<<params[7]  <<std::endl;
                            }
                            if(Tz.is_open())
                            {
                                Tz<<smov<<" "<<params[8]  <<std::endl;
                            }
                        }

                        //***********************************************

                    }
                    this->UpdateInfos();
                    m_X.fill(0.0);

                }
                m_BestError[i][smov] = m_CurrentError;


            }
        }
        Rx.close();
        Ry.close();
        Rz.close();
        Tx.close();
        Ty.close();
        Tz.close();
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
    std::fill(Mean.begin(),Mean.end(), 0.0);
    std::fill(Variance.begin(),Variance.end(), 0.0);
    std::fill(StdDeviation.begin(),StdDeviation.end(), 0.0);
    std::fill(MedianValues.begin(),MedianValues.end(), 0.0);

    typename Transform::ParametersType initialParams;
    initialParams.set_size(9);
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
        for(unsigned slice = 0; slice< m_Images[im]->GetLargestPossibleRegion().GetSize()[2]; slice++)
        {
            //if(m_BestError[im][slice] > StdDeviation[im])
            if(m_BestError[im][slice] > 1.25 * MedianValues[im])
            {
                typename Transform::ParametersType tmpP;
                tmpP = m_Transforms[im]->GetSliceParameters(slice);
                initialParams[3] = tmpP[3];
                initialParams[4] = tmpP[4];
                initialParams[5] = tmpP[5];

                m_Transforms[im]->SetSliceParameters(slice,initialParams);
            }

        }
    }

    std::cout<<" Done !"<<std::endl;


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
typename MotionCorrectionByIntersection<TImage>::
Transform::ParametersType MotionCorrectionByIntersection<TImage>::SimulateMotion(double _Rmin, double _Rmax, double _Tmin, double _Tmax)
{
    typename Transform::ParametersType params;
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

}







#endif // BTKMOTIONCORRECTIONBYINTERSECTION_TXX
