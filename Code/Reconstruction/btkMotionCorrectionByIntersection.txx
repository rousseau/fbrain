#ifndef BTKMOTIONCORRECTIONBYINTERSECTION_TXX
#define BTKMOTIONCORRECTIONBYINTERSECTION_TXX

#include "btkMotionCorrectionByIntersection.h"

#include "itkFRPROptimizer.h"

namespace btk
{
//-------------------------------------------------------------------------------------------------
template<typename TImage>
MotionCorrectionByIntersection<TImage>::MotionCorrectionByIntersection():m_VerboseMode(true),m_MaxLoop(3),m_VerboseDbg(false)
  ,m_CurrentError(0.0),m_UseSliceExclusion(true),m_NumberOfParameters(6)
{
    // Activate VerboseDbg when runngin debug mode
#ifndef NDEBUG
        m_VerboseDbg = true; //true
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
        // initialization
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

            // center of transform
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



        // Set the reference slice. the reference slice never move.
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
    // initialization for the random
    srand(time(NULL));


    typename TransformType::ParametersType params;
    params.set_size(m_NumberOfParameters);

    // Default value for parameters
    typename TransformType::ParametersType initialParams, simulationParams, zeroParams, MinBounds, MaxBounds;
    // Active parameters is usefull for activate/desactivate optimization of some parameters (Tz for example)
    typename TransformType::ParametersType ActiveParameters;
    initialParams.set_size(m_NumberOfParameters);
    ActiveParameters.set_size(m_NumberOfParameters);
    zeroParams.set_size(m_NumberOfParameters);
    zeroParams.Fill(0.0);
    initialParams.Fill(0.0);
    simulationParams.Fill(0.0);
    ActiveParameters.Fill(1.0);

    // Min, Max values of the parameters (for btk::SmartStepGradientDescentOptimizer)
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
    // Terminal Display
    (m_UseSliceExclusion) ? Exclusion = "ON" : Exclusion = "OFF";
    std::cout<<" * ---------------------------------------- * "<<std::endl;
    std::cout<<" * Registration by intersection of slices * "<<std::endl;
    std::cout<<" * ---------------------------------------- * "<<std::endl<<std::endl;

    std::cout<<" * ----------------! INFOS !----------------- * "<<std::endl;
    std::cout<<"    Number of loops : "<<m_MaxLoop<<std::endl;
    std::cout<<"    Processing Slices Exclusion : "<<Exclusion<<std::endl;
    std::cout<<" * ---------------------------------------- * "<<std::endl<<std::endl;
    std::cout<<"Processing ..."<<std::endl;


    for(unsigned int i = 0; i< m_NumberOfImages; i++)
    {
        m_Transforms[i]->Initialize();
        m_Transforms[i]->GetInverse(m_InverseTransforms[i]);
    }

    // Loop over number max of loop
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
          // Uncomment this part if you want to write error per slice in a file (for matlab)
//        std::ofstream Rx, Ry, Rz, Tx, Ty, Tz;
//        Rx.open("RxError.txt");
//        Ry.open("RyError.txt");
//        Rz.open("RzError.txt");
//        Tx.open("TxError.txt");
//        Ty.open("TyError.txt");
//        Tz.open("TzError.txt");
        //************************************

        std::cout<<"\rLoop N°: "<<loop+1<<" "<<std::flush;
        // Loop Over All Images
        for(unsigned int i = 0; i < m_NumberOfImages; i++)
        {
            if(m_VerboseMode)
            {
                std::cout<<"**************** "<<std::endl;
                std::cout<<"Moving Image n° : "<<i<<std::flush;
            }



            typename ImageType::SizeType sizeMov = m_Images[i]->GetLargestPossibleRegion().GetSize();

            unsigned int smov = 0;
            // Loop over slice
            for(smov = 0; smov < sizeMov[2]; smov++ )
            {
                //Initialize parameters with the previous value (identity at the first step)
                typename TransformType::ParametersType p = m_Transforms[i]->GetSliceParameters(smov);
                initialParams[0] = p[0];
                initialParams[1] = p[1];
                initialParams[2] = p[2];
                // Euler 3D
                initialParams[3] = p[3];
                initialParams[4] = p[4];
                initialParams[5] = p[5];

                m_CurrentError = 0.0;// Initialize current error

                if(m_VerboseMode)
                {
                    std::cout<<"Moving slice n° : "<<smov<<std::endl;
                }

                // if smov isn't the reference slice
                if(!(smov == m_ReferenceSlice && i == m_ReferenceStack))
                {
                    // Initialization of the cost function
                    typename btk::SlicesIntersectionITKCostFunction<ImageType>::Pointer f = btk::SlicesIntersectionITKCostFunction<ImageType>::New();
                    f->SetNumberOfParameters(m_NumberOfParameters);
                    f->SetVerboseMode(m_VerboseDbg);
                    //f->SetImages(m_Images);
                    f->SetImages(m_BlurredImages);
                    f->SetMasks(m_Masks);
                    f->SetTransforms(m_Transforms);
                    f->SetInverseTransforms(m_InverseTransforms);
                    f->SetMovingImageNum(i);
                    f->SetMovingSliceNum(smov);
                    f->Initialize();//Don't forget

                    btk::SmartStepGradientDescentOptimizer::Pointer optimizer = btk::SmartStepGradientDescentOptimizer::New();// first optimizer (gradient descent)
                    itk::FRPROptimizer::Pointer optimizer2 = itk::FRPROptimizer::New();//second optimizer to refine result of the first one (this is a test)
                    //TODO: Test with others itk::Optimizers (powell, amoeba...)

                    if(loop <=1) // the first and the second loop we use first optimizer
                    {

                        optimizer->SetCostFunction(f.GetPointer());
                        optimizer->SetNumberOfIterations(300);
                        optimizer->SetMaxStep(1.0);
                        optimizer->SetMinStep(0.01);
                        optimizer->SetMinBounds(MinBounds);
                        optimizer->SetMaxBounds(MaxBounds);
                        optimizer->SetUseBounds(false);
                        optimizer->SetOptimizedParameters(ActiveParameters);
                    }
                    else // otherwise we use second one
                    {
                        optimizer2->SetCostFunction(f.GetPointer());
                        optimizer2->SetStepLength(0.01);
                        optimizer2->SetMaximize(false);
                        optimizer2->SetMaximumIteration( 200 );
                        //optimizer2->SetToFletchReeves();
                        optimizer2->SetToPolakRibiere();

                    }


                    optimizer->SetInitialPosition( initialParams );
                    optimizer2->SetInitialPosition( initialParams );

                    // Before optimization we check if there are intersections
                    double initialError = f->GetValue(zeroParams);

                    // Get if the last evaluation of cost function has an intersection
                    if(f->GetIntersection())
                    {
                        // if yes we can optimize
                        if(loop <=1)//first optimizer
                        {
                            try
                            {
                                optimizer->StartOptimization();
                            }
                            catch(itk::ExceptionObject &obj)
                            {
                                btkCoutMacro("Error : "<<obj);
                            }
                            m_CurrentError = optimizer->GetValue(); // Get The current error
                            m_X = optimizer->GetCurrentPosition();// Get the current parameters
                        }
                        else // second optimizer
                        {
                            try
                            {
                                optimizer2->StartOptimization();
                            }
                            catch(itk::ExceptionObject &obj)
                            {
                                btkCoutMacro("Error : "<<obj);
                            }
                            m_CurrentError = optimizer2->GetValue();// Get the current error
                            m_X = optimizer2->GetCurrentPosition();// Get the current parameters
                        }


                        if(m_VerboseMode)
                        {
                            std::cout<<"Best error : "<<m_CurrentError<<std::endl;
                        }

                        params = m_Transforms[i]->GetSliceParameters(smov);

                        if(m_CurrentError == initialError) // if error is exactly the same as initial error (0 for all parameters)
                        {
                            // We do nothing and we take the initial parameters
                            for(unsigned int x = 0; x< params.size(); x++)
                            {
                                m_X[x] = params[x];
                            }
                        }

                        if(m_VerboseMode)
                        {
                            std::cout<<"Final Parameters : "<<m_X<<std::endl;
                        }

                        //We copy m_X into params and we convert degrees into radians for the sliceBySliceTransform
                        for(unsigned int x = 0; x< params.size(); x++)
                        {

                            params[x] = m_X[x];
                            if(x < 3)
                            {
                                params[x] = MathFunctions::DegreesToRadians(m_X[x]);
                            }


                        }
                        // We set the parameters into the transformation
                        m_Transforms[i]->SetSliceParameters(smov,params);
                        // DEBUG PART: uncomment this if you want to write error in files
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
                    // UpdateInfos compute inverse of m_Transforms after each slice.
                    // if we don't cost function is out of date and results are wrong
                    this->UpdateInfos();
                    m_X.fill(0.0); // reinitialize m_X

                }
                m_BestError[i][smov] = m_CurrentError; // Save the current error on a vector

            }
        }
        //DEBUG: Uncomment this part if you want save errors in files
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
        // after each loop display the sum of errors
        std::cout<<"**** ERROR : "<<error<<" ****"<<std::endl;

    }
    // if you want to use Outliers exclusion
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
        //Compute the inverse
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

    // Compute the Mean Error per images :
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

    // Compute the variance and the standard deviation per images :
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
    // for each image
    for(unsigned int im = 0; im< m_NumberOfImages; im++)
    {
        m_Outliers[im].resize(m_Images[im]->GetLargestPossibleRegion().GetSize()[2]);
        std::fill(m_Outliers[im].begin(),m_Outliers[im].end(),false);
        // for each slices
        for(unsigned slice = 0; slice< m_Images[im]->GetLargestPossibleRegion().GetSize()[2]; slice++)
        {
            //TODO: choose a test for detecting outiers
            //if(m_BestError[im][slice] > StdDeviation[im])
            if(m_BestError[im][slice] > 1.25 * MedianValues[im]) // like K.Kim 2010
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
                    m_Transforms[im]->SetSliceParameters(slice,initialParams); // TODO: reinitialize the parameters of the outlier slice or not ?
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
    // Currently Unused to be removed
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
    //Filter image with a gaussian kernel at each loop.
    // Result with this are not significaly better...
    //To be removed

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
