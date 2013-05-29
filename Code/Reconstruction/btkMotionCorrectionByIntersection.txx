#ifndef BTKMOTIONCORRECTIONBYINTERSECTION_TXX
#define BTKMOTIONCORRECTIONBYINTERSECTION_TXX

#include "btkMotionCorrectionByIntersection.h"
#include "itkFRPROptimizer.h"
#include "itkConjugateGradientOptimizer.h"


namespace btk
{
//-------------------------------------------------------------------------------------------------
template<typename TImage>
MotionCorrectionByIntersection<TImage>::MotionCorrectionByIntersection():m_VerboseMode(true),m_MaxLoop(3),m_VerboseDbg(false)
  ,m_CurrentError(0.0),m_UseSliceExclusion(true),m_NumberOfParameters(6)
{
    // Activate VerboseDbg when running debug mode
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
        m_ReferenceStack = 0;
        m_ReferenceSlice = m_Images[m_ReferenceStack]->GetLargestPossibleRegion().GetSize()[2]/2;

        std::cout<<"Reference image : "<<m_ReferenceStack<<std::endl;
        std::cout<<"Reference slice : "<<m_ReferenceSlice<<std::endl;


    }


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
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::Update()
{

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

    MinBounds[0] = MinBounds[1] = MinBounds[2] = -50; //degrees
    MinBounds[3] = MinBounds[4] = MinBounds[5] = -50;//mm

    MaxBounds[0] = MaxBounds[1] = MaxBounds[2] = 50; //degrees
    MaxBounds[3] = MaxBounds[4] = MaxBounds[5] = 50;//mm



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



    unsigned int loop = 0;

    // Loop Over All Images


    unsigned int nbGroup = 1;
    bool stopGrouping = false;
    bool incrementLoop=false;
    while(!stopGrouping)
    {
        std::cout<<"\rLoop N°: "<<loop+1<<" "<<std::flush;
        ActiveParameters.Fill(0.0);

        for(unsigned int i = 0; i < m_NumberOfImages; i++)
        {
            std::cout<<"Image "<<i<<std::endl;

            if(m_VerboseMode)
            {
                std::cout<<"**************** "<<std::endl;
                std::cout<<"Moving Image n° : "<<i<<std::endl;
            }



            typename ImageType::SizeType sizeMov = m_Images[i]->GetLargestPossibleRegion().GetSize();
            std::vector< unsigned int > SlicesGroup;
            SlicesGroup.resize(sizeMov[2]);



            if(nbGroup >= sizeMov[2])
            {
                nbGroup = sizeMov[2];
                incrementLoop = true;
                std::cout<<"start looping"<<std::endl;
                std::cout<<"NB Group * : "<<nbGroup<<std::endl;
            }
            else
            {
                incrementLoop = false;
                std::cout<<"NB Group : "<<nbGroup<<std::endl;
            }

            if(loop >= m_MaxLoop)
            {
                stopGrouping = true;
            }



            unsigned int smov = 0;
            unsigned int referenceGroup = 1000; //avoid random initialization
            // Loop over slice
            for(smov = 0; smov < sizeMov[2]; smov++ )
            {
                SlicesGroup[smov] = smov%nbGroup;
                if(i == m_ReferenceStack && smov == m_ReferenceSlice)
                {
                    //check the group of the reference slice
                    referenceGroup = smov%nbGroup;
                    std::cout<<"reference slice  is "<<smov<<" in image : "<<i<<" in group : "<<smov%nbGroup<<std::endl;
                }

                //std::cout<<"slice num : "<<smov<<"going to group : "<<smov%nbGroup<<std::endl;

            }
            for(unsigned int g = 0; g< nbGroup; g++)
            {
                //Initialize parameters with the previous value (identity at the first step)
                for(unsigned int slice = 0; slice< sizeMov[2]; slice++)
                {
                    if(SlicesGroup[slice] == g)
                    {

                        typename TransformType::ParametersType p = m_Transforms[i]->GetSliceParameters(slice);
                        initialParams[0] = btk::MathFunctions::RadiansToDegrees(p[0]);
                        initialParams[1] = btk::MathFunctions::RadiansToDegrees(p[1]);
                        initialParams[2] = btk::MathFunctions::RadiansToDegrees(p[2]);
                        // Euler 3D
                        initialParams[3] = p[3];
                        initialParams[4] = p[4];
                        initialParams[5] = p[5];
                        break;
                    }

                }


                typename btk::SlicesIntersectionITKCostFunction<ImageType>::Pointer f = btk::SlicesIntersectionITKCostFunction<ImageType>::New();
                f->SetNumberOfParameters(m_NumberOfParameters);
                f->SetVerboseMode(m_VerboseDbg);
                f->SetImages(m_Images);
                //f->SetImages(m_BlurredImages);
                f->SetMasks(m_Masks);
                f->SetTransforms(m_Transforms);
                f->SetInverseTransforms(m_InverseTransforms);
                f->SetMovingImageNum(i);
                //f->SetMovingSliceNum(smov);
                f->SetGroupNum(g);
                f->SetSlicesGroup(SlicesGroup);
                f->Initialize();//Don't forget

                //itk::PowellOptimizer::Pointer optimizer = itk::PowellOptimizer::New();
                //itk::ConjugateGradientOptimizer::Pointer optimizer = itk::ConjugateGradientOptimizer::New();
               //itk::RegularStepGradientDescentOptimizer::Pointer optimizer = itk::RegularStepGradientDescentOptimizer::New();
//                itk::GradientDescentOptimizer::Pointer optimizer = itk::GradientDescentOptimizer::New();
//                itk::FRPROptimizer::Pointer optimizer = itk::FRPROptimizer::New();
//                optimizer->SetToPolakRibiere();
//                optimizer->SetMaximumStepLength (0.8);
//                optimizer->SetMinimumStepLength (0.1);
                //optimizer->SetCostFunction(f.GetPointer());
                //optimizer->SetMaximumIteration(500);
//                optimizer->SetLearningRate(0.1);

                btk::SmartStepGradientDescentOptimizer::Pointer optimizer = btk::SmartStepGradientDescentOptimizer::New();// first optimizer (gradient descent)
                optimizer->SetCostFunction(f.GetPointer());
                optimizer->SetNumberOfIterations(1000);
                optimizer->SetMaxStep(5.0);
                optimizer->SetMinStep(0.05);
                optimizer->SetMinBounds(MinBounds);
                optimizer->SetMaxBounds(MaxBounds);
                optimizer->SetUseBounds(false);

                //optimizer->SetOptimizedParameters(ActiveParameters);
                m_CurrentError = 0.0;// Initialize current error
                optimizer->SetVerboseMode(false);



//DEBUG
//                if(i == 1 && g == 30 && loop == 1)
//                {
//                    std::cout<<i<<std::endl;
//                    std::cout<<g<<std::endl;
//                    std::cout<<loop<<std::endl;
//                    std::cout<<incrementLoop<<std::endl;

//                    std::cout<<"Optimizer set to a random value "<<std::endl;
//                    initialParams[5] = 5.0;
//                    //initialParams.Fill(10.0);
//                }
                optimizer->SetInitialPosition( initialParams );


                // Before optimization we check if there are intersections

                double initialError =  f->GetValue(initialParams);


                bool DoOptimization = true;

                if(i == m_ReferenceStack && g == referenceGroup)
                {
                    DoOptimization = false;
                    std::cout<<"No Optimization this time "<<std::endl;
                }


                // Get if the last evaluation of cost function has an intersection
                if(f->GetIntersection() && DoOptimization)
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


                    if(m_VerboseMode)
                    {
                        std::cout<<"slice : "<<g<<std::endl;
                        std::cout<<"initial Error : "<<initialError<<std::endl;
                        std::cout<<"Best error : "<<m_CurrentError<<std::endl;
                        std::cout<<"Parameters : "<<m_X<<std::endl;
                    }


                    //Save previous parameter values
                    //typename TransformType::ParametersType previousP = m_Transforms[i]->GetSliceParameters(smov);

                    //Set new values to evaluate the overall cost function
                    params[0] = MathFunctions::DegreesToRadians(m_X[0]);
                    params[1] = MathFunctions::DegreesToRadians(m_X[1]);
                    params[2] = MathFunctions::DegreesToRadians(m_X[2]);
                    params[3] = m_X[3];
                    params[4] = m_X[4];
                    params[5] = m_X[5];

                    for(unsigned int slice = 0; slice< sizeMov[2]; slice++)
                    {
                        if(SlicesGroup[slice] == g)
                        {
                            m_Transforms[i]->SetSliceParameters(slice,params);
                            m_BestError[i][slice] = m_CurrentError;
                        }

                    }

                }
                // UpdateInfos compute inverse of m_Transforms after each slice.
                // if we don't cost function is out of date and results are wrong
                this->UpdateInfos();
                m_X.fill(0.0); // reinitialize m_X
            }

        }

        nbGroup*=2;





        // perform the sum of errors :
        double error = 0.0;
        for(unsigned int i = 0; i< m_BestError.size(); i++)
        {
            for(unsigned j = 0; j< m_BestError[i].size(); j++)
            {
                error +=m_BestError[i][j];
            }
        }

        std::cout<<"Cumulated Error : "<<error<<std::endl;



        if(incrementLoop)
        {
            loop++;
        }

    }






    std::cout<<" Done !"<<std::endl;
}




}//namespace



#endif // BTKMOTIONCORRECTIONBYINTERSECTION_TXX
