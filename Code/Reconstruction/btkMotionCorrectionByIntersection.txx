#ifndef BTKMOTIONCORRECTIONBYINTERSECTION_TXX
#define BTKMOTIONCORRECTIONBYINTERSECTION_TXX

#include "btkMotionCorrectionByIntersection.h"

namespace btk
{
//-------------------------------------------------------------------------------------------------
template<typename TImage>
MotionCorrectionByIntersection<TImage>::MotionCorrectionByIntersection():m_VerboseMode(false),m_MaxLoop(5),m_VerboseDbg(false)
  ,m_CurrentError(0.0),m_UseSliceExclusion(true)
{

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


        for(int i = 0; i< m_NumberOfImages; i++)
        {
            m_Transforms[i] = Transform::New();
            m_InverseTransforms[i] = Transform::New();
            m_Transforms[i]->SetImage(m_Images[i]);
            m_Transforms[i]->Initialize();
            typename ImageType::SizeType size = m_Images[i]->GetLargestPossibleRegion().GetSize();
            m_Transforms[i]->GetInverse(m_InverseTransforms[i]);
            m_BestError[i].resize(m_Images[i]->GetLargestPossibleRegion().GetSize()[2]);


        }

        m_X.set_size(6);
        m_X.fill(0);
        m_ScaleX.set_size(6);
        m_ScaleX.fill(1.0);

        //Like all itk Example

//        const double rotationScale = 0.05;
//        const double translationScale = 0.001;

        const double rotationScale = 180.0;
        const double translationScale = 100.0;


//        //Euler angles
//        m_ScaleX[0] = 1.0/rotationScale;
//        m_ScaleX[1] = 1.0/rotationScale;
//        m_ScaleX[2] = 1.0/rotationScale;
//        //translation
//        m_ScaleX[3] = 1.0/translationScale;
//        m_ScaleX[4] = 1.0/translationScale;
//        m_ScaleX[5] = 1.0/translationScale;

        //Euler angles
        m_ScaleX[0] = rotationScale;
        m_ScaleX[1] = rotationScale;
        m_ScaleX[2] = rotationScale;
        //translation
        m_ScaleX[3] = translationScale;
        m_ScaleX[4] = translationScale;
        m_ScaleX[5] = translationScale;

        // reference :

        m_ReferenceStack = 0;
        m_ReferenceSlice = m_Images[m_ReferenceStack]->GetLargestPossibleRegion().GetSize()[2]/2;


        #ifndef NDEBUG
        m_VerboseDbg = false;
        #else
        m_VerboseDbg = false;
        #endif

    }


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void MotionCorrectionByIntersection<TImage>::Update()
{

    typename Transform::ParametersType params;
    params.set_size(6);

    typename Transform::ParametersType initialParams;
    typename Transform::ParametersType DeltaSimplex;
    initialParams.set_size(6);
    DeltaSimplex.set_size(6);
    initialParams.Fill(0);
    //    DeltaSimplex[0] = 10.0/180.0;
    //    DeltaSimplex[1] = 10.0/180.0;
    //    DeltaSimplex[2] = 30.0/180.0;

    //    DeltaSimplex[3] = 10.0/100.0;
    //    DeltaSimplex[4] = 10.0/100.0;
    //    DeltaSimplex[5] = 0.0/100.0;

    std::cout<<"Scales : "<<m_ScaleX<<std::endl;

    std::cout<<" * Registration by intersection of slices * "<<std::endl;
    std::cout<<"Processing ..."<<std::endl;
    for(unsigned loop = 0; loop < 1/*m_MaxLoop*/; loop++)
    {

        std::cout<<"Loop N°: "<<loop+1<<std::endl;

        for(unsigned int i = 0; i < m_NumberOfImages; i++)
        {
            if(m_VerboseMode)
            {

                std::cout<<"**************** "<<std::endl;
                std::cout<<"Loop : "<<loop<<std::endl;
                std::cout<<"Moving Image n° : "<<i<<std::endl;
            }



            typename ImageType::SizeType sizeMov = m_Images[i]->GetLargestPossibleRegion().GetSize();
            //TODO : Parallelize this
            for(unsigned int smov = 0; smov < sizeMov[2]; smov++ )
            {
                m_CurrentError = 0.0;
                if(m_VerboseMode)
                {
                    std::cout<<"Moving slice n° : "<<smov<<std::endl;
                }


                if(!(smov == m_ReferenceSlice && i == m_ReferenceStack))
                {

                    typename btk::SlicesIntersectionITKCostFunction<ImageType>::Pointer f = btk::SlicesIntersectionITKCostFunction<ImageType>::New();
                    f->SetNumberOfParameters(6);

                    //itk::PowellOptimizer::Pointer powell = itk::PowellOptimizer::New();
                    itk::AmoebaOptimizer::Pointer powell = itk::AmoebaOptimizer::New();

                    f->SetVerboseMode(m_VerboseDbg);
                    f->SetImages(m_Images);
                    f->SetMasks(m_Masks);
                    f->SetTransforms(m_Transforms);
                    f->SetInverseTransforms(m_InverseTransforms);
                    f->SetMovingImageNum(i);
                    f->SetMovingSliceNum(smov);

                    f->Initialize();


                    //powell.minimize(m_X);
                    powell->SetCostFunction(f.GetPointer());
                    powell->SetOptimizeWithRestarts(false);
                    powell->SetMaximize(false);
                    //powell->SetMaximumIteration( 1000 );
                    //powell->SetMaximumLineIteration(1000);
                    //powell->SetMetricWorstPossibleValue(250.0);
//                    std::cout<<powell->GetValueTolerance()<<std::endl;
//                    std::cout<<powell->GetStepLength()<<std::endl;

                    powell->SetInitialPosition( initialParams );

                    powell->SetMaximumNumberOfIterations(10000);

                    //powell->AutomaticInitialSimplexOff();
                    //powell->SetInitialSimplexDelta(DeltaSimplex);
                    //std::cout<<"Delta Simplex : "<<powell->GetInitialSimplexDelta()<<std::endl;

                    powell->SetScales(m_ScaleX);
                    double initialError = f->GetValue(initialParams);

                    if(initialError != DBL_MAX)
                    {
                        try
                        {
                            powell->StartOptimization();
                        }
                        catch(itk::ExceptionObject &obj)
                        {
                            btkCoutMacro("Error : "<<obj);
                        }

                        // if error = max(double) there are no intersection and we do nothing with this slice
                        //m_CurrentError = powell->GetCurrentCost();
                        m_CurrentError = powell->GetValue();

                        if(m_VerboseMode)
                        {
                            std::cout<<"Best error : "<<m_CurrentError<<std::endl;
                        }

                        m_X = powell->GetCurrentPosition();

                        if(m_CurrentError == DBL_MAX || m_CurrentError == initialError)
                        {
                            m_X.fill(0);
                        }


                        if(m_VerboseMode)
                        {
                            std::cout<<"Final Parameters : "<<m_X<<std::endl;
                        }


                        for(unsigned int x = 0; x< m_X.size(); x++)
                        {
                            if(x <3)
                            {
                              params[x] = MathFunctions::DegreesToRadians(m_X[x]);
                            }
                            else
                            {
                                params[x] = m_X[x];
                            }

                        }

                        m_Transforms[i]->SetSliceParameters(smov,params);

                    }
                    this->UpdateInfos();
                    m_X.fill(0);

                }
                m_BestError[i][smov] = m_CurrentError;
            }
        }
    }

    if(m_UseSliceExclusion)
    {
        this->SlicesExclusion();
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
    initialParams.set_size(6);
    initialParams.Fill(0);

    std::cout<<"Processing slice exclusion..."<<std::endl;

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
                  m_Transforms[im]->SetSliceParameters(slice,initialParams);
            }

        }
    }


}
//-------------------------------------------------------------------------------------------------



}







#endif // BTKMOTIONCORRECTIONBYINTERSECTION_TXX
