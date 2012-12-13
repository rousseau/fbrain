#ifndef BTKSlicesIntersectionVNLCostFunction_TXX
#define BTKSlicesIntersectionVNLCostFunction_TXX

#include "btkSlicesIntersectionVNLCostFunction.hxx"
namespace btk
{

//-------------------------------------------------------------------------------------------------
template<class TImage>
SlicesIntersectionVNLCostFunction<TImage>::SlicesIntersectionVNLCostFunction(unsigned int dim)
    :vnl_cost_function(dim),m_NumberOfPointsOfLine(1000),m_VerboseMode(false),m_MovingImageNum(0),m_MovingSliceNum(0)
{
}
//-------------------------------------------------------------------------------------------------
template<class TImage>
SlicesIntersectionVNLCostFunction<TImage>::~SlicesIntersectionVNLCostFunction()
{
}
//-------------------------------------------------------------------------------------------------
template<class TImage>
void SlicesIntersectionVNLCostFunction<TImage>::Initialize()
{
    if(m_Images.empty() || m_Masks.empty() || m_Transforms.empty() || m_InverseTransforms.empty())
    {
        //TODO: Add More details (m_MovingSlice, ...)
        btkException("Initialization error ! Some Needed arguments are missing (check Images, masks, transforms or movingslice");

    }
    else
    {
        m_NumberOfImages = m_Images.size();
        m_Interpolators.resize(m_NumberOfImages);
        for(unsigned int i = 0; i<m_NumberOfImages; i++)
        {
            m_Interpolators[i] = Interpolator::New();
            m_Interpolators[i]->SetInputImage(m_Images[i]);
        }

        //----
        m_Transform = TransformType::New();
        //---
        m_Transform->SetIdentity();
        typename ImageType::SizeType size;
        size =  m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize();
        m_CenterOfMovingSlice.resize(3);

        // For the center of the transform !
        ContinuousIndexType centerIndex;
        centerIndex[0] = (size[0]-1)/2.0 ;
        centerIndex[1] = (size[1]-1)/2.0 ;
        centerIndex[2] = m_MovingSliceNum ;

        typename ImageType::PointType center;

        m_Images[m_MovingImageNum]->TransformContinuousIndexToPhysicalPoint(centerIndex,center);

        m_Transform->SetCenter(center);

    }


}
//-------------------------------------------------------------------------------------------------
template<class TImage>
double SlicesIntersectionVNLCostFunction<TImage>::f(const vnl_vector<double> &x) const
{

    double CostFunction = 0;
    unsigned int NumberOfIntersectedVoxels = 0;
    double SumOfIntersectedVoxels = 0;



    TransformType::ParametersType params;
   // params.SetSize(x.size());
    //------
    params.SetSize(9);
    //------
    //params.Fill(0.0);


    params = m_Transform->GetParameters();
    for(unsigned int i = 0; i< params.size(); i++)
    {
        if(i < 3)
        {
           // params[i] = MathFunctions::DegreesToRadians(x[i]);
            params[i] = x[i];
        }
        else if(i > 5)
        {
            params[i] = x[i-3];
        }
        else
        {
            params[i] = params[i];
        }

    }

    if(m_VerboseMode)
    {
        std::cout<<"Parameters : "<<params<<std::endl;
    }

    m_Transform->SetParameters(params);
    //m_Transform->SetSliceParameters(m_MovingSliceNum,params); // Parameters found by otpimizers and to apply to the current moving slice


    TransformType::Pointer InverseX = TransformType::New();
    InverseX->SetIdentity();
    InverseX->SetCenter(m_Transform->GetCenter());
    InverseX->SetFixedParameters(m_Transform->GetFixedParameters());

    m_Transform->GetInverse(InverseX);


    //TODO: paralelize this...
    // 39min to 27min when parallalized
    unsigned int ifixed = 0;

    for( ifixed = 0; ifixed< m_Images.size(); ifixed ++)
    {
        if(ifixed != m_MovingImageNum)
        {
            unsigned int numberOfFixedSlices = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[2];

            // TODO : Check if it is possible to paralelize this :
            unsigned int sfixed = 0;
            #pragma omp parallel for private(sfixed) schedule(dynamic)
            for(sfixed = 0; sfixed < numberOfFixedSlices; sfixed++)
            {


                //      slice
                //     1------2
                //     |      |
                //     |      |
                //     4------3

                // 4 corner of moving slice:
                typename ImageType::IndexType MovingCorner1, MovingCorner2, MovingCorner3, MovingCorner4;

                // 4 corners of fixed slice :
                typename ImageType::IndexType FixedCorner1, FixedCorner2, FixedCorner3, FixedCorner4;

                typename ImageType::RegionType MovingRegion, FixedRegion;
                typename ImageType::SizeType MovingRgSize, FixedRgSize;


                MovingRgSize = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize();
                MovingRgSize[2] = 1;

                FixedRgSize = m_Images[ifixed]->GetLargestPossibleRegion().GetSize();
                FixedRgSize[2] = 1;


                MovingCorner1[0] = 0; MovingCorner1[1] = 0; MovingCorner1[2] = m_MovingSliceNum;
                MovingCorner2[0] = 0; MovingCorner2[1] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[1] - 1; MovingCorner2[2] = m_MovingSliceNum;
                MovingCorner3[0] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[0] - 1; MovingCorner3[1] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[1] - 1;MovingCorner3[2] = m_MovingSliceNum;
                MovingCorner4[0] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[0] - 1;MovingCorner4[1] = 0; MovingCorner4[2] = m_MovingSliceNum;

                FixedCorner1[0] = 0; FixedCorner1[1] = 0; FixedCorner1[2] = sfixed;
                FixedCorner2[0] = 0; FixedCorner2[1] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[1] - 1; FixedCorner2[2] = sfixed;
                FixedCorner3[0] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[0] - 1; FixedCorner3[1] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[1] - 1;FixedCorner3[2] = sfixed;
                FixedCorner4[0] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[0] - 1;FixedCorner4[1] = 0; FixedCorner4[2] = sfixed;

                MovingRegion.SetIndex(MovingCorner1);
                FixedRegion.SetIndex(FixedCorner1);

                MovingRegion.SetSize(MovingRgSize);
                FixedRegion.SetSize(FixedRgSize);




                typename ImageType::PointType Point1, Point2, CurrentPoint, TCurrentPoint, CorrespondingPoint;
                bool FoundP1, FoundP2, intersection;

                intersection = FoundP1 = FoundP2 = false;

                typename ImageType::IndexType  CurrentIndex;
                typename Interpolator::ContinuousIndexType CorrespondingIndex;

                typedef itk::LineConstIterator<ImageType> LineIterator;

                std::vector<LineIterator> IteratorTab;


                LineIterator itM12(m_Images[m_MovingImageNum], MovingCorner1, MovingCorner2);
                LineIterator itM23(m_Images[m_MovingImageNum], MovingCorner2, MovingCorner3);
                LineIterator itM34(m_Images[m_MovingImageNum], MovingCorner3, MovingCorner4);
                LineIterator itM41(m_Images[m_MovingImageNum], MovingCorner4, MovingCorner1);


                // for the fixed slice we take the opposite
                LineIterator itF12(m_Images[ifixed], FixedCorner1, FixedCorner4);
                LineIterator itF23(m_Images[ifixed], FixedCorner4, FixedCorner3);
                LineIterator itF34(m_Images[ifixed], FixedCorner3, FixedCorner2);
                LineIterator itF41(m_Images[ifixed], FixedCorner2, FixedCorner1);

                IteratorTab.push_back(itM12);
                IteratorTab.push_back(itM23);
                IteratorTab.push_back(itM34);
                IteratorTab.push_back(itM41);
                IteratorTab.push_back(itF12);
                IteratorTab.push_back(itF23);
                IteratorTab.push_back(itF34);
                IteratorTab.push_back(itF41);


                std::vector<typename ImageType::PointType> P1, P2;
                for(unsigned int i = 0; i< IteratorTab.size(); i++)
                {
                    LineIterator it = IteratorTab[i];
                    bool ConsecutivePoint = true;
                    unsigned int countPoint = 0;
                   for(it.GoToBegin(); !it.IsAtEnd(); ++it)
                    {

                        CurrentIndex = it.GetIndex();
                        // i 0-3 for moving slice, and i 4-7 for fixed slice
                        if(i<4)// Moving Slice
                        {
                            m_Images[m_MovingImageNum]->TransformIndexToPhysicalPoint(CurrentIndex, CurrentPoint);
                            // Apply X
                            TCurrentPoint = m_Transform->TransformPoint(CurrentPoint);
                            // Apply Inverse of the Optimized founded transform (if not identity)
                            CorrespondingPoint = m_InverseTransforms[ifixed]->TransformPoint(TCurrentPoint);
                            // Found the Corresponding Index in the fixed slice !
                            m_Images[ifixed]->TransformPhysicalPointToContinuousIndex(CorrespondingPoint,CorrespondingIndex );

                            // Is in the fixed slice ?
                            if(FixedRegion.IsInside(CorrespondingIndex))
                            {
                                countPoint++;

                                ConsecutivePoint = true;
                                if(FoundP1 == false && FoundP2 == false)
                                {
                                    if(countPoint > 5)
                                    {
                                        P1.clear();
                                        break;
                                    }

                                    P1.push_back(TCurrentPoint);


                                }
                                else if(FoundP1 == true && FoundP2 == false)
                                {
                                    if(countPoint > 5)
                                    {
                                        P2.clear();
                                        break;
                                    }
                                    P2.push_back(TCurrentPoint);

                                }
                                else if(FoundP1 && FoundP2)
                                {
                                    // You Should not be here !
                                    break;
                                }
                            }
                            else
                            {
                                ConsecutivePoint = false;
                            }

                            if(!ConsecutivePoint) // When we have all candidate for P1 or P2
                            {
                                if(P1.size() > 0) // if we have found some candidates for P1
                                {
                                    // check the candidate who have the Z coordinate closest to the moving slice num
                                    double Z = (double)sfixed;
                                    double min = 99999.9;
                                    typename Interpolator::ContinuousIndexType ix;
                                    for(unsigned int j = 0; j< P1.size(); j++)
                                    {
                                        m_Images[ifixed]->TransformPhysicalPointToContinuousIndex(P1[j],ix);

                                        if(fabs(ix[2]-Z)< min)
                                        {
                                            min = fabs(ix[2]-Z);

                                            Point1[0] = P1[j][0];
                                            Point1[1] = P1[j][1];
                                            Point1[2] = P1[j][2];
                                        }
                                    }
                                    FoundP1 = true;

                                    P1.clear();
                                    break;
                                }
                                else if(FoundP1 && P2.size() > 0)
                                {
                                    // check the candidate who have the Z coordinate closest to the moving slice num
                                    double Z = (double)sfixed;
                                    double min = 99999.9;
                                    typename Interpolator::ContinuousIndexType ix;
                                    for(unsigned int j = 0; j< P2.size(); j++)
                                    {
                                        m_Images[ifixed]->TransformPhysicalPointToContinuousIndex(P2[j],ix);
                                        if(fabs(ix[2]-Z) < min)
                                        {
                                            min = fabs(ix[2]-Z);
                                            Point2[0] = P2[j][0];
                                            Point2[1] = P2[j][1];
                                            Point2[2] = P2[j][2];
                                        }
                                    }
                                    FoundP2 = true;

                                    P2.clear();
                                    break;
                                }
                            }

                        } //moving slice

                        //TODO : Duplication of code...
                        else // fixed slice
                        {
                            m_Images[ifixed]->TransformIndexToPhysicalPoint(CurrentIndex, CurrentPoint);
                            //Apply  of the Optimized founded transform (if not identity)
                            TCurrentPoint = m_Transforms[ifixed]->TransformPoint(CurrentPoint);
                            //Apply inverse of X for return in the moving slice space
                            CorrespondingPoint = InverseX->TransformPoint(TCurrentPoint);

                            m_Images[m_MovingImageNum]->TransformPhysicalPointToContinuousIndex(CorrespondingPoint, CorrespondingIndex);
                            if(MovingRegion.IsInside(CorrespondingIndex))
                            {
                                countPoint++;

                                ConsecutivePoint = true;
                                if(FoundP1 == false && FoundP2 == false)
                                {
                                    if(countPoint > 5)
                                    {
                                        P1.clear();
                                        break;
                                    }
                                    P1.push_back(TCurrentPoint);

                                }
                                else if(FoundP1 == true && FoundP2 == false)
                                {
                                    if(countPoint > 5)
                                    {
                                        P2.clear();
                                        break;
                                    }
                                    P2.push_back(TCurrentPoint);

                                }
                                else if(FoundP1 && FoundP2)
                                {
                                    // You Should not be here !
                                    break;
                                }
                            }
                            else
                            {
                                ConsecutivePoint = false;
                            }

                            if(!ConsecutivePoint) // When we have all candidate for P1 or P2
                            {
                                if(P1.size() > 0) // if we have found some candidates for P1
                                {
                                    // check the candidate who have the Z coordinate closest to the moving slice num
                                    double Z = (double)m_MovingSliceNum;
                                    double min = 99999.9;
                                    typename Interpolator::ContinuousIndexType ix;
                                    for(unsigned int j = 0; j< P1.size(); j++)
                                    {
                                        m_Images[m_MovingImageNum]->TransformPhysicalPointToContinuousIndex(P1[j],ix);
                                        if(fabs(ix[2]-Z)< min)
                                        {
                                            min = fabs(ix[2]-Z);
                                            Point1[0] = P1[j][0];
                                            Point1[1] = P1[j][1];
                                            Point1[2] = P1[j][2];
                                        }
                                    }
                                    FoundP1 = true;

                                    P1.clear();
                                    break;
                                }
                                else if(FoundP1 && P2.size()>0)
                                {
                                    // check the candidate who have the Z coordinate closest to the moving slice num
                                    double Z = (double)m_MovingSliceNum;
                                    double min = 99999.9;
                                    typename Interpolator::ContinuousIndexType ix;
                                    for(unsigned int j = 0; j< P2.size(); j++)
                                    {
                                        m_Images[m_MovingImageNum]->TransformPhysicalPointToContinuousIndex(P2[j],ix);
                                        if(fabs(ix[2]-Z) < min)
                                        {
                                            min = fabs(ix[2]-Z);
                                            Point2[0] = P2[j][0];
                                            Point2[1] = P2[j][1];
                                            Point2[2] = P2[j][2];
                                        }
                                    }
                                    FoundP2 = true;

                                    P2.clear();
                                    break;
                                }
                            }

                        }


                    }

                    if(FoundP1 && FoundP2)
                    {
                        intersection = true;
                        break;
                    }
                }

                //****
                //TODO : Cleaning
                IteratorTab.clear();
                P1.clear();
                P2.clear();
                //****
                if(intersection)
                {



                    typename ImageType::PointType P, Pfixed, Pmoving;
                    itk::Vector<double,3> P12;
                    typename Interpolator::ContinuousIndexType IndexFixed, IndexMoving;
                    typename ImageType::IndexType MaskIndexFx,MaskIndexMv;

                    P12[0] = Point2[0] - Point1[0];
                    P12[1] = Point2[1] - Point1[1];
                    P12[2] = Point2[2] - Point1[2];


                    for(unsigned int i = 0; i<m_NumberOfPointsOfLine; i++)
                    {

                        P = Point1 + (P12 * i/(m_NumberOfPointsOfLine-1));


                        Pfixed = m_InverseTransforms[ifixed]->TransformPoint(P);

                        Pmoving = InverseX->TransformPoint(P);



                        m_Images[ifixed]->TransformPhysicalPointToContinuousIndex(Pfixed, IndexFixed);
                        m_Images[m_MovingImageNum]->TransformPhysicalPointToContinuousIndex(Pmoving,IndexMoving);

                        m_Images[ifixed]->TransformPhysicalPointToIndex(Pfixed, MaskIndexFx);
                        m_Images[m_MovingImageNum]->TransformPhysicalPointToIndex(Pmoving,MaskIndexMv);


                        // Test if still in the image, while EvaluateAtContinuousIndex method don't check if pixel is inside bounding box
                        if(m_Interpolators[ifixed]->IsInsideBuffer(IndexFixed) && m_Interpolators[m_MovingImageNum]->IsInsideBuffer(IndexMoving))
                        {
                            VoxelType fixedVoxel, movingVoxel;

                            if(m_Masks[ifixed]->GetPixel(MaskIndexFx) == 1 && m_Masks[m_MovingImageNum]->GetPixel(MaskIndexMv) == 1)
                            {
                                movingVoxel = m_Interpolators[m_MovingImageNum]->EvaluateAtContinuousIndex(IndexMoving);
                                fixedVoxel = m_Interpolators[ifixed]->EvaluateAtContinuousIndex(IndexFixed);


                                double SquareDiff = SquareDifference(fixedVoxel,movingVoxel);


                                #pragma omp critical
                                //#pragma omp atomic
                                {
                                    SumOfIntersectedVoxels += SquareDiff;
                                    NumberOfIntersectedVoxels++;
                                }
                            }
                            else
                            {

                            }


                        }//interpolate

                    }


                }
                else
                {

                }

            }
        }
    }
    //Normalize the sum
    if(NumberOfIntersectedVoxels != 0)
    {
        CostFunction = (SumOfIntersectedVoxels/(VoxelType)NumberOfIntersectedVoxels) ;
    }
    else
    {
        CostFunction =  MAX_COSTFUNCTION_VALUE;
    }
    if(m_VerboseMode)
    {
        std::cout<<"Number of Intersected Voxels : "<<NumberOfIntersectedVoxels<<std::endl;
        std::cout<<"Sum of difference between intersected voxels : "<<SumOfIntersectedVoxels<<std::endl;
        std::cout<<"CostFunction Value : "<<CostFunction<<std::endl;
    }

    return CostFunction;


}
//-------------------------------------------------------------------------------------------------
template<class TImage>
vnl_vector<double> SlicesIntersectionVNLCostFunction<TImage>::GetGradient(const vnl_vector<double> &x, double stepsize) const
{
    //finite difference
    vnl_vector<double> tx = x;
    vnl_vector<double> gradient(x.size());
    double h = stepsize;
    for (int i = 0; i < dim; ++i)
    {
        double tplus = x[i] + h;
        tx[i] = tplus;
        double fplus = this->f(tx);

        double tminus = x[i] - h;
        tx[i] = tminus;
        double fminus = this->f(tx);

        gradient[i] = (fplus - fminus) / (tplus - tminus);
        tx[i] = x[i];
    }
    return gradient;
}
//-------------------------------------------------------------------------------------------------
}
#endif
