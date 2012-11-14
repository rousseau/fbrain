#ifndef BTKSlicesIntersectionVNLCostFunction_TXX
#define BTKSlicesIntersectionVNLCostFunction_TXX

#include "btkSlicesIntersectionVNLCostFunction.hxx"
namespace btk
{
static int count = 0;
//-------------------------------------------------------------------------------------------------
template<class TImage>
SlicesIntersectionVNLCostFunction<TImage>::SlicesIntersectionVNLCostFunction(unsigned int dim)
    :vnl_cost_function(dim),m_NumberOfPointsOfLine(1000),m_VerboseMode(false)
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


        m_Transform = TransformType::New();
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
    int NumberOfIntersectedVoxels = 0;
    double SumOfIntersectedVoxels = 0;



    TransformType::ParametersType params, refparams;
    params.SetSize(x.size());
    refparams.SetSize(x.size());


    for(unsigned int i = 0; i< x.size(); i++)
    {
        if(i < 3)
        {
            params[i] = MathFunctions::DegreesToRadians(x[i]);
        }
        else
        {
            params[i] = x[i];
        }

    }

    if(m_VerboseMode)
    {
        std::cout<<"Parameters : "<<x<<std::endl;
    }

    m_Transform->SetParameters(params); // Parameters found by vnl_otpimizers (powell) and to apply to the current moving slice


    TransformType::Pointer InverseX = TransformType::New();
    m_Transform->GetInverse(InverseX);

    //TODO: paralelize this...

    for(unsigned int ifixed = 0; ifixed< m_Images.size(); ifixed ++)
    {
        if(ifixed != m_MovingImageNum)
        {
            unsigned int numberOfFixedSlices = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[2];

            // TODO : Check if it is possible to paralelize this :
            for(unsigned int sfixed = 0; sfixed < numberOfFixedSlices; sfixed++)
            {


                //      slice
                //     1------2
                //     |      |
                //     |      |
                //     4------3
                // 4 corner of moving image:
                typename ImageType::IndexType MovingCorner1, MovingCorner2, MovingCorner3, MovingCorner4;

                // 4 corners of fixed image :
                typename ImageType::IndexType FixedCorner1, FixedCorner2, FixedCorner3, FixedCorner4;

                MovingCorner1[0] = 0; MovingCorner1[1] = 0; MovingCorner1[2] = m_MovingSliceNum;
                MovingCorner2[0] = 0; MovingCorner2[1] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[1] - 1; MovingCorner2[2] = m_MovingSliceNum;
                MovingCorner3[0] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[0] - 1; MovingCorner3[1] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[1] - 1;MovingCorner3[2] = m_MovingSliceNum;
                MovingCorner4[0] = m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize()[0] - 1;MovingCorner4[1] = 0; MovingCorner4[2] = m_MovingSliceNum;

                FixedCorner1[0] = 0; FixedCorner1[1] = 0; FixedCorner1[2] = sfixed;
                FixedCorner2[0] = 0; FixedCorner2[1] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[1] - 1; FixedCorner2[2] = sfixed;
                FixedCorner3[0] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[0] - 1; FixedCorner3[1] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[1] - 1;FixedCorner3[2] = sfixed;
                FixedCorner4[0] = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[0] - 1;FixedCorner4[1] = 0; FixedCorner4[2] = sfixed;


                typename ImageType::PointType Point1, Point2, CurrentPoint,CorrespondingPoint;
                bool FoundP1, FoundP2, intersection;

                intersection = FoundP1 = FoundP2 = false;

                typename ImageType::IndexType  CurrentIndex;


                // FIXME : Use of uninitialized value of type 8 (valgrind)
                typedef itk::LineConstIterator<ImageType> LineIterator;

                std::vector<LineIterator> IteratorTab;


                LineIterator itM12(m_Images[m_MovingImageNum], MovingCorner1, MovingCorner2);
                LineIterator itM23(m_Images[m_MovingImageNum], MovingCorner2, MovingCorner3);
                LineIterator itM34(m_Images[m_MovingImageNum], MovingCorner3, MovingCorner4);
                LineIterator itM41(m_Images[m_MovingImageNum], MovingCorner4, MovingCorner1);

                LineIterator itF12(m_Images[ifixed], FixedCorner1, FixedCorner2);
                LineIterator itF23(m_Images[ifixed], FixedCorner2, FixedCorner3);
                LineIterator itF34(m_Images[ifixed], FixedCorner3, FixedCorner4);
                LineIterator itF41(m_Images[ifixed], FixedCorner4, FixedCorner1);

                IteratorTab.push_back(itM12);
                IteratorTab.push_back(itM23);
                IteratorTab.push_back(itM34);
                IteratorTab.push_back(itM41);
                IteratorTab.push_back(itF12);
                IteratorTab.push_back(itF23);
                IteratorTab.push_back(itF34);
                IteratorTab.push_back(itF41);


                for(unsigned int i = 0; i< IteratorTab.size(); i++)
                {


                    LineIterator it = IteratorTab[i];
                    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
                    {
                        CurrentIndex = it.GetIndex();

                        // Iterate over corners of moving slice
                        if(i<4)
                        {
                            m_Images[m_MovingImageNum]->TransformIndexToPhysicalPoint(CurrentIndex, CurrentPoint);
                            // Apply X
                            CurrentPoint = m_Transform->TransformPoint(CurrentPoint);
                            // Apply Inverse of the Optimized founded transform (if not identity)
                            CorrespondingPoint = m_InverseTransforms[ifixed]->TransformPoint(CurrentPoint);

                            if(m_Interpolators[ifixed]->IsInsideBuffer(CorrespondingPoint))
                            {

                                if(FoundP1 == false && FoundP2 == false)
                                {
                                    Point1 = CurrentPoint;
                                    FoundP1 = true;
                                    break;
                                }
                                else if(FoundP1 == true && FoundP2 == false)
                                {
                                    Point2 = CurrentPoint;
                                    FoundP2 = true;
                                    break;
                                }
                                else
                                {

                                }
                            }


                        }
                        //Iterate over corners of fixed slice
                        else
                        {
                            m_Images[ifixed]->TransformIndexToPhysicalPoint(CurrentIndex, CurrentPoint);
                            //Apply  of the Optimized founded transform (if not identity)
                            CurrentPoint = m_Transforms[ifixed]->TransformPoint(CurrentPoint);
                            //Apply inverse of X for return in the moving slice space
                            CorrespondingPoint = InverseX->TransformPoint(CurrentPoint);
                            if(m_Interpolators[m_MovingImageNum]->IsInsideBuffer(CorrespondingPoint))
                            {

                                if(FoundP1 == false && FoundP2 == false)
                                {
                                    Point1 = CurrentPoint;
                                    FoundP1 = true;
                                    break;
                                }
                                else if(FoundP1 == true && FoundP2 == false)
                                {
                                    Point2 = CurrentPoint;
                                    FoundP2 = true;

                                    break;
                                }
                                else
                                {

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

                        P = Point1 + (P12 * i/m_NumberOfPointsOfLine);

                        //NOTE : Maybe there is a bug with the inverse transform of a sliceBySlice Transform !!
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


                                VoxelType SquareDiff = SquareDifference(fixedVoxel,movingVoxel);

                                SumOfIntersectedVoxels += SquareDiff;
                                NumberOfIntersectedVoxels++;
                            }
                            else
                            {

                            }


                        }//interpolate

                    }


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
}
#endif
