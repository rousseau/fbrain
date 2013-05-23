#ifndef BTKSlicesIntersectionVNLCostFunction_TXX
#define BTKSlicesIntersectionVNLCostFunction_TXX

#include "btkSlicesIntersectionVNLCostFunction.hxx"
#include "btkSliceIntersectionFilter.h"

namespace btk
{

//-------------------------------------------------------------------------------------------------
template<class TImage>
SlicesIntersectionVNLCostFunction<TImage>::SlicesIntersectionVNLCostFunction(unsigned int dim)
    :vnl_cost_function(dim),m_NumberOfPointsOfLine(100),m_VerboseMode(false),m_MovingImageNum(0),m_MovingSliceNum(0)
    ,m_Intersection(false)
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
        unsigned int linearSize =  0;
        for(unsigned int i = 0; i<m_NumberOfImages; i++)
        {
            m_Interpolators[i] = Interpolator::New();
            m_Interpolators[i]->SetInputImage(m_Images[i]);
            linearSize+= m_Images[i]->GetLargestPossibleRegion().GetSize()[2];

        }
        m_Stack.resize(linearSize);


        m_X = TransformType::New();
        m_X->SetIdentity();

        m_InverseX = TransformType::New();
        m_InverseX->SetIdentity();

        typename ImageType::SizeType size;
        size =  m_Images[m_MovingImageNum]->GetLargestPossibleRegion().GetSize();

        // For the center of the transform !
        ContinuousIndexType centerIndex;
        centerIndex[0] = (size[0]-1.0)/2.0 ;
        centerIndex[1] = (size[1]-1.0)/2.0 ;
        centerIndex[2] = m_MovingSliceNum ;

        typename ImageType::PointType center;
        m_Images[m_MovingImageNum]->TransformContinuousIndexToPhysicalPoint(centerIndex,center);
        m_X->SetCenter(center);

        std::pair<unsigned int, unsigned int> RefImSlice;
        RefImSlice.first = 1;
        RefImSlice.second = m_Images[1]->GetLargestPossibleRegion().GetSize()[2]/2;
        m_ReferenceSlice = this->ReturnStackIndexFromImageAndSlice(RefImSlice);

//        btkTicTocInit();
//        btkTic();
//        this->ComputeIntersectionsOfAllSlices();
//        btkToc();



    }


}
//-------------------------------------------------------------------------------------------------
template<class TImage>
vnl_vector<double> SlicesIntersectionVNLCostFunction<TImage>::GetGradient(const vnl_vector<double> &x, double stepsize) const
{
//    //finite difference
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
template< class TImage>
bool
SlicesIntersectionVNLCostFunction<TImage>::
FoundIntersectionPoints(unsigned int _fixedImage, unsigned int _fixedSlice, unsigned int _movingImage, unsigned int _movingSlice, typename ImageType::PointType &Point1, typename ImageType::PointType &Point2) const
{


    // Lets consider the 4 corners of a slice like this :

    //      slice
    //     1------2
    //     |      |
    //     |      |
    //     4------3


    Point1[0] = 0;
    Point1[1] = 0;
    Point1[2] = 0;

    Point2[0] = 0;
    Point2[1] = 0;
    Point2[2] = 0;

    // 4 corner of moving slice:
    typename ImageType::IndexType MovingCorner1, MovingCorner2, MovingCorner3, MovingCorner4;

    // 4 corners of fixed slice :
    typename ImageType::IndexType FixedCorner1, FixedCorner2, FixedCorner3, FixedCorner4;

    typename ImageType::RegionType MovingRegion, FixedRegion;
    typename ImageType::SizeType MovingRgSize, FixedRgSize;


    MovingRgSize = m_Images[_movingImage]->GetLargestPossibleRegion().GetSize();
    MovingRgSize[2] = 1;

    FixedRgSize = m_Images[_fixedImage]->GetLargestPossibleRegion().GetSize();
    FixedRgSize[2] = 1;


    MovingCorner1[0] = 0; MovingCorner1[1] = 0; MovingCorner1[2] = _movingSlice;
    MovingCorner2[0] = 0; MovingCorner2[1] = m_Images[_movingImage]->GetLargestPossibleRegion().GetSize()[1] - 1; MovingCorner2[2] = _movingSlice;
    MovingCorner3[0] = m_Images[_movingImage]->GetLargestPossibleRegion().GetSize()[0] - 1; MovingCorner3[1] = m_Images[_movingImage]->GetLargestPossibleRegion().GetSize()[1] - 1;MovingCorner3[2] = _movingSlice;
    MovingCorner4[0] = m_Images[_movingImage]->GetLargestPossibleRegion().GetSize()[0] - 1;MovingCorner4[1] = 0; MovingCorner4[2] = _movingSlice;

    FixedCorner1[0] = 0; FixedCorner1[1] = 0; FixedCorner1[2] = _fixedSlice;
    FixedCorner2[0] = 0; FixedCorner2[1] = m_Images[_fixedImage]->GetLargestPossibleRegion().GetSize()[1] - 1; FixedCorner2[2] = _fixedSlice;
    FixedCorner3[0] = m_Images[_fixedImage]->GetLargestPossibleRegion().GetSize()[0] - 1; FixedCorner3[1] = m_Images[_fixedImage]->GetLargestPossibleRegion().GetSize()[1] - 1;FixedCorner3[2] = _fixedSlice;
    FixedCorner4[0] = m_Images[_fixedImage]->GetLargestPossibleRegion().GetSize()[0] - 1;FixedCorner4[1] = 0; FixedCorner4[2] = _fixedSlice;

    MovingRegion.SetIndex(MovingCorner1);
    FixedRegion.SetIndex(FixedCorner1);

    MovingRegion.SetSize(MovingRgSize);
    FixedRegion.SetSize(FixedRgSize);


    typename ImageType::PointType CurrentPoint, TCurrentPoint, CorrespondingPoint;
    bool FoundP1, FoundP2, intersection;

    intersection = FoundP1 = FoundP2 = false;

    typename ImageType::IndexType  CurrentIndex;
    typename Interpolator::ContinuousIndexType CorrespondingIndex;

    typedef itk::LineConstIterator<ImageType> LineIterator;

    std::vector<LineIterator> IteratorTab;


    LineIterator itM12(m_Images[_movingImage], MovingCorner1, MovingCorner2);
    LineIterator itM23(m_Images[_movingImage], MovingCorner2, MovingCorner3);
    LineIterator itM34(m_Images[_movingImage], MovingCorner3, MovingCorner4);
    LineIterator itM41(m_Images[_movingImage], MovingCorner4, MovingCorner1);


    // for the fixed slice we take the opposite
    LineIterator itF12(m_Images[_fixedImage], FixedCorner1, FixedCorner4);
    LineIterator itF23(m_Images[_fixedImage], FixedCorner4, FixedCorner3);
    LineIterator itF34(m_Images[_fixedImage], FixedCorner3, FixedCorner2);
    LineIterator itF41(m_Images[_fixedImage], FixedCorner2, FixedCorner1);

    IteratorTab.push_back(itM12);
    IteratorTab.push_back(itM23);
    IteratorTab.push_back(itM34);
    IteratorTab.push_back(itM41);
    IteratorTab.push_back(itF12);
    IteratorTab.push_back(itF23);
    IteratorTab.push_back(itF34);
    IteratorTab.push_back(itF41);


    std::vector<typename ImageType::PointType> P1, P2;
    //TODO: This part is very complex, is there a way to simplify it ?
    // We iterate over each edge of slices for found starting and ending points of intersection
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
                m_Images[_movingImage]->TransformIndexToPhysicalPoint(CurrentIndex, CurrentPoint);
                // Apply X
                TCurrentPoint = m_X->TransformPoint(CurrentPoint);// use m_Transform[_movingImage] instead
                //TCurrentPoint = m_Transforms[_movingImage]->TransformPoint(CurrentPoint);
                // Apply Inverse of the Optimized founded transform (if not identity)
                CorrespondingPoint = m_InverseTransforms[_fixedImage]->TransformPoint(TCurrentPoint);
                // Found the Corresponding Index in the fixed slice !
                m_Images[_fixedImage]->TransformPhysicalPointToContinuousIndex(CorrespondingPoint,CorrespondingIndex );


                // Is in the fixed slice ?
                if(FixedRegion.IsInside(CorrespondingIndex))
                {
                    countPoint++;

                    ConsecutivePoint = true;
                    if(FoundP1 == false && FoundP2 == false)
                    {

                        P1.push_back(TCurrentPoint);


                    }
                    else if(FoundP1 == true && FoundP2 == false)
                    {
//
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
                        // check the candidate who have the Z coordinate closest to the fixed slice num
                        double Z = (double)_fixedSlice;
                        double min = 99999.9;
                        typename Interpolator::ContinuousIndexType ix;
                        for(unsigned int j = 0; j< P1.size(); j++)
                        {
//                            Point1[0] += P1[j][0];
//                            Point1[1] += P1[j][1];
//                            Point1[2] += P1[j][2];
                            m_Images[_fixedImage]->TransformPhysicalPointToContinuousIndex(P1[j],ix);

                            if(fabs(ix[2]-Z)< min)
                            {
                                min = fabs(ix[2]-Z);

                                Point1[0] = P1[j][0];
                                Point1[1] = P1[j][1];
                                Point1[2] = P1[j][2];
                            }
                        }

//                        Point1[0] = Point1[0]/P1.size();
//                        Point1[1] = Point1[1]/P1.size();
//                        Point1[2] = Point1[2]/P1.size();
                        FoundP1 = true;

                        P1.clear();
                        break;
                    }
                    else if(FoundP1 && P2.size() > 0)
                    {
                        // check the candidate who have the Z coordinate closest to the fixed slice num
                        double Z = (double)_fixedSlice;
                        double min = 99999.9;
                        typename Interpolator::ContinuousIndexType ix;
                        for(unsigned int j = 0; j< P2.size(); j++)
                        {

//                            Point2[0] += P2[j][0];
//                            Point2[1] += P2[j][1];
//                            Point2[2] += P2[j][2];


                            m_Images[_fixedImage]->TransformPhysicalPointToContinuousIndex(P2[j],ix);
                            if(fabs(ix[2]-Z) < min)
                            {
                                min = fabs(ix[2]-Z);
                                Point2[0] = P2[j][0];
                                Point2[1] = P2[j][1];
                                Point2[2] = P2[j][2];
                            }
                        }
//                        Point2[0] = Point2[0]/P2.size();
//                        Point2[1] = Point2[1]/P2.size();
//                        Point2[2] = Point2[2]/P2.size();
                        FoundP2 = true;

                        P2.clear();
                        break;
                    }
                }

            } //moving slice

            //TODO : Duplication of code...
            else // fixed slice
            {
                m_Images[_fixedImage]->TransformIndexToPhysicalPoint(CurrentIndex, CurrentPoint);
                //Apply transform of the fixedImage (if not identity)
                TCurrentPoint = m_Transforms[_fixedImage]->TransformPoint(CurrentPoint);
                //Apply inverse of X for return in the moving slice space
                CorrespondingPoint = m_InverseX->TransformPoint(TCurrentPoint); // use m_InverseTransforms instead
                //CorrespondingPoint = m_InverseTransforms[_movingImage]->TransformPoint(TCurrentPoint);
                m_Images[_movingImage]->TransformPhysicalPointToContinuousIndex(CorrespondingPoint, CorrespondingIndex);

                if(MovingRegion.IsInside(CorrespondingIndex))
                {
                    countPoint++;

                    ConsecutivePoint = true;
                    if(FoundP1 == false && FoundP2 == false)
                    {

                        P1.push_back(TCurrentPoint);

                    }
                    else if(FoundP1 == true && FoundP2 == false)
                    {

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
                        double Z = (double)_movingSlice;
                        double min = 99999.9;// for initialization
                        typename Interpolator::ContinuousIndexType ix;
                        for(unsigned int j = 0; j< P1.size(); j++)
                        {
//                            Point1[0] += P1[j][0];
//                            Point1[1] += P1[j][1];
//                            Point1[2] += P1[j][2];

                            m_Images[_movingImage]->TransformPhysicalPointToContinuousIndex(P1[j],ix);
                            if(fabs(ix[2]-Z)< min)
                            {
                                min = fabs(ix[2]-Z);
                                Point1[0] = P1[j][0];
                                Point1[1] = P1[j][1];
                                Point1[2] = P1[j][2];
                            }
                        }
//                        Point1[0] = Point1[0]/P1.size();
//                        Point1[1] = Point1[1]/P1.size();
//                        Point1[2] = Point1[2]/P1.size();
                        FoundP1 = true;

                        P1.clear();
                        break;
                    }
                    else if(FoundP1 && P2.size()>0)
                    {
                        // check the candidate who have the Z coordinate closest to the moving slice num
                        double Z = (double)_movingSlice;
                        double min = 99999.9;// for init
                        typename Interpolator::ContinuousIndexType ix;
                        for(unsigned int j = 0; j< P2.size(); j++)
                        {
//                            Point2[0] += P2[j][0];
//                            Point2[1] += P2[j][1];
//                            Point2[2] += P2[j][2];

                            m_Images[_movingImage]->TransformPhysicalPointToContinuousIndex(P2[j],ix);
                            if(fabs(ix[2]-Z) < min)
                            {
                                min = fabs(ix[2]-Z);
                                Point2[0] = P2[j][0];
                                Point2[1] = P2[j][1];
                                Point2[2] = P2[j][2];
                            }
                        }
//                        Point2[0] = Point2[0]/P2.size();
//                        Point2[1] = Point2[1]/P2.size();
//                        Point2[2] = Point2[2]/P2.size();
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
    // Cleaning
    IteratorTab.clear();
    P1.clear();
    P2.clear();

    if(intersection)
    {
        return true;
    }
    else
    {
        return false;
    }


}
//-------------------------------------------------------------------------------------------------
template<class TImage>
void
SlicesIntersectionVNLCostFunction<TImage>::
GetTransformsWithParams(const vnl_vector<double> &x) const
{
    TransformType::ParametersType params;

    params.SetSize(6);

    for(unsigned slice= 0; slice< m_Stack.size(); slice++)
    {
        std::pair<unsigned int ,unsigned int> imageSlice = this->ReturnImageAndSliceFromStackIndex(slice);
        unsigned int im = imageSlice.first;
        unsigned int sl = imageSlice.second;

        for(unsigned int p = 0; p< 6; p++)
        {
            if(p< 3)
            {
                params[p] = btk::MathFunctions::DegreesToRadians(x[slice *6 +p]);
            }
            else
            {
                params[p] = x[slice *6 +p];
            }
        }

        m_Transforms[im]->SetSliceParameters(sl,params);
    }
}

//-------------------------------------------------------------------------------------------------
template<class TImage>
double SlicesIntersectionVNLCostFunction<TImage>::f(const vnl_vector<double> &x) const

{


    // Initials values
    double        CostFunction = 0.0;
    unsigned long NumberOfIntersectedVoxels = 0;
    double        SumOfIntersectedVoxels = 0.0;
    unsigned int  NumberOfIntersectedSlices = 0;

    TransformType::ParametersType params;
    params.SetSize(x.size());
    params = m_X->GetParameters();
    // Convert Angles in degrees to radians
    for(unsigned int i = 0; i< params.size(); i++)
    {

        if(i < 3)
        {
            params[i] = MathFunctions::DegreesToRadians(x[i]);
        }
        else
          params[i] = x[i];

    }

    if(m_VerboseMode)
    {
        std::cout<<"Input parameters : "<<params<<std::endl;
    }

    for(unsigned int s = 0; s< m_SlicesGroup.size(); s++)
    {
        if(m_SlicesGroup[s] == m_GroupNum)
        {
            m_Transforms[m_MovingImageNum]->SetSliceParameters(s, params);
            m_Transforms[m_MovingImageNum]->GetInverse(m_InverseTransforms[m_MovingImageNum]);
            // Set Optimizer's parameters in transformation
            m_X->SetParameters(params);

            //Get the inverse for the process
            m_InverseX->SetCenter(m_X->GetCenter());
            m_InverseX->SetFixedParameters(m_X->GetFixedParameters());

            m_X->GetInverse(m_InverseX);

            unsigned int ifixed = 0;


            // For all others images
            for( ifixed = 0; ifixed< m_Images.size(); ifixed ++)
            {
                unsigned int numberOfFixedSlices = m_Images[ifixed]->GetLargestPossibleRegion().GetSize()[2];

                // Test if images are orthogonals, if not don't use it !
                bool AreOrthos = btk::ImageHelper<ImageType>::AreOrthos(m_Images[m_MovingImageNum],m_Images[ifixed]);

                // if current image is different of moving image AND if they are orthogonals
                if(ifixed != m_MovingImageNum && AreOrthos)
                {
                    unsigned int sfixed = 0;

                    /* PARALLELIZATION : Use of reduction to avoid use of #pragma omp critical
                     * Each thread have a local copy of both NumberOfIntersectedVoxels and SumOfIntersectedVoxels
                     * At the end we sum all local copys into a global one.
                     **/
                    #pragma omp parallel for private(sfixed) schedule(dynamic)\
                    reduction(+:NumberOfIntersectedVoxels) reduction(+:SumOfIntersectedVoxels)

                    for(sfixed = 0; sfixed < numberOfFixedSlices; sfixed++)
                    {
                        // Check for starting point and ending point
                        typename ImageType::PointType Point1, Point2;

                        /* Point 1 is the starting point, point 2 the ending one, if there is no intersection the function return false
                        the function is looking for an intersection between fixed slice (sfixed) in fixed image (ifixed)
                        and the moving one (m_MovingImageNum,m_MovingSliceNum) */

                        //TODO: we can improve computed time by precomputing intersections points outside this loop (Point1, Point2 and intersection will be vectors)
                        //Not sure if we loose time or not...
                        //ITK
                        bool intersection = this->FoundIntersectionPoints(ifixed,sfixed,m_MovingImageNum,s,Point1,Point2);


                        //VTK
//                        bool intersection = false;

//                            typedef SliceIntersectionFilter<ImageType> intersectionFilter;
//                            typename intersectionFilter::Pointer filter = intersectionFilter::New();
//                            filter->SetImage1(m_Images[ifixed]);
//                            filter->SetImage2(m_Images[m_MovingImageNum]);
//                            filter->SetSlice1(sfixed);
//                            filter->SetSlice2(s);
//                            filter->SetTransform1(m_Transforms[ifixed]);
//                            filter->SetTransform2(m_X);

//                            filter->Update();

//                            intersection = filter->GetIsIntersection();
//                            Point1 = filter->GetPoint1();
//                            Point2 = filter->GetPoint2();






                        // If we have found an intersection
                        if(intersection)
                        {

                            NumberOfIntersectedSlices++;

                            typename ImageType::PointType P, Pfixed, Pmoving,Pfirst,Plast;
                            itk::Vector<double,3> P12;
                            typename Interpolator::ContinuousIndexType IndexFixed, IndexMoving;
                            typename ImageType::IndexType MaskIndexFx,MaskIndexMv;

        //                    P12[0] = Point2[0] - Point1[0];
        //                    P12[1] = Point2[1] - Point1[1];
        //                    P12[2] = Point2[2] - Point1[2];

                            //Looking for the first point and the last point :
                            double p1 = std::sqrt((Point1[0]* Point1[0]) + (Point1[1]* Point1[1]) + (Point1[2]* Point1[2]));
                            double p2 = std::sqrt((Point2[0]* Point2[0]) + (Point2[1]* Point2[1]) + (Point2[2]* Point2[2]));

                            if(p1 > p2)
                            {
                                Pfirst = Point2;
                                Plast = Point1;
                               //
                                P12[0] = Point1[0] - Point2[0];
                                P12[1] = Point1[1] - Point2[1];
                                P12[2] = Point1[2] - Point2[2];

                            }
                            else
                            {
                                Pfirst = Point1;
                                Plast = Point2;
                                //
                                P12[0] = Point2[0] - Point1[0];
                                P12[1] = Point2[1] - Point1[1];
                                P12[2] = Point2[2] - Point1[2];

                            }

                            // Compute the distance between the 2 extrema (points)
                            double EuclideanDist = std::sqrt((P12[0]* P12[0]) + (P12[1]* P12[1]) + (P12[2]* P12[2]));
                            // Set the number of points equal to this distance times the spacing in mm (in x direction)
                            int NumberOfPoints = std::floor(EuclideanDist);// * m_Images[ifixed]->GetSpacing()[0];
                            //std::cout<<"Using "<<NumberOfPoints<<" intersection points."<<std::endl;


                            // Loop Over intersection line points
                            //for(unsigned int i = 0; i<m_NumberOfPointsOfLine; i++)// user set a fixed number of points
                            for(unsigned int i=0; i<NumberOfPoints;i++)// a point each mm
                            //for(float i=0.0; i<NumberOfPoints; i=i+space)// a point each spacing
                            {
                                //P = Point1 + (P12 * i/(m_NumberOfPointsOfLine-1));//if m_NumberOfPointsOfLine is used

                                //Compute the point
                                //P = Point1 + (P12 * i/(NumberOfPoints));
                                P = Pfirst + (P12 * i/(NumberOfPoints));


                                //std::cout<<"P-->"<<P<<std::endl;

                                //On peut accélerer ce calcul en travaillant dans l'espace image (?) de l'image moving
                                //en optimisant également le nombre de test sur les masques

                                // Point in fixed image
                                Pfixed = m_InverseTransforms[ifixed]->TransformPoint(P);
                                //Point in moving image
                                Pmoving = m_InverseX->TransformPoint(P);

                                m_Images[ifixed]->TransformPhysicalPointToContinuousIndex(Pfixed, IndexFixed);
                                m_Images[m_MovingImageNum]->TransformPhysicalPointToContinuousIndex(Pmoving,IndexMoving);

                                m_Images[ifixed]->TransformPhysicalPointToIndex(Pfixed, MaskIndexFx);
                                m_Images[m_MovingImageNum]->TransformPhysicalPointToIndex(Pmoving,MaskIndexMv);

                                // Test if still in the image, while EvaluateAtContinuousIndex method don't check if pixel is inside bounding box
                                if(m_Interpolators[ifixed]->IsInsideBuffer(IndexFixed) && m_Interpolators[m_MovingImageNum]->IsInsideBuffer(IndexMoving))
                                {
                                    VoxelType fixedVoxel, movingVoxel;

                                    //Test if pixel is inside both fixed maks and moving mask
                                    if(m_Masks[ifixed]->GetPixel(MaskIndexFx) > 0 && m_Masks[m_MovingImageNum]->GetPixel(MaskIndexMv) > 0)
                                    {
                                        // Get The values with linear interpolator
                                        movingVoxel = m_Interpolators[m_MovingImageNum]->EvaluateAtContinuousIndex(IndexMoving);
                                        fixedVoxel = m_Interpolators[ifixed]->EvaluateAtContinuousIndex(IndexFixed);

                                        //SumOfIntersectedVoxels += SquaredDifference(fixedVoxel,movingVoxel); // MSE
                                        SumOfIntersectedVoxels += AbsoluteDifference(fixedVoxel,movingVoxel); //MAE
                                        //SumOfIntersectedVoxels += RootSquaredDifference(fixedVoxel,movingVoxel); //RSE
                                        NumberOfIntersectedVoxels++;// we count each point

                                        //std::cout<<"Number of Points:"<<NumberOfIntersectedVoxels<<std::endl;


                                    }//End if inside masks
                                }//End if IsInsideBuffer
                            }//End Loop over POINT

                        }//End If intersection
                    }//end of loop over SLICES
                }// end if image = reference
            }// end loop over IMAGES

            //Normalize the sum
            if(NumberOfIntersectedVoxels != 0)
            {
                CostFunction += (SumOfIntersectedVoxels/(double)NumberOfIntersectedVoxels *1.0) ;

            }

            if(m_VerboseMode)
            {
                std::cout<<"Number of Intersected Slices : "<<NumberOfIntersectedSlices<<std::endl;
                std::cout<<"Number of Intersected Voxels : "<<NumberOfIntersectedVoxels<<std::endl;
                std::cout<<"Sum of difference between intersected voxels : "<<SumOfIntersectedVoxels<<std::endl;
                std::cout<<"CostFunction Value : "<<CostFunction<<std::endl;
            }


        }
    }

    if(NumberOfIntersectedVoxels == 0)
    {
        CostFunction = 10e6; //cost function can not be equal to 0, it is the min value !!!
        m_Intersection = false;
    }
    else
    {
        m_Intersection = true;
    }
    return CostFunction;
    //return CostFunction * CostFunction; //squared


}
//-------------------------------------------------------------------------------------------------
template<class TImage>
void SlicesIntersectionVNLCostFunction<TImage>::ComputeIntersectionsOfAllSlices()
{
    m_Points1.resize(m_Stack.size());
    m_Points2.resize(m_Stack.size());
    m_Intersections.resize(m_Stack.size());

   // std::cout<<""<<m_Stack.size()<<std::endl;

    typename TransformType::Pointer IdTransform = TransformType::New();
    IdTransform->SetIdentity();

    //Initialization
    for(unsigned int i =0; i< m_Stack.size(); i++)
    {
        m_Points1[i].resize(m_Stack.size());
        m_Points2[i].resize(m_Stack.size());
        m_Intersections[i].resize(m_Stack.size());
        std::fill(m_Intersections[i].begin(), m_Intersections[i].end(),false);

    }



    for(unsigned int i = 0; i< m_Stack.size(); i++)
    {


        for(unsigned int j = 0; j< m_Stack.size(); j++)
        {
            std::pair< unsigned int, unsigned int > ImageSlice1;
            ImageSlice1 = this->ReturnImageAndSliceFromStackIndex(i);
            std::pair< unsigned int, unsigned int > ImageSlice2;
            ImageSlice2 = this->ReturnImageAndSliceFromStackIndex(j);

            bool AreOrthos =  btk::ImageHelper< ImageType>::AreOrthos(m_Images[ImageSlice1.first],m_Images[ImageSlice2.first]);
//            std::cout<<"Intersection between : "<<i<<" and : "<<j<<std::endl;

//            std::cout<<i<<" is image : "<<ImageSlice1.first<<" and slice : "<<ImageSlice1.second<<std::endl;
//            std::cout<<j<<" is image : "<<ImageSlice2.first<<" and slice : "<<ImageSlice2.second<<std::endl;

            if(i != j && AreOrthos && !m_Intersections[i][j])
            {
                bool intersection = false;
                typedef SliceIntersectionFilter<ImageType> intersectionFilter;
                typename intersectionFilter::Pointer filter = intersectionFilter::New();

                filter->SetImage1(m_Images[ImageSlice1.first]);
                filter->SetImage2(m_Images[ImageSlice2.first]);
                filter->SetSlice1(ImageSlice1.second);
                filter->SetSlice2(ImageSlice2.second);
                filter->SetTransform1(IdTransform);
                filter->SetTransform2(IdTransform);

                filter->Update();

                intersection = filter->GetIsIntersection();
                if(intersection)
                {
                    m_Intersections[i][j] = true;
                    m_Points1[i][j] = filter->GetPoint1();
                    m_Points2[i][j] = filter->GetPoint2();

                    //Do the inverse

                    m_Intersections[j][i] = true;

                    m_Points1[j][i] = filter->GetPoint1();
                    m_Points2[j][i] = filter->GetPoint2();

                    //std::cout<<"There is an intersection with points : "<<m_Points1[i][j]<<" and "<<m_Points2[i][j]<<std::endl;
                }
                else
                {
                    m_Intersections[i][j] = false;
                   // std::cout<<"There is no intersection"<<std::endl;
                }

            }
            else
            {
                m_Intersections[i][j]=false;
               // std::cout<<"No intersection while image and slice are the same or Orthos  = "<<AreOrthos<<std::endl;
            }
        }
    }
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
}
#endif
