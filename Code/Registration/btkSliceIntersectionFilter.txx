#include "btkSliceIntersectionFilter.h"
#include "vtkTriangleFilter.h"
#include "vtkPlaneSource.h"
#include "vtkIdList.h"

#include "vtkIntersectionPolyDataFilter.h" // Only in version >= vtk 5.9


namespace btk
{
//-------------------------------------------------------------------------------------------------
template<typename TImage>
SliceIntersectionFilter<TImage>::SliceIntersectionFilter()
{
    m_IsIntersection = false;
    m_Slice1 = 0;
    m_Slice2 = 0;
    m_Plane1 = NULL;
    m_Plane2 = NULL;


    m_Line = NULL;

}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void
SliceIntersectionFilter<TImage>::Initialize()
{
    if(m_Image1.IsNull() || m_Image2.IsNull())
    {
        btkException("Error : Image 1 and/or image 2 is NULL !");
    }
    if(m_Transform1 == NULL || m_Transform2 == NULL)
    {
        btkException("Error : Transformations are NULL !");
    }

    //Initialization of objects :

    m_Plane1 = vtkSmartPointer< vtkPolyData >::New();
    m_Plane2 = vtkSmartPointer< vtkPolyData >::New();

    m_Line = vtkSmartPointer< vtkPolyData >::New();

    m_Point1[0] = 0.0;
    m_Point1[1] = 0.0;
    m_Point1[2] = 0.0;

    m_Point2[0] = 0.0;
    m_Point2[1] = 0.0;
    m_Point2[2] = 0.0;



}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void
SliceIntersectionFilter<TImage>::Update()
{
    this->Initialize();
    this->CreatePlane1();
    this->CreatePlane2();

    vtkSmartPointer<vtkIntersectionPolyDataFilter> intersectionPolyDataFilter =
            vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
    intersectionPolyDataFilter->SetInput( 0, m_Plane1 );
    intersectionPolyDataFilter->SetInput( 1, m_Plane2 );
    intersectionPolyDataFilter->Update();

    m_Line = intersectionPolyDataFilter->GetOutput(0); // Output(0) is the intersection Line


    vtkIdType numberOfPointsInCell;
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();

    for(int i = 0; i<m_Line->GetNumberOfCells();i++)
    {
        m_Line->GetCellPoints(i, pointIdList);

        numberOfPointsInCell = pointIdList->GetNumberOfIds();

        for(int j = 0; j< numberOfPointsInCell; j++)
        {
            double point[3];
            std::vector<double> vpoint(3);
            m_Line->GetPoint(pointIdList->GetId(j),point);

            vpoint[0] = point[0];
            vpoint[1] = point[1];
            vpoint[2] = point[2];

            m_LinePoints.push_back(vpoint);

        }

    }

    if(m_LinePoints.size() == 0)
    {
        m_IsIntersection = false;
    }
    else
    {
        m_IsIntersection = true;
        // Check the starting and ending point of the line :
        std::vector< double > refPoint;
        refPoint.resize(3);
        refPoint[0] = m_LinePoints[0][0];
        refPoint[1] = m_LinePoints[0][1];
        refPoint[2] = m_LinePoints[0][2];

        double distMax = 0.0;
        double dist = 0.0;
        int indexOfRefPoints = 0;

        // looking for point wich has the max distance with refPoint
        for(int i = 0; i< m_LinePoints.size(); i++)
        {
            dist = (refPoint[0] - m_LinePoints[i][0]) * (refPoint[0] - m_LinePoints[i][0])+
                    (refPoint[1] - m_LinePoints[i][1]) * (refPoint[1] - m_LinePoints[i][1]) +
                    (refPoint[2] - m_LinePoints[i][2]) * (refPoint[2] - m_LinePoints[i][2]);
            if(dist >= distMax)
            {
                indexOfRefPoints = i;
                distMax = dist;
            }
        }

        std::vector< double > BeginPoint, EndPoint;
        BeginPoint.resize(3);
        EndPoint.resize(3);

       // std::cout<<"index of point : "<<indexOfRefPoints<<std::endl;
        BeginPoint[0] = m_LinePoints[indexOfRefPoints][0];
        BeginPoint[1] = m_LinePoints[indexOfRefPoints][1];
        BeginPoint[2] = m_LinePoints[indexOfRefPoints][2];

        distMax = 0.0;
        dist = 0.0;
        // looking for point wich has the max distance with refPoint
        for(int i = 0; i< m_LinePoints.size(); i++)
        {
            dist = (BeginPoint[0] - m_LinePoints[i][0]) * (BeginPoint[0] - m_LinePoints[i][0])+
                    (BeginPoint[1] - m_LinePoints[i][1]) * (BeginPoint[1] - m_LinePoints[i][1]) +
                    (BeginPoint[2] - m_LinePoints[i][2]) * (BeginPoint[2] - m_LinePoints[i][2]);

            //std::cout<<dist<<std::endl;
            if(dist >= distMax)
            {
                EndPoint[0] = m_LinePoints[i][0];
                EndPoint[1] = m_LinePoints[i][1];
                EndPoint[2] = m_LinePoints[i][2];
                distMax = dist;
            }
        }

        //std::cout<<"Begin Point : "<<BeginPoint[0]<<" "<<BeginPoint[1]<<" "<<BeginPoint[2]<<std::endl;
        //std::cout<<"Ending Point : "<<EndPoint[0]<<" "<<EndPoint[1]<<" "<<EndPoint[2]<<std::endl;

        m_Point1[0] = BeginPoint[0];
        m_Point1[1] = BeginPoint[1];
        m_Point1[2] = BeginPoint[2];

        m_Point2[0] = EndPoint[0];
        m_Point2[1] = EndPoint[1];
        m_Point2[2] = EndPoint[2];

    }


}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void
SliceIntersectionFilter<TImage>::CreatePlane1()
{
    // Create PolyDataPlane for image 1 :
    typename ImageType::SizeType size = m_Image1->GetLargestPossibleRegion().GetSize();
    typename ImageType::IndexType corner1, corner2, corner3, corner4, normal;

    corner1[0] = 0; corner1[1] = 0; corner1[2] = m_Slice1;
    corner2[0] = 0; corner2[1] = size[1]-1; corner2[2] = m_Slice1;
    corner3[0] = size[0]-1; corner3[1] = size[1] -1 ; corner3[2] = m_Slice1;
    corner4[0] = size[0]-1; corner4[1] = 0; corner4[2] = m_Slice1;
    normal[0] = 0; normal[1] = 0;normal[2] = 1;

    typename ImageType::PointType point1, point2, point3, point4, Vnormal;
    typename ImageType::PointType Tpoint1, Tpoint2, Tpoint3, Tpoint4;

    m_Image1->TransformIndexToPhysicalPoint(corner1,point1);
    m_Image1->TransformIndexToPhysicalPoint(corner2,point2);
    m_Image1->TransformIndexToPhysicalPoint(corner3,point3);
    m_Image1->TransformIndexToPhysicalPoint(corner4,point4);
    m_Image1->TransformIndexToPhysicalPoint(normal,Vnormal);

    Tpoint1 = m_Transform1->TransformPoint(point1);
    Tpoint2 = m_Transform1->TransformPoint(point2);
    Tpoint3 = m_Transform1->TransformPoint(point3);
    Tpoint4 = m_Transform1->TransformPoint(point4);

    vtkSmartPointer< vtkPlaneSource > planeSource = vtkSmartPointer< vtkPlaneSource >::New();
    planeSource->SetOrigin(Tpoint1[0],Tpoint1[1],Tpoint1[2]);
    planeSource->SetPoint1(Tpoint2[0],Tpoint2[1],Tpoint2[2]);
    planeSource->SetPoint2(Tpoint4[0],Tpoint4[1],Tpoint4[2]);
    //planeSource->SetNormal(Vnormal[0],Vnormal[1],Vnormal[2]);
    //planeSource->Update();

    //m_Plane1 = planeSource->GetOutput();

    // transform in triangle
    vtkSmartPointer< vtkTriangleFilter > triangleFilter = vtkSmartPointer< vtkTriangleFilter >::New();
    triangleFilter->SetInput(planeSource->GetOutput());
    triangleFilter->Update();
    m_Plane1 = triangleFilter->GetOutput();

}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void
SliceIntersectionFilter<TImage>::CreatePlane2()
{
    // Create PolyDataPlane for image 2 :
    typename ImageType::SizeType size = m_Image2->GetLargestPossibleRegion().GetSize();
    typename ImageType::IndexType corner1, corner2, corner3, corner4, normal;

    corner1[0] = 0; corner1[1] = 0; corner1[2] = m_Slice2;
    corner2[0] = 0; corner2[1] = size[1]-1; corner2[2] = m_Slice2;
    corner3[0] = size[0]-1; corner3[1] = size[1] -1 ; corner3[2] = m_Slice2;
    corner4[0] = size[0]-1; corner4[1] = 0; corner4[2] = m_Slice2;
    normal[0] = 0; normal[1] = 0;normal[2] = 1;

    typename ImageType::PointType point1, point2, point3, point4, Vnormal;
    typename ImageType::PointType Tpoint1, Tpoint2, Tpoint3, Tpoint4;

    m_Image2->TransformIndexToPhysicalPoint(corner1,point1);
    m_Image2->TransformIndexToPhysicalPoint(corner2,point2);
    m_Image2->TransformIndexToPhysicalPoint(corner3,point3);
    m_Image2->TransformIndexToPhysicalPoint(corner4,point4);
    m_Image2->TransformIndexToPhysicalPoint(normal,Vnormal);

    Tpoint1 = m_Transform2->TransformPoint(point1);
    Tpoint2 = m_Transform2->TransformPoint(point2);
    Tpoint3 = m_Transform2->TransformPoint(point3);
    Tpoint4 = m_Transform2->TransformPoint(point4);

    vtkSmartPointer< vtkPlaneSource > planeSource = vtkSmartPointer< vtkPlaneSource >::New();
    planeSource->SetOrigin(Tpoint1[0],Tpoint1[1],Tpoint1[2]);
    planeSource->SetPoint1(Tpoint2[0],Tpoint2[1],Tpoint2[2]);
    planeSource->SetPoint2(Tpoint4[0],Tpoint4[1],Tpoint4[2]);
    //planeSource->SetNormal(Vnormal[0],Vnormal[1],Vnormal[2]);
    //planeSource->Update();

    //m_Plane2 = planeSource->GetOutput();

    // transform in triangle
    vtkSmartPointer< vtkTriangleFilter > triangleFilter = vtkSmartPointer< vtkTriangleFilter >::New();
    triangleFilter->SetInput(planeSource->GetOutput());
    triangleFilter->Update();
    m_Plane2 = triangleFilter->GetOutput();
}
}
