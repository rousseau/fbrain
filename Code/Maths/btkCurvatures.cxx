
// BTK includes
#include "btkCurvatures.h"

// VTK includes
#include "vtkVersion.h"
#include "vtkPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "vtkTriangle.h"
#include "vtkMath.h"
#include "vtkTensor.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"

// VTK standard use for filters implementation
vtkCxxRevisionMacro(btkCurvatures, "$Revision: 1.15 $");
vtkStandardNewMacro(btkCurvatures);

btkCurvatures::btkCurvatures()
{
    this->CurvatureType = CURVATURE_TENSOR;
    this->Tolerance = 1.0e-3;
}
//---------------------------------------------------------
void btkCurvatures::GetCurvatureTensor(vtkPolyData *input)
{
    vtkDebugMacro("Start btkCurvatures::GetCurvatureTensor");

    // empty array check
    if( input->GetNumberOfPolys() == 0 || input->GetNumberOfPoints()==0)
    {
        vtkErrorMacro("No points or cells to operate on");
        return;
    }

    // Create and allocate data for computation
    input->BuildLinks();
    int numberOfVerticies = input->GetNumberOfPoints();
    vtkSmartPointer<vtkIdList> neighboursList = vtkSmartPointer<vtkIdList>::New();
    unsigned short numberOfCellsContainingVertex;
    vtkIdType* cellIds, *points, numberOfPointsInCell;
    double vertex_i[3], vertex_j[3], vertex_k[3];
    vtkDataArray *normals = input->GetPointData()->GetNormals();
    double normalVectorToVertex[3];
    double triangleArea = 0.0;
    double totalTriangleAreas = 0.0;
    double edge_ij[3];
    double projectionVector_ij[3];
    double dotProductOfNormalsToVertex;
    double directionalCurv_ij = 0.0;
    double** M_vi; // matrix of contribution of each edge
    M_vi = new double*[3];
    M_vi[0] = new double[3];
    M_vi[1] = new double[3];
    M_vi[2] = new double[3];
    double *eigenvalues;
    eigenvalues = new double[3];
    double** eigenvectors;
    eigenvectors = new double*[3];
    eigenvectors[0] = new double[3];
    eigenvectors[1] = new double[3];
    eigenvectors[2] = new double[3];
    double eigenvec0[3];
    double eigenvec1[3];
    double eigenvec2[3];
    double c0[3], c1[3], c2[3];
    double lambda1 = 0.0, lambda2=0.0; // eigenvalues
    double K_1 = 0.0; // maximal direction curvature
    double K_2 = 0.0; // minimal direction curvature
    double S_i[3][3]; // curvature tensor
    double T_1[3], T_2[3]; // the two remaining eigenvectors of M_vi defining the principal directions
    vtkDoubleArray* curvatureTensor = vtkDoubleArray::New();
    curvatureTensor->SetName("Curvature Tensor");
    curvatureTensor->SetNumberOfComponents(9);
    curvatureTensor->SetNumberOfTuples(numberOfVerticies);
    vtkDoubleArray* maxCurvature = vtkDoubleArray::New();
    maxCurvature->SetName("Maximum Curvature Vector");
    maxCurvature->SetNumberOfComponents(3);
    maxCurvature->SetNumberOfTuples(numberOfVerticies);
    vtkDoubleArray* minCurvature = vtkDoubleArray::New();
    minCurvature->SetName("Minimum Curvature Vector");
    minCurvature->SetNumberOfComponents(3);
    minCurvature->SetNumberOfTuples(numberOfVerticies);
    vtkDoubleArray* principalCurvatures = vtkDoubleArray::New();
    principalCurvatures->SetName("Principal Curvatures");
    principalCurvatures->SetNumberOfComponents(2);
    principalCurvatures->SetNumberOfTuples(numberOfVerticies);

    // main loop
    for(int i=0; i<numberOfVerticies; i++)
    {

        input->GetPointCells(i,numberOfCellsContainingVertex, cellIds); // get all cells containing i
        input->GetPoint(i,vertex_i);
        normals->GetTuple(i,normalVectorToVertex);

        // Reset the neighbours list and the Matrix M_vi for each vertex i
        neighboursList->Reset();
        for(int n=0; n<3; n++)
        {
            for(int m=0; m<3; m++)
            {
                M_vi[n][m] = 0.0;
            }
        }

        // find the neighbours of i
        for(int c=0; c<numberOfCellsContainingVertex; c++)
        {
            input->GetCellPoints(cellIds[c], numberOfPointsInCell, points);
            for(int j=0; j<numberOfPointsInCell; j++)
            {
                if(points[j] != i)
                {
                    neighboursList->InsertNextId(points[j]);
                }
            }
        }

        // compute the edge projection on the tangeant plane to normal to vertex_i
        for(int j=0; j<(neighboursList->GetNumberOfIds()); j++ )
        {
            if(j == neighboursList->GetNumberOfIds()-1)
            {
                input->GetPoint(neighboursList->GetId(j), vertex_j);
                input->GetPoint(neighboursList->GetId(0), vertex_k);
            }
            else
            {
                input->GetPoint(neighboursList->GetId(j), vertex_j);
                input->GetPoint(neighboursList->GetId(j+1), vertex_k);
            }
            triangleArea = vtkTriangle::TriangleArea(vertex_i, vertex_j, vertex_k);
            totalTriangleAreas += triangleArea;
            edge_ij[0] = vertex_j[0]-vertex_i[0];
            edge_ij[1] = vertex_j[1]-vertex_i[1];
            edge_ij[2] = vertex_j[2]-vertex_i[2];
            // Projection of edge_ij onto the tangeant plane at the normal to vertex_i
            dotProductOfNormalsToVertex = vtkMath::Dot(normalVectorToVertex,edge_ij);
            projectionVector_ij[0] = edge_ij[0] - normalVectorToVertex[0]*dotProductOfNormalsToVertex;
            projectionVector_ij[1] = edge_ij[1] - normalVectorToVertex[1]*dotProductOfNormalsToVertex;
            projectionVector_ij[2] = edge_ij[2] - normalVectorToVertex[2]*dotProductOfNormalsToVertex;
            vtkMath::Normalize(projectionVector_ij);
            // normal curvature in direction ij
            directionalCurv_ij = 2.0*dotProductOfNormalsToVertex/vtkMath::Norm(edge_ij);
            // matrix M_vi = directionalCurv_ij*totalArea_ij*projectionVector*projectionVector^t
            M_vi[0][0] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[0]*projectionVector_ij[0];
            M_vi[1][0] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[1]*projectionVector_ij[0];
            M_vi[2][0] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[2]*projectionVector_ij[0];

            M_vi[0][1] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[0]*projectionVector_ij[1];
            M_vi[1][1] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[1]*projectionVector_ij[1];
            M_vi[2][1] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[2]*projectionVector_ij[1];

            M_vi[0][2] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[0]*projectionVector_ij[2];
            M_vi[1][2] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[1]*projectionVector_ij[2];
            M_vi[2][2] += totalTriangleAreas*directionalCurv_ij*projectionVector_ij[2]*projectionVector_ij[2];
        }


        // looking for eigenvectors and eigenvalues of M_vi, diagonalizing M_vi (Householder Matrix of M_vi has one eigenvector colinear to N_vi)
        vtkMath::Jacobi(M_vi,eigenvalues,eigenvectors); // gives normalized eigenvals and eigenvecs of the sym matrix M_vi in decreasing order
        eigenvec0[0] = eigenvectors[0][0]; eigenvec0[1] = eigenvectors[1][0]; eigenvec0[2] = eigenvectors[2][0];
        eigenvec1[0] = eigenvectors[0][1]; eigenvec1[1] = eigenvectors[1][1]; eigenvec1[2] = eigenvectors[2][1];
        eigenvec2[0] = eigenvectors[0][2]; eigenvec2[1] = eigenvectors[1][2]; eigenvec2[2] = eigenvectors[2][2];
        // find the colinear vector to the normal vector at vertex_i, the cross product should be close to zero if the eigenvector is parallel to N_vi
        vtkMath::Cross(eigenvec0,normalVectorToVertex,c0);
        vtkMath::Cross(eigenvec1,normalVectorToVertex,c1);
        vtkMath::Cross(eigenvec2,normalVectorToVertex,c2);
        if(vtkMath::Norm(c0)<Tolerance)
        {
            lambda1 = eigenvalues[1];
            lambda2 = eigenvalues[2];
            T_1[0] = eigenvec1[0]; T_1[1]=eigenvec1[1]; T_1[2]=eigenvec1[2];
            T_2[0] = eigenvec2[0]; T_2[1]=eigenvec2[1]; T_2[2]=eigenvec2[2];
        }
        else if(vtkMath::Norm(c1)<Tolerance)
        {
            lambda1 = eigenvalues[0];
            lambda2 = eigenvalues[2];
            T_1[0] = eigenvec0[0]; T_1[1]=eigenvec0[1]; T_1[2]=eigenvec0[2];
            T_2[0] = eigenvec2[0]; T_2[1]=eigenvec2[1]; T_2[2]=eigenvec2[2];
        }
        else if(vtkMath::Norm(c2)<Tolerance)
        {
            lambda1 = eigenvalues[0];
            lambda2 = eigenvalues[1];
            T_1[0] = eigenvec0[0]; T_1[1]=eigenvec0[1]; T_1[2]=eigenvec0[2];
            T_2[0] = eigenvec1[0]; T_2[1]=eigenvec1[1]; T_2[2]=eigenvec1[2];
        }
        else
        {
            cout << "(!) No eigenvector parallel to N_vi at point i = : "<< i << " ." << endl;
        }
        // define the principal curvatures as functions of the non zero eigenvalues of M_vi
        K_1 = (3.0*lambda1)-lambda2;
        K_2 = (3.0*lambda2)-lambda1;

        // the curvature tensor is the symmetric tensor with eigenvectors: {N T_1 T_2} and eigenvalues {0 K_1 K_2}
        S_i[0][0] = K_1*T_1[0]*T_1[0] + K_2*T_2[0]*T_2[0];
        S_i[1][1] = K_1*T_1[1]*T_1[1] + K_2*T_2[1]*T_2[1];
        S_i[2][2] = K_1*T_1[2]*T_1[2] + K_2*T_2[2]*T_2[2];

        S_i[0][1] = K_1*T_1[0]*T_1[1] + K_2*T_2[0]*T_2[1];
        S_i[1][0] = S_i[0][1];

        S_i[0][2] = K_1*T_1[0]*T_1[2] + K_2*T_2[0]*T_2[2];
        S_i[2][0] = S_i[0][2];

        S_i[1][2] = K_1*T_1[1]*T_1[2] + K_2*T_2[1]*T_2[2];
        S_i[2][1] = S_i[1][2];

        curvatureTensor->SetTuple9(i,S_i[0][0],S_i[0][1],S_i[0][2],S_i[1][0],S_i[1][1],S_i[1][2],	S_i[2][0],S_i[2][1],S_i[2][2]);
        maxCurvature->SetTuple3(i,T_1[0],T_1[1],T_1[2]);
        minCurvature->SetTuple3(i,T_2[0],T_2[1],T_2[2]);
        principalCurvatures->SetTuple2(i,K_1,K_2);


    }// loop over all verticies

    input->GetPointData()->SetTensors(curvatureTensor);
    input->GetPointData()->SetVectors(maxCurvature);
    input->GetPointData()->SetVectors(minCurvature);
    input->GetPointData()->SetScalars(principalCurvatures);


    // cleaning
    if (principalCurvatures)    { principalCurvatures->Delete(); }
    if (maxCurvature)           { maxCurvature->Delete(); }
    if (minCurvature)           { minCurvature->Delete(); }
    if (curvatureTensor)        { curvatureTensor->Delete();}
    if (eigenvalues || M_vi || eigenvectors)
    {
        delete [] eigenvalues; eigenvalues = 0;
        delete [] M_vi; M_vi = 0;
        delete [] eigenvectors; eigenvectors = 0;
    }

}

//-----------------------------------------------------------------------------------------------------------------
void btkCurvatures::GetGaussCurvature(vtkPolyData * input)
{
    vtkDebugMacro("Start btkCurvatures::GetGaussCurvature()");
    // Data Array initialization
    this->GetCurvatureTensor(input);
    int numberOfVerticies = input->GetNumberOfPoints();
    vtkDoubleArray* gaussCurvature = vtkDoubleArray::New();
    gaussCurvature->SetName("Gauss_Curvature");
    gaussCurvature->SetNumberOfComponents(1);
    gaussCurvature->SetNumberOfTuples(numberOfVerticies);
    double *gaussCurvatureData = gaussCurvature->GetPointer(0);
    vtkDoubleArray *principalCurvatures = vtkDoubleArray::New();
    principalCurvatures = static_cast<vtkDoubleArray*>(input->GetPointData()->GetArray("Principal Curvatures"));
    double K_1 = 0.0;
    double K_2 = 0.0;

    // compute gauss curvature H = K_1*K_2
    for(vtkIdType i=0; i<numberOfVerticies; i++)
    {
        K_1 = principalCurvatures->GetComponent(i,0);
        K_2 = principalCurvatures->GetComponent(i,1);
        gaussCurvatureData[i] = K_1*K_2;
    }

    input->GetPointData()->AddArray(gaussCurvature);
    input->GetPointData()->SetActiveScalars("Gauss_Curvature");

    // clean
    if(gaussCurvature)    { gaussCurvature->Delete();}

}

//---------------------------------------------------------------------------------------------------------------------
void btkCurvatures::GetMeanCurvature(vtkPolyData *input)
{
    vtkDebugMacro("Start btkCurvatures::GetMeanCurvature()");
    // Data Array initialization
    this->GetCurvatureTensor(input);
    int numberOfVerticies = input->GetNumberOfPoints();
    vtkDoubleArray* meanCurvature = vtkDoubleArray::New();
    meanCurvature->SetName("Mean_Curvature");
    meanCurvature->SetNumberOfComponents(1);
    meanCurvature->SetNumberOfTuples(numberOfVerticies);
    double *meanCurvatureData = meanCurvature->GetPointer(0);
    vtkDoubleArray *principalCurvatures = static_cast<vtkDoubleArray*>(input->GetPointData()->GetArray("Principal Curvatures"));

    double K_1 = 0.0;
    double K_2 = 0.0;
    // compute mean curvature H = 1/2(K_1+K_2)
    for(vtkIdType i=0; i<numberOfVerticies; i++)
    {
        K_1 = principalCurvatures->GetComponent(i,0);
        K_2 = principalCurvatures->GetComponent(i,1);
        meanCurvatureData[i] = 0.5*(K_1+K_2);
    }
    input->GetPointData()->AddArray(meanCurvature);
    input->GetPointData()->SetActiveScalars("Mean_Curvature");

    // clean
    if(meanCurvature)    { meanCurvature->Delete();}
}

//---------------------------------------------------------------------------------------------------------------------
void btkCurvatures::GetBarCurvature(vtkPolyData *input)
{
    vtkDebugMacro("Start btkCurvatures::GetBarCurvature()");
    // Data Array initialization
    this->GetCurvatureTensor(input);
//    input->BuildLinks();
    int numberOfVerticies = input->GetNumberOfPoints();
    double barycenterOfNeighbourhood[3];
    double vectorIG[3];
    // Data used for computing neighbourhood of each vertex i
    vtkSmartPointer<vtkIdList> neighboursList = vtkSmartPointer<vtkIdList>::New();
    unsigned short numberOfCellsContainingVertex;
    vtkIdType* cellIds, *points, numberOfPointsInCell;
    double vertex_i[3], vertex_j[3];
    vtkDataArray *normals = input->GetPointData()->GetNormals();
    double normalVectorToVertex[3];

    vtkDoubleArray* barCurvature = vtkDoubleArray::New();
    barCurvature->SetName("Bar_Curvature");
    barCurvature->SetNumberOfComponents(1);
    barCurvature->SetNumberOfTuples(numberOfVerticies);
    double *barCurvatureData = barCurvature->GetPointer(0);

    // compute bar curvature B = Projection of the vector IG onto the normal to vertex_i (G = barycenter of local point neighbours)
    for(int i=0; i<numberOfVerticies; i++)
    {
        input->GetPointCells(i,numberOfCellsContainingVertex, cellIds); // get all cells containing i
        input->GetPoint(i,vertex_i);
        normals->GetTuple(i,normalVectorToVertex);
        // initialize barycenter computation
        barycenterOfNeighbourhood[0] = 0.0;
        barycenterOfNeighbourhood[1] = 0.0;
        barycenterOfNeighbourhood[2] = 0.0;

        // Reset the neighbours list and the Matrix M_vi for each vertex i
        neighboursList->Reset();
        // find the list of neighbours of vertex_i
        for(int c=0; c<numberOfCellsContainingVertex; c++)
        {
            input->GetCellPoints(cellIds[c], numberOfPointsInCell, points);
            for(int j=0; j<numberOfPointsInCell; j++)
            {
                if(points[j] != i)
                {
                    neighboursList->InsertNextId(points[j]);
                }
            }
        }
        // compute the barycenter of these points
        for(int j=0; j<neighboursList->GetNumberOfIds(); j++)
        {
            input->GetPoint(neighboursList->GetId(j),vertex_j);
            barycenterOfNeighbourhood[0] += vertex_j[0];
            barycenterOfNeighbourhood[1] += vertex_j[1];
            barycenterOfNeighbourhood[2] += vertex_j[2];
        }
        barycenterOfNeighbourhood[0] = barycenterOfNeighbourhood[0]/neighboursList->GetNumberOfIds();
        barycenterOfNeighbourhood[1] = barycenterOfNeighbourhood[1]/neighboursList->GetNumberOfIds();
        barycenterOfNeighbourhood[2] = barycenterOfNeighbourhood[2]/neighboursList->GetNumberOfIds();

        // compute the vector IG coordinates
        vectorIG[0] = barycenterOfNeighbourhood[0]-vertex_i[0];
        vectorIG[1] = barycenterOfNeighbourhood[1]-vertex_i[1];
        vectorIG[2] = barycenterOfNeighbourhood[2]-vertex_i[2];

        // projection of IG on Normal to vertex_i
        barCurvatureData[i] = vtkMath::Dot(vectorIG,normalVectorToVertex);

    }

    input->GetPointData()->AddArray(barCurvature);
    input->GetPointData()->SetActiveScalars("Bar_Curvature");

    if(barCurvature)    {barCurvature->Delete();}
}

//--------------------------------------------------------------------------------------------------------
int btkCurvatures::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Null input check
  if (!input)
    {
    return 0;
    }

  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetFieldData()->PassData(input->GetFieldData());

  //-------------------------------------------------------//
  //    Set Curvatures as PointData  Scalars               //
  //-------------------------------------------------------//

  if ( this->CurvatureType == CURVATURE_TENSOR )
  {
      this->GetCurvatureTensor(output);
  }
  else if ( this->CurvatureType == CURVATURE_GAUSS )
  {
      this->GetGaussCurvature(output);
  }
  else if ( this->CurvatureType == CURVATURE_MEAN )
  {
      this->GetMeanCurvature(output);
  }
  else if ( this->CurvatureType ==  CURVATURE_BAR )
  {
      this->GetBarCurvature(output);
  }
  else
  {
      vtkErrorMacro("Only Gauss, Mean, and Bar Curvature type available");
      return 1;
  }

  return 1;
}

/*-------------------------------------------------------*/
void btkCurvatures::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "CurvatureType: " << this->CurvatureType << "\n";
}

