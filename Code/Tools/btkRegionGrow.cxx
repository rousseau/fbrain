/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/07/2013
  Author(s): Aïcha Bentaieb (abentaieb@unistra.fr)

  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.

==========================================================================*/
#include "btkRegionGrow.h"


//*************************************
//AICHA : ces fichiers devraient etre dans btkRegionGrow.h ???
//*************************************
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataReader.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataNormals.h"
#include "vtkMath.h"

vtkCxxRevisionMacro(btkRegionGrow, "$Revision: 1.15 $");
vtkStandardNewMacro(btkRegionGrow);

btkRegionGrow::btkRegionGrow()
{
    m_Size = 0;
    // init
}

//*************************************
//AICHA : COMMENTER CETTE FONCTION ?
//*************************************
btkRegionGrow::IdVector btkRegionGrow::getNeighbours(vtkIdType point)
{
    IdVector  neighboursList;
    unsigned short numberOfCellsContainingVertex;
    vtkIdType* cellIds, *points, numberOfPointsInCell;

    this->m_Polydata->GetPointCells(point,numberOfCellsContainingVertex, cellIds); // get all cells containing i
    for(int c=0; c<numberOfCellsContainingVertex; c++)
    {
        this->m_Polydata->GetCellPoints(cellIds[c], numberOfPointsInCell, points);
        for(int j=0; j<numberOfPointsInCell; j++)
        {
            if( (points[j] != point) && (this->isInside(points[j], neighboursList) == false) )
            {
                neighboursList.push_back(points[j]);
            }
        }
    }
    return neighboursList;
}

// ---------------------------------------------------------------------------------------------------------------------
std::vector<double> btkRegionGrow::getNormalVector(vtkIdType point)
{
    vtkDataArray *normals = this->m_Polydata->GetPointData()->GetNormals();
    double normalVector[3];
    normals->GetTuple(point,normalVector);
    std::vector<double> NormalToPoint;
    NormalToPoint.push_back(normalVector[0]);
    NormalToPoint.push_back(normalVector[1]);
    NormalToPoint.push_back(normalVector[2]);

    return NormalToPoint;

}

// ---------------------------------------------------------------------------------------------------------------------
double btkRegionGrow::getCurvValue(vtkIdType pointId)
{
    vtkDoubleArray *curvArray = vtkDoubleArray::New();
    curvArray = static_cast<vtkDoubleArray*>(this->m_Polydata->GetPointData()->GetArray("Bar_Curvature"));
    double curv_I = curvArray->GetComponent(pointId,0);

    return curv_I;

}
// ---------------------------------------------------------------------------------------------------------------------
//*************************************
//AICHA : A QUOI SERT CETTE FONCTION ?
//*************************************

bool btkRegionGrow::isInside(vtkIdType value, IdVector vec)
{
    for(int i = 0; i<vec.size(); i++)
    {
        if(vec[i] == value)
        {
            return true;
        }
        else
            continue;
    }
    return false;
}

// ---------------------------------------------------------------------------------------------------------------------
void btkRegionGrow::growing(vtkIdType m_SeedPoint , int m_numberOfPoints)
{

    // initialize
    vtkIdType seedPoint = m_SeedPoint;
    m_Vec.push_back(seedPoint);
    m_Size++;
    double curv_i = 0.0;
    double curv = getCurvValue(seedPoint);
    double vec1[3];
    double vec2[3];

    // get normal vector to seed point
    std::vector<double> normalSP;
    normalSP = this->getNormalVector(seedPoint);
    vec1[0] = normalSP[0]; vec1[1]=normalSP[1]; vec1[2]=normalSP[2];

    // get the neighbours
    IdVector voisins;
    voisins = this->getNeighbours(seedPoint);

    std::vector<double> normalVertex;


    // region grow
    unsigned int it = 0;

    while( it<100 && m_Size<m_numberOfPoints )
    {

        for(int i =0; i<voisins.size(); i++)
        {
        
//*************************************
//AICHA : COMMENTER CETTE BOUCLE?
//*************************************
            vtkIdType currentNeighbour = voisins[i];

            curv_i = getCurvValue(currentNeighbour);

            normalVertex = this->getNormalVector(currentNeighbour);

            vec2[0]=normalVertex[0]; vec2[1]=normalVertex[1]; vec2[2]=normalVertex[2];

            normalVertex.clear();
            normalSP.clear();

            double prod = vtkMath::Dot(vec1,vec2)/(vtkMath::Norm(vec1) * vtkMath::Norm(vec2));
            //double angle = std::acos(prod)*(180/PI); // add an other criterion based on the angle between normals of points


            if( (std::fabs(curv_i) <= std::fabs(curv)) && (this->isInside(currentNeighbour, m_Vec) == false) /*&& (angle < 15)*/)
            {
                m_SeedPoint = currentNeighbour;
                growing(m_SeedPoint, m_numberOfPoints);
            }

            if(m_Size == m_numberOfPoints)
                break;
        it++;

        }

    }

}

// ---------------------------------------------------------------------------------------------------------------------
void btkRegionGrow::Update()
{
    this->growing(m_SeedPoint, m_numberOfPoints);
}
