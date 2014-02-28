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

// STL includes
#include "cmath"

// VTK includes used in functions declared below
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

// ---------------------------------------------------------------------------------------------------------------------
// function used to find the neighbours of any point from its ID
btkRegionGrow::IdVector btkRegionGrow::getNeighbours(vtkIdType point)
{
    IdVector  neighboursList; // store the neighbours in a vector of Ids
    unsigned short numberOfCellsContainingVertex;
    vtkIdType* cellIds, *points, numberOfPointsInCell;

    this->m_Polydata->GetPointCells(point,numberOfCellsContainingVertex, cellIds); // get all cells containing i

    for(int c=0; c<numberOfCellsContainingVertex; c++)
    {
        this->m_Polydata->GetCellPoints(cellIds[c], numberOfPointsInCell, points); // get all the points present in the current cell
        for(int j=0; j<numberOfPointsInCell; j++)
        {
            // if the point find in the current cell isn't the considered initial point, we store its ID in the neighbours list
            if( (points[j] != point) && (this->isInside(points[j], neighboursList) == false) ) // check if the point wasn't already added to the list as points can appear many times in different cells
            {
                neighboursList.push_back(points[j]);
            }
        }
    }
    return neighboursList;
}

// ---------------------------------------------------------------------------------------------------------------------
// get the normal vector of a point by giving its ID
std::vector<double> btkRegionGrow::getNormalVector(vtkIdType point)
{
    vtkDataArray *normals = this->m_Polydata->GetPointData()->GetNormals(); // normals are stored in the input polydata as a data array
    double normalVector[3];
    normals->GetTuple(point,normalVector);
    std::vector<double> NormalToPoint; // store the vector normal
    NormalToPoint.push_back(normalVector[0]);
    NormalToPoint.push_back(normalVector[1]);
    NormalToPoint.push_back(normalVector[2]);

    return NormalToPoint;

}

// ---------------------------------------------------------------------------------------------------------------------
// get the curvature value of a point
double btkRegionGrow::getCurvValue(vtkIdType pointId)
{
    // the curvature value is storred in the polydata as a vtkDoubleArray (according to btkCurvatures.cxx)
    vtkDoubleArray *curvArray = vtkDoubleArray::New();
    curvArray = static_cast<vtkDoubleArray*>(this->m_Polydata->GetPointData()->GetArray("Bar_Curvature"));// copy the curvature array from the polydata
    double curv_I = curvArray->GetComponent(pointId,0);

    return curv_I;

}
// ---------------------------------------------------------------------------------------------------------------------
// Function that checks if a value of type ID is present inside a vector of IDs
bool btkRegionGrow::isInside(vtkIdType value, IdVector vec)
{
    for(unsigned int i = 0; i<vec.size(); i++) // go through all the elements of the vector to check if the value is present
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
void btkRegionGrow::growing(vtkIdType m_SeedPoint , unsigned int m_numberOfPoints)
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

    while( it<m_Nb_Loop && m_Size<m_numberOfPoints ) // while the size criterion and the number of itterations isnt reached we compute a region grow trying to find the neighbouring verticies of a seedPoint that have a smaller or equivalent curvature
    {

        for(unsigned int i =0; i<voisins.size(); i++) // for all the neighbours of the seedPoint we check the curvature value
        {
        
            vtkIdType currentNeighbour = voisins[i];

            curv_i = getCurvValue(currentNeighbour); // get the curvature of the neighbour point computed

            normalVertex = this->getNormalVector(currentNeighbour); // get the normal vector at this neighbour vertex

            vec2[0]=normalVertex[0]; vec2[1]=normalVertex[1]; vec2[2]=normalVertex[2]; // store the normal vector

            normalVertex.clear(); // reset the normal vector for the next point computed
            normalSP.clear();// reset the normal to the seedPoint vector for the next seedPoint

            /* adding a criterion based on the angle between every normals of the preselected points - ROI points could be a line of parallel normals to point vectors*/
            /* ADD THESE LINES IF ANGLE UNCOMMENTED */
            //double prod = vtkMath::Dot(vec1,vec2)/(vtkMath::Norm(vec1) * vtkMath::Norm(vec2));
            //double angle = std::acos(prod)*(180/PI); // add an other criterion based on the angle between normals of points


            // Check if the curvature value between the seedPoint and the neighbouring points is filling the criterion and decrease -> in that case add the value to the list of point IDs of the ROI
            if( (std::abs(curv_i) <= std::abs(curv)) && (this->isInside(currentNeighbour, m_Vec) == false) /*&& (angle < 15)*/)
            {
                m_SeedPoint = currentNeighbour;
                growing(m_SeedPoint, m_numberOfPoints); // compute the region grow from the new seedPoint (first neighbour selected)
            }

            if(m_Size == m_numberOfPoints)
                break;
        it++;

        }

    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Update function -> Updates the list of points selected after region growing
void btkRegionGrow::Update()
{
    this->growing(m_SeedPoint, m_numberOfPoints);
}
