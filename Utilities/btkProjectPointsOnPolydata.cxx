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


/* input1 = polydata to project points on
 * input2 = points selected after region grow
 * output = corresponding projection points
 */

#include <tclap/CmdLine.h>
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkMath.h"
#include "vtkVertex.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkActor.h"
#include "vtkIdTypeArray.h"
#include "vtkPointData.h"

// ---------------------------------------------------------------------------------------------------------------------

// function checking wheter a value is inside a vector of IDs
bool isInside(vtkIdType value, std::vector<vtkIdType> vec)
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

// Function finding the neighbours of a point
std::vector<vtkIdType> getNeighbours(vtkPolyData *Polydata, vtkIdType point)
{
    std::vector<vtkIdType>  neighboursList;
    unsigned short numberOfCellsContainingVertex;
    vtkIdType* cellIds, *points, numberOfPointsInCell;

    Polydata->GetPointCells(point,numberOfCellsContainingVertex, cellIds); // get all cells containing i
    for(int c=0; c<numberOfCellsContainingVertex; c++)
    {
        Polydata->GetCellPoints(cellIds[c], numberOfPointsInCell, points);
        for(int j=0; j<numberOfPointsInCell; j++)
        {
            if( (points[j] != point) && (isInside(points[j], neighboursList) == false) )
            {
                neighboursList.push_back(points[j]);
            }
        }
    }
    return neighboursList;
}
// ---------------------------------------------------------------------------------------------------------------------

// Main function of projection according to the nearest neighbour
int main(int argc, char *argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Project points on a polydata",' ',"v1.0",true);
        TCLAP::ValueArg<std::string> inputPolydata("i","input_poly","vtk polydata to project points on",true,"","string");
        cmd.add(inputPolydata);
        TCLAP::ValueArg<std::string> inputPoints("m","input_points","vtk points",true,"","string");
        cmd.add(inputPoints);
        TCLAP::ValueArg<std::string> output("o","output_poly","vtk output highlighted projection points",true,"","string");
        cmd.add(output);
        cmd.parse(argc,argv);
        std::string input_points = inputPoints.getValue();
        std::string input_poly = inputPolydata.getValue();
        std::string output_poly = output.getValue();

        //Read inputs
        vtkSmartPointer<vtkPolyDataReader> readerPoints = vtkSmartPointer<vtkPolyDataReader>::New();
        readerPoints->SetFileName(input_points.c_str());
        readerPoints->Update();
        vtkSmartPointer<vtkPolyDataReader> readerPd = vtkSmartPointer<vtkPolyDataReader>::New();
        readerPd->SetFileName(input_poly.c_str());
        readerPd->Update();

        //data
        vtkSmartPointer<vtkPolyData> hull = vtkSmartPointer<vtkPolyData>::New();
        hull = readerPd->GetOutput();
        hull->Update();
        vtkSmartPointer<vtkPolyData> points = vtkSmartPointer<vtkPolyData>::New();
        points = readerPoints->GetOutput();
        points->GetPoints();
        points->Update();

        // Variables declaration
        hull->BuildLinks();
        points->BuildLinks();
        int N = points->GetNumberOfPoints();

        double dist[N][3];
        double projection[3];
        double initial[3];
        std::vector<double> distance;
        std::vector<vtkIdType> closestPointId;

        for(int i = 0; i<N; i++)
        {
            points->GetPoint(i,initial);// get each point to project coordinates
            for(int j=0; j<hull->GetNumberOfPoints(); j++)
            {
                hull->GetPoint(j, projection);
                distance.push_back(vtkMath::Distance2BetweenPoints(initial,projection)); // compute the distance between points of the surface and the point to project
            }

            std::vector<double>::iterator result = std::min_element(distance.begin(),distance.end());
            int d = std::distance(distance.begin(), result);

            closestPointId.push_back(d);

            // store the corresponding projection point
            if (isInside(d, closestPointId)==true)
            {

                std::vector<vtkIdType> voisins;
                voisins = getNeighbours(hull, d);
                d = voisins[1];
                closestPointId.pop_back();
                closestPointId.push_back(d);
            }

            distance.clear();
        }


        for(int i=0; i<closestPointId.size();i++)
        {
            hull->GetPoint(i,projection);
            dist[i][0] = projection[0];
            dist[i][1] = projection[1];
            dist[i][2] = projection[2];
        }



        vtkSmartPointer<vtkPoints> projectionPoints = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
        vtkSmartPointer<vtkCellArray> verticies = vtkSmartPointer<vtkCellArray>::New();// topology
        // store the coordinates of the projections
        for(int i=0; i<N; i++)
        {
            projectionPoints->InsertNextPoint(dist[i][0],dist[i][1],dist[i][2]);
        }
        // find the verticies corresponding to the projections
        vertex->GetPointIds()->SetNumberOfIds(projectionPoints->GetNumberOfPoints());
        for(int i=0; i<projectionPoints->GetNumberOfPoints(); i++)
        {
            vertex->GetPointIds()->SetId(i,i);
        }
        verticies->InsertNextCell(vertex);

        vtkSmartPointer<vtkPolyData> outputPD = vtkSmartPointer<vtkPolyData>::New();
        outputPD->SetPoints(projectionPoints);
        outputPD->SetVerts(verticies);

        // write the output vtk file
        vtkSmartPointer<vtkPolyDataWriter> writer =  vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(output_poly.c_str());
        writer->SetInput(outputPD);
        writer->Write();


        return EXIT_SUCCESS;

    }catch (TCLAP::ArgException &e) // catch any exceptions
    {std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}//done
