/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 08/04/2013
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

/* This program takes as input any mesh type made of triangle strips or polygons and
 *converts it into a polyhedral mesh of triangles exported in output file.
**/

/* VTK Includes */
#include "vtkVersion.h"
#include "vtkPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkSmartPointer.h"
#include "vtkGenericDataObjectReader.h"
#include "vtkCellTypes.h"
#include "vtkPolyDataWriter.h"

/* OTHERS */
#include "iostream"
#include "map"
#include <tclap/CmdLine.h>


int main(int argc, char *argv[])
{
    try
    {
        // Command line arguments.
        TCLAP::CmdLine cmd("Transform any input mesh into a closed polyhedral mesh of triangles using vtkTriangleFilter.", ' ', "1.0");

        TCLAP::ValueArg<std::string> inputMeshArg("i", "input_file", "Input mesh file (vtk file)", true, "", "string");
        cmd.add(inputMeshArg);
        TCLAP::ValueArg<std::string> outputMeshArg("o", "output_file","Output vtk file (triangular mesh, i.e. one cell = one triangle)", true,"","string");
        cmd.add(outputMeshArg);

        // Parse Args.
        cmd.parse(argc,argv);

        // Get the value parsed by each arg.
        std::string input_file  = inputMeshArg.getValue();
        std::string output_file = outputMeshArg.getValue();

        // Read the input mesh file
        std::cout<< "Reading input file" << std::endl;
        vtkSmartPointer<vtkGenericDataObjectReader> fileReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        fileReader->SetFileName(input_file.c_str());
        fileReader->Update();
        // Get the polydata from the input file
        vtkSmartPointer<vtkPolyData> inputPolydata = vtkSmartPointer<vtkPolyData>::New();
        inputPolydata = fileReader->GetPolyDataOutput();
        inputPolydata->Update();
        if(fileReader->IsFilePolyData())
        {
            std::cout << "Input file is of type vtkPolydata" << std::endl;
            std::cout << "Input file has: " << inputPolydata->GetNumberOfPoints() << " points." << std::endl;
        }
        

        // Cell container used to store the different cell types in the input polydata
        typedef std::map<int,int> CellContainer;
        CellContainer cellMap;
        for (int i=0; i<inputPolydata->GetNumberOfCells(); i++)
        {
            cellMap[inputPolydata->GetCellType(i)]++; // store each cell type encountered in the input mesh in a cell Container
        }
        // Print the different cell types before applying conversion to triangle cells
        CellContainer::const_iterator it = cellMap.begin();
        while(it != cellMap.end())
        {
            std::cout << "Cell type = " << vtkCellTypes::GetClassNameFromTypeId(it->first) << " occurs " << it->second << " times."<< std::endl;
            ++it;
        }

        // Generate triangles from input
        std::cout << "Conversion of celltypes to triangles ...  " << std::endl;
        vtkSmartPointer<vtkTriangleFilter> triangleMesh = vtkSmartPointer<vtkTriangleFilter>::New();
        triangleMesh->SetInput(fileReader->GetOutput());
        triangleMesh->PassLinesOff();  // No individual line segments
        triangleMesh->PassVertsOff(); // No individual vertices
        triangleMesh->Update();

        // Store the mesh in a polydata structure
        std::cout << "Writting vtk polydata file ... " << std::endl;
        vtkSmartPointer<vtkPolyData> outputPolydata = vtkSmartPointer<vtkPolyData>::New();
        outputPolydata = triangleMesh->GetOutput();
        outputPolydata->Update();

        // Write output file as a vtk polydata file
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInput(outputPolydata);
        writer->SetFileName(output_file.c_str());
        writer->Write();

        // Check the type of cells created
        CellContainer outCellMap;
        for (int i=0; i<outputPolydata->GetNumberOfCells(); i++)
        {
            outCellMap[outputPolydata->GetCellType(i)]++;
        }
        CellContainer::const_iterator i = outCellMap.begin();
        while(i != outCellMap.end())
        {
            std::cout << "New Cell type = " << vtkCellTypes::GetClassNameFromTypeId(i->first) << " occurs " << i->second << " times."<< std::endl;
            ++i;
        }

        return 1;

    } catch (TCLAP::ArgException &e) // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
