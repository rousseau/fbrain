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

/* Includes */
#include <iostream>
#include <string>
#include <fstream>
#include <tclap/CmdLine.h>
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"

int main(int argc, char *argv[])
{
    try
    {

        TCLAP::CmdLine cmd("txt file of scalar curv values.", ' ', "1.0");
        TCLAP::ValueArg<std::string> input("i", "input_file", "Input mesh curvature file (.vtk) ",true," ","string" );
        TCLAP::ValueArg<std::string> inputPoints("m", "input_points", "points of interest to project on the surface (.vtk)",true," ","string" );
        TCLAP::ValueArg<std::string> outputFile("o", "output", "points of interest corresponding curv (.txt)",true," ","string" );


        cmd.add(input);
        cmd.add(inputPoints);
        cmd.add(outputFile);
        cmd.parse(argc,argv);

        std::string input_file = input.getValue();
        std::string input_points = inputPoints.getValue();
        std::string output= outputFile.getValue();


        // read file
        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(input_file.c_str());
        reader->Update();
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata = reader->GetOutput();
        polydata->Update();

        // get curv array
        vtkDoubleArray *curvArray = vtkDoubleArray::New();
        curvArray = static_cast<vtkDoubleArray*>(polydata->GetPointData()->GetArray("Bar_Curvature"));
        double curv_I = 0.0;

        // get points of interest (file)
        vtkSmartPointer<vtkPolyDataReader> readerpoints = vtkSmartPointer<vtkPolyDataReader>::New();
        readerpoints->SetFileName(input_points.c_str());
        readerpoints->Update();

        vtkSmartPointer<vtkPolyData> polydataPoints = vtkSmartPointer<vtkPolyData>::New();
        polydataPoints = readerpoints->GetOutput();
        polydataPoints->Update();

        // open txt output file
        std::ofstream fichier(output.c_str(), std::ios::trunc);

        if(fichier)
        {
            // get curv array for each point
            for(int i=0; i<polydataPoints->GetNumberOfPoints(); i++)
            {
                curv_I = curvArray->GetComponent(i,0);
                fichier << i <<" , "<< curv_I << std::endl;
            }
            fichier.close();
        }
        else
                std::cerr << "ERROR, impossible to open file !" << std::endl;

        return EXIT_SUCCESS;

    }catch (TCLAP::ArgException &e) // catch any exception
    {std::cerr << "error: " << e.error() << "for arg " << e.argId() << std::endl;}
}

