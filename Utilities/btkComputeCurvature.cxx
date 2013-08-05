/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 26/06/2012
  Author(s): Aïcha BenTaieb (abentaieb@unistra.fr)

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
#include <tclap/CmdLine.h>

/* BTK */
#include "btkCurvatures.h"

/* VTK */
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"


int main(int argc, char *argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Compute curvature on a triangulated mesh.", ' ', "1.0");
        TCLAP::ValueArg<std::string> input("i", "input_file","Input VTK triangle mesh polydata.", true, " ", "string");
        cmd.add(input);
        TCLAP::ValueArg<std::string> output("o", "output_file", "Output VTK polydata (including curvature, tensors, normal, principal curvatures).", true, " ", "string");
        cmd.add(output);

        cmd.parse(argc,argv);

        std::string input_file = input.getValue();
        std::string output_file = output.getValue();

        // read input file (triangle + smoothed mesh)
        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(input_file.c_str());
        reader->Update();

        // get the polydata
        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata = reader->GetOutput();
        polydata->Update();

        btkCurvatures *curv = btkCurvatures::New();
        curv->SetInput(polydata);

        std::cout<<"Choose a curvature type: \n";
        std::cout<<"0: barycentric curvature \n";
        std::cout<<"1: gaussian curvature \n";
        std::cout<<"2: mean curvature \n";
        std::cout<<"3: principal curvatures."<<std::endl;

        int choice = -1;
        std::cin>>choice;

        switch(choice)
        {
        case 0:
            std::cout<<"Curvature type = BARYCENTER"<<std::endl;
            curv->SetCurvatureTypeToBar();
            break;
        case 1:
            std::cout<<"Curvature type = GAUSSIAN"<<std::endl;
            curv->SetCurvatureTypeToGaussian();
            break;
        case 2:
            std::cout<<"Curvature type = MEAN"<<std::endl;
            curv->SetCurvatureTypeToMean();
            break;
        case 3:
            std::cout<<"Curvature type = PRINCIPAL CURVATURES"<<std::endl;
            curv->SetCurvatureTypeToTensor();
        }

        curv->Update();

        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetInput(curv->GetOutput());
        writer->SetFileName(output_file.c_str());
        writer->Write();

        return EXIT_SUCCESS;

    }catch (TCLAP::ArgException &e) // catch any exception
    {std::cerr << "error: " << e.error() << "for arg " << e.argId() << std::endl;}
}

