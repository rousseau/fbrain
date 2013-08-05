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


/* filter to select a seed point with mouse click
 * using vtkInteractorStyle
 */

#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCommand.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPointPicker.h"
#include "vtkObjectFactory.h"
#include "vtkRendererCollection.h"
#include "vtkProp3DCollection.h"
#include "vtkCallbackCommand.h"
#include "vtkWorldPointPicker.h"
#include "vtkDataSetMapper.h"
#include "vtkCellPicker.h"
#include "vtkIdTypeArray.h"
#include "vtkSelectionNode.h"
#include "vtkSelection.h"
#include "vtkExtractSelection.h"
#include "vtkUnstructuredGrid.h"
#include "vtkProperty.h"

#include <stdio.h>
#include <tclap/CmdLine.h>

// Catch mouse events
class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
  public:
  static MouseInteractorStyle* New();

  MouseInteractorStyle()
  {
    selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
    selectedActor = vtkSmartPointer<vtkActor>::New();
  }

    virtual void OnLeftButtonDown()
    {
      // Get the location of the click (in window coordinates)
      int* pos = this->GetInteractor()->GetEventPosition();

      vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
      picker->SetTolerance(0.0005);

      // Pick from this location.
      picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

      double* worldPosition = picker->GetPickPosition();
      std::cout << "Cell id is: " << picker->GetCellId() << std::endl;

      if(picker->GetCellId() != -1)
        {

            std::cout << "Pick position is: " << worldPosition[0] << " " << worldPosition[1]
                  << " " << worldPosition[2] << endl;

            vtkSmartPointer<vtkIdTypeArray> ids = vtkSmartPointer<vtkIdTypeArray>::New();
            ids->SetNumberOfComponents(1);
            ids->InsertNextValue(picker->GetCellId());

            vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
            selectionNode->SetFieldType(vtkSelectionNode::CELL);
            selectionNode->SetContentType(vtkSelectionNode::INDICES);
            selectionNode->SetSelectionList(ids);

            vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
            selection->AddNode(selectionNode);

            vtkSmartPointer<vtkExtractSelection> extractSelection = vtkSmartPointer<vtkExtractSelection>::New();
            extractSelection->SetInput(0, this->Data);
            extractSelection->SetInput(1, selection);
            extractSelection->Update();

            // In selection
            vtkSmartPointer<vtkUnstructuredGrid> selected =  vtkSmartPointer<vtkUnstructuredGrid>::New();
            selected->ShallowCopy(extractSelection->GetOutput());

            std::cout << "There are " << selected->GetNumberOfPoints()
                      << " points in the selection." << std::endl;
            std::cout << "There are " << selected->GetNumberOfCells()
                      << " cells in the selection." << std::endl;


            selectedMapper->SetInputConnection(selected->GetProducerPort());

            selectedActor->SetMapper(selectedMapper);
            selectedActor->GetProperty()->EdgeVisibilityOn();
            selectedActor->GetProperty()->SetEdgeColor(1,0,0);
            selectedActor->GetProperty()->SetLineWidth(3);

            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);

        }
      // Forward events
      vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }

    vtkSmartPointer<vtkPolyData> Data;
    vtkSmartPointer<vtkDataSetMapper> selectedMapper;
    vtkSmartPointer<vtkActor> selectedActor;

};
vtkStandardNewMacro(MouseInteractorStyle);


// Main function
int main(int argc, char *argv[])
{
    try
    {
        /* 1. input image = triangle mesh polydata */
        TCLAP::CmdLine cmd("Select a cell using mouse clic on a visual rendering", ' ', "1.0");
        TCLAP::ValueArg<std::string> inputFile("i", "input_file", "Input triangulated mesh file (vtk file)", true, "", "string");
        cmd.add(inputFile);
        cmd.parse(argc,argv);
        std::string input_file = inputFile.getValue();


        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(input_file.c_str());
        reader->Update();
        vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
        input = reader->GetOutput();
        input->Update();


		//A QUOI SERT UN MAPPER, ETC ?
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInput(input);
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
        renderer->AddActor(actor);
        renderer->ResetCamera();

        vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();
        renwin->AddRenderer(renderer);

        vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        iren->SetRenderWindow(renwin);

        vtkSmartPointer<MouseInteractorStyle> style = vtkSmartPointer<MouseInteractorStyle>::New();
        style->SetDefaultRenderer(renderer);
        style->Data = input;


        iren->SetInteractorStyle(style);

        renwin->Render();
        iren->Start();
		// COMMENTER CHAQUE LIGNE DE CETTE SECTION (DEPUIS MAPPER...)


        return EXIT_SUCCESS;

    }catch (TCLAP::ArgException &e) // catch any exception
    {std::cerr << "error: " << e.error() << "for arg " << e.argId() << std::endl;}

}

