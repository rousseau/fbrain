/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 10/01/2013 (last updated: 23/05/2013)
  Author(s): Larbi Boubchir (boubchir at unistra dot fr)
  
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

#include <sstream>

// VTK includes
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkParticleReader.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkLine.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyLine.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>

// ---------------------------------
// Usage : btkGenerateVtkFileFromFiberDataTextFiles [number_of_fibers]
//
// 	Input argument:
// 	[number_of_fibers]: the number of fibers where each fiber data is stored into a text file (.txt) containing the coordinates X Y Z of all the points of a fiber
//
// 	Output:
// 	Generate and visualize a vtk file containing all fibers
//
// Example : btkGenerateVtkFileFromFiberDataTextFiles 16
// This code can be tested to generate the vtk file of 'Fiber Cup phantom' which contains 16 fibers 
// (16 files {Fiber1.txt, Fiber2.txt, ... , Fiber16.txt}).
// Fiber Cup phantom data can be found here: http://www.lnao.fr/spip.php?rubrique79
// More information can be found in [P. Fillard et al.,"Quantitative evaluation of 10 tractography algorithms on a realistic diffusion MR phantom", NeuroImage, pp. 220–234, Vol. 56(1), 2011.]
// the output data is generated into a vtk file nammely 'All_Fibers.vtk' 
// 
// ---------------------------------

//
// main
int main(int argc, char* argv[])
{

// Verify input arguments
if ( argc != 2 )
   {
    std::cout << "\nUsage: " << argv[0]
              << "btkGenerateVtkFileFromFiberDataTextFiles [number_of_fibers]\n" << std::endl;
    return EXIT_FAILURE;
   }
   
// Get all data from the file
int filenumber = atoi( argv[1] );
std::cout << "Number of fibers: " << filenumber << std::endl;

vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New(); 

uint nb=0;
for(int f=1;f<=filenumber;f++)
{
  std::ostringstream oss;
  oss << "Fiber" << f << ".txt"; // Each fiber data is stored into a text file namely Fiber%.txt where % corresponds to 'filenumber'.
  std::string filename = oss.str();
  std::ifstream fin(filename.c_str());
 
  std::string line;
  
  uint nb_points = 0;
  while(std::getline(fin, line))
    {
    double x,y,z;
    std::stringstream linestream;
    linestream << line;
    linestream >> x >> y >> z; // import the coordinates X Y Z of each point
    points->InsertNextPoint(x, y, z);
    nb_points++; 
    }
  fin.close();
   
  polyLine->GetPointIds()->SetNumberOfIds(nb_points);
   
  for(unsigned int i = 0; i < nb_points; i++)
    {
      //std::cout << i+(f-1)*nb_points << std::endl;
      polyLine->GetPointIds()->SetId(i,nb); nb++;
    }

  cells->InsertNextCell(polyLine);
 }
  // Add the points to the dataset
  polydata->SetPoints(points);
  // Add the lines to the dataset
  polydata->SetLines(cells);
 
  //polydata->BuildLinks();

  // Save output data into VTK file
  vtkSmartPointer<vtkPolyDataWriter> writer = 
    vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetInput(polydata);
  
  writer->SetFileName("All_Fibers.vtk"); // save the output data
  writer->SetFileTypeToASCII();
  writer->Write();

// ----------------------------------- 
// Visualization
// ----------------------------------- 
   vtkSmartPointer<vtkVertexGlyphFilter> glyphFilter =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  glyphFilter->SetInputConnection(polydata->GetProducerPort());
#else
  glyphFilter->SetInputData(polydata);
#endif
  glyphFilter->Update();
 
  // Visualize
  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(glyphFilter->GetOutputPort());
 
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
 
  //Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);
 
  //Add the actor to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(.3, .3, .3); // Background color green
 
  //Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}
