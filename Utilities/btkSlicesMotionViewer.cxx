/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 08/03/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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
/* ITK */
#include "itkImage.h"

/* BTK */
#include "btkImageHelper.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkIOTransformHelper.h"
#include "btkSlicesToPolydata.h"
#include "btkRenderPolydata.h"
/* VTK */
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"


/* OTHERS */
#include "iostream"
#include "sstream"
#include <tclap/CmdLine.h>


const unsigned int Dimension = 3;
typedef float PixelType;
//typedef  double PixelType;
typedef itk::Image< PixelType, Dimension > itkImage;
typedef btk::EulerSliceBySliceTransform< double, Dimension, PixelType > Transform;
typedef btk::SlicesToPolyData< itkImage> PolyDataFilter;
typedef btk::RenderPolyData RenderP;

int main( int argc, char *argv[] )
{
    TCLAP::CmdLine cmd("Construct a polyData which outline the slices for each input images", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","input image ",true,"string",cmd);
    TCLAP::MultiArg<std::string> tranArg("t","transform","transforms to apply",true,"string",cmd);
    TCLAP::MultiArg<std::string> outputArg("o","output","vtk polydata",true,"string",cmd);
   TCLAP::SwitchArg  RenderArg("r","render","Open a render window with the polydatas", cmd, false);


    std::vector< std::string > input;
    std::vector< std::string > outputName;
    std::vector<std::string> transfoNames;

    std::vector< itkImage::Pointer > inputsImages;

    // Parse the argv array.
    cmd.parse( argc, argv );
    input = inputArg.getValue();
    outputName = outputArg.getValue();
    transfoNames = tranArg.getValue();
    bool render = RenderArg.getValue();

    std::vector< Transform::Pointer > transforms;
    transforms.resize(inputsImages.size());

    inputsImages = btk::ImageHelper<itkImage>::ReadImage(input);
    transforms =  btk::IOTransformHelper< Transform >::ReadTransform(transfoNames);
    std::vector< vtkSmartPointer< vtkPolyData> > output(inputsImages.size());

    RenderP::Pointer renderer = RenderP::New();
    renderer->SetNumberOfPolyData(inputsImages.size());

    for(unsigned int i = 0; i< inputsImages.size(); i++)
    {
        output[i] = vtkSmartPointer< vtkPolyData >::New();
        PolyDataFilter::Pointer filter = PolyDataFilter::New();
        transforms[i]->SetImage(inputsImages[i]);
        filter->SetInput(inputsImages[i]);
        filter->SetTransform(transforms[i]);
        filter->Update();
        output[i] = filter->GetOutput();
        std::vector<double> color(3);
        color[0] =color[1]=color[2]= 0.0;
        color[i] = 255.0;
        renderer->SetNthPolyData(i,output[i]);
        renderer->SetNthPolyDataColor(i,color);

        vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
        writer->SetInput(output[i]);
        writer->SetFileName(outputName[i].c_str());
        writer->Write();

    }
    if(render)
    {
        renderer->Render();
    }
}
