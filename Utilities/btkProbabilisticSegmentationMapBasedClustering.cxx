/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 10/01/2013 (last updated: 28/05/2013)
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

// TCLAP : Templatized C++ Command Line Parser includes
#include <tclap/CmdLine.h>

// STL includes
#include "cstdlib"
#include "string"
#include "iomanip"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// VTK includes
#include "vtkGenericDataObjectReader.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkPointData.h"

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkNumericTraits.h"
#include "itkPoint.h"
#include "itkContinuousIndex.h"

typedef float PixelType;
const unsigned int Dimension = 4;
typedef itk::Image< PixelType,Dimension > Image4DType;
typedef Image4DType::Pointer Image4DPointer;
typedef itk::ImageFileReader< Image4DType >  Reader4DType;

typedef itk::Image<PixelType,Dimension-1> Image3DType;
typedef Image3DType::Pointer Image3DPointer;

// Usage   : btkProbabilisticSegmentationMapBasedClustering -b inputfile.vtk -p prob_map -o outputfile.vtk
// Example : btkProbabilisticSegmentationMapBasedClustering -b data.vtk -p 01019-natbrain4D.nii.gz -o clustering1-data.vtk

// main
int main ( int argc, char *argv[] )
{

  // Define command line parser
  TCLAP::CmdLine cmd("White-Matter Fiber Tracts Clustering based on a probabilistic segmentation Map of brain tissue");
 
  // Define command line arguments
  TCLAP::ValueArg<std::string> bundleFileNameArg("b", "bundles", "Fibers bundles filename (vtk file)", true, "", "string", cmd);
  TCLAP::ValueArg<std::string> mapFileNameArg("p", "probability-map", "Probability segmentation Map 4D image(nifti file)", true, "", "string", cmd);  
  TCLAP::ValueArg<std::string> outputFileNameArg("o", "output", "fibers clustering (vtk file)", false, "", "string", cmd);

  // Parse arguments
  cmd.parse(argc, argv);

  // Define command line variables
  // Get back arguments' values        
  std::string bundleFileName = bundleFileNameArg.getValue();
  std::string mapFileName    = mapFileNameArg.getValue();
  std::string outputFileName = outputFileNameArg.getValue();
 
  // Load bundle (vtk file)
  vtkSmartPointer<vtkPolyDataReader> bundleReader =  vtkSmartPointer<vtkPolyDataReader>::New();
  bundleReader->SetFileName(bundleFileName.c_str());
  bundleReader->Update();
  
  vtkSmartPointer<vtkPolyData> bundle = bundleReader->GetOutput();
  vtkSmartPointer<vtkCellArray> lines = bundle->GetLines();

  // Load probability map (4D image)  
  Reader4DType::Pointer reader = Reader4DType::New();
  reader->SetFileName(mapFileName);
  reader->Update();
  Image4DPointer mapImage = reader->GetOutput();
    
  Image4DType::RegionType input4DRegion = mapImage->GetLargestPossibleRegion();
  Image4DType::SizeType input4DSize = input4DRegion.GetSize();  
  Image4DType::IndexType index = input4DRegion.GetIndex();
  
  uint numberOf3Dimages = input4DSize[3];
  
  // Number of fibers
  float number_fiber = bundle->GetNumberOfLines();
  std::cout << "Number of fibers: " << number_fiber << std::endl;
  
  // For each fibers
  vtkIdType numberOfPoints, *pointIds;
  unsigned int i = 0;  

  vtkSmartPointer<vtkDoubleArray> label_data = vtkSmartPointer<vtkDoubleArray>::New();
  label_data->SetNumberOfComponents(1);

  std::cout << "Processing: start..." << std::endl;
  //----------------------------------------------------------------------------
  
  while(lines->GetNextCell(numberOfPoints, pointIds) != 0)
  {
  std::vector<float> label_vector(numberOf3Dimages);
  for(int l=0;l < numberOf3Dimages;l++)
    label_vector[l] = 0;

  for(int p=0; p < (int)numberOfPoints; p++) //For all the points of each fiber
    {
    // Get current point's coordinates
    double worldCoordinates[3];
    bundle->GetPoint(pointIds[p], worldCoordinates);

    // Convert world coordinates to index coordinates (from worldcoordinate to indexIn3DImage)
    Image4DType::PointType worldPoint;
    worldPoint[0] = -worldCoordinates[0]; 
    worldPoint[1] = -worldCoordinates[1]; 
    worldPoint[2] =  worldCoordinates[2];

    Image3DType::Pointer indexIn3DImage = Image3DType::New();

	for(int l=0;l < numberOf3Dimages;l++)
	  {     //Continuous index 
		worldPoint[3] = l;
		Image4DType::IndexType index;
		mapImage->TransformPhysicalPointToIndex(worldPoint, index);

           //Get the value at location (x,y,z,l) from image4D
	   label_vector[l] += mapImage->GetPixel(index);
	  }
    }  
  double maxlabel_point = std::max_element(label_vector.begin(), label_vector.end()) - label_vector.begin();
   
  for(int p=0; p < (int)numberOfPoints; p++)
    label_data->InsertNextTupleValue(&maxlabel_point);   
   
  i++; //next bundle 
  }

//------------------------------------------------------------------------------

// Create the output vtk data
bundle->GetPointData()->SetScalars(label_data);

// Save the output data into VTK file
vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
writer->SetInput(bundle);
// If no filename is given for output, set it up with input name
    if(outputFileName.empty())
    {
        std::stringstream filename;
        filename << "clustering1-"+bundleFileName;
        writer->SetFileName(filename.str().c_str());
        writer->Write();
    }
    else // !outputFileName.empty()
    {
        writer->SetFileName(outputFileName.c_str());
        writer->SetFileTypeToASCII();
        writer->Write();
    }

std::cout << "Save the clustering result --> done." << std::endl;

return EXIT_SUCCESS;
}
