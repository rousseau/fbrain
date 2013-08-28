/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 30/05/2013 (last updated: 27/06/2013)
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
typedef itk::ImageFileWriter< Image4DType >  Writer4DType;

typedef itk::Image<PixelType,Dimension-1> Image3DType;
typedef Image3DType::Pointer Image3DPointer;
typedef itk::ImageFileReader< Image3DType >  Reader3DType;
typedef itk::ImageFileWriter< Image3DType >  Writer3DType;


//
// Usage: btkQuantifyClusteredRegions -i fiber_clustering.vtk -r 4D-empty-image.nii.gz -o labeled_voxels.nii.gz
//


// main
int main ( int argc, char *argv[] )
{

  // Define command line parser
  TCLAP::CmdLine cmd("Fiber Clustering Quantification");
 
  // Define command line arguments
  TCLAP::ValueArg<std::string> inputFileNameArg("i", "input", "Fiber tracts clustering data (vtk file)", true, "", "string", cmd);
  TCLAP::ValueArg<std::string> dataFileNameArg("r", "reference-image", "4D empty image (nifti file)", true, "", "string", cmd);  
  TCLAP::ValueArg<std::string> outputFileNameArg("o", "output", "Labeled voxels (nifti file)", false, "", "string", cmd);

  // Parse arguments
  cmd.parse(argc, argv);

  // Define command line variables
  // Get back arguments' values        
  std::string inputFileName  =  inputFileNameArg.getValue();
  std::string dataFileName   =  dataFileNameArg.getValue();
  std::string outputFileName =  outputFileNameArg.getValue();
 
  // Load fibers (vtk file)-----------------------------------
  vtkSmartPointer<vtkPolyDataReader> fiberReader =  vtkSmartPointer<vtkPolyDataReader>::New();
  fiberReader-> SetFileName(inputFileName.c_str());
  fiberReader-> Update();
  
  vtkSmartPointer<vtkPolyData>  fibers  = fiberReader-> GetOutput();
  vtkSmartPointer<vtkCellArray> lines   = fibers-> GetLines();

  // Number of fibers
  float number_fiber = fibers-> GetNumberOfLines();
  std::cout << "\nNumber of fibers: " << number_fiber << std::endl;

  // Reading reference image  (4D image)--------------------------  
  // This 4D empty image is generated using:
  // btkCreate4DEmptyImageFromReferenceImage -i 01019-natbrain4D.nii.gz -o 4D-empty-image.nii.gz
 
  Reader4DType::Pointer reader = Reader4DType::New();
  reader->SetFileName(dataFileName);
  reader->Update();
  Image4DPointer dataImage = reader->GetOutput();
  
  Image4DType::RegionType input4DRegion = dataImage->GetLargestPossibleRegion();
  Image4DType::SizeType input4DSize = input4DRegion.GetSize();  
  Image4DType::IndexType index = input4DRegion.GetIndex();
  
  uint numberOf3Dimages = input4DSize[3];
  
  for(int n=0;n < numberOf3Dimages;n++)
 {  
	for(int nx=0; nx < input4DSize[0];nx++)
	for(int ny=0; ny < input4DSize[1];ny++)
	for(int nz=0; nz < input4DSize[2];nz++)
	{Image4DType::IndexType index;
        index[0]=nx; index[1]=ny;index[2]=nz;index[3]=n;
	dataImage->SetPixel(index,0);
	}
 }
  // For each fibers
  vtkIdType numberOfPoints, *pointIds;
  unsigned int i = 0;  

  //----------------------------------------------------------------------------
  std::cout << "Processing...";
  int j=0, label;
  while(lines->GetNextCell(numberOfPoints, pointIds) != 0)
  {
  for(int p=0; p < (int)numberOfPoints; p++) //For all the points of each fiber
    {
    // Get current point's coordinates
    double worldCoordinates[3];
    fibers-> GetPoint(pointIds[p], worldCoordinates);

    // Convert world coordinates to index coordinates (from worldcoordinate to indexIn3DImage)
    Image4DType::PointType worldPoint;
    worldPoint[0] = -worldCoordinates[0]; 
    worldPoint[1] = -worldCoordinates[1]; 
    worldPoint[2] =  worldCoordinates[2];

    label = fibers->GetPointData()->GetScalars()->GetTuple1(j);j++;
    
    Image3DType::Pointer indexIn3DImage = Image3DType::New();

	for(int l=0;l < numberOf3Dimages;l++)
	  {     //Continuous index 
		worldPoint[3] = l;
		Image4DType::IndexType index;
		dataImage->TransformPhysicalPointToIndex(worldPoint, index);
		dataImage->SetPixel(index,label+1);
	  }
    }  

  i++; //next bundle 
  }
std::cout << " done." << std::endl;
//------------------------------------------------------------------------------
   //Write the result
    Writer4DType::Pointer writer = Writer4DType::New();
    writer->SetFileName( outputFileName );
    writer->SetInput( dataImage );
    writer->Update(); 


return EXIT_SUCCESS;
}
