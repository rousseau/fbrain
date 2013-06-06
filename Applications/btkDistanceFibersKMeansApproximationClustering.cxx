/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 14/02/2013 (last updated: 28/05/2013)
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
#include <cstdlib>
#include <string>
#include <iomanip>
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
#include "vtkTable.h"
#include "vtkKMeansStatistics.h"

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
const unsigned int Dimension = 3;
typedef itk::Image<PixelType,Dimension> Image3DType;
typedef Image3DType::Pointer Image3DPointer;

//
// Function declaration
//
float Distance_measure(const std::vector<float>& f1,const std::vector<float>& f2,const char& type);
float Orientation_measure(const std::vector<float>& f1,const std::vector<float>& f2);

//
// Distance measure (Asymmetric/symmetric Chamfer or Hausdorff distance)
float Distance_measure(const std::vector<float>& f1,const std::vector<float>& f2,const char& type)
{
    float distance;
	
    if(f1==f2) return 0; //same fiber --> distance = 0
    else	
    {
	std::vector<float> di(f1.size()/3);
	for(int i=0; i< f1.size()/3; i++)
	{	
		std::vector<float> dj(f2.size()/3);	
		for(int j=0; j< f2.size()/3; j++)
		{
		dj[j] = sqrt(pow(f1[3*i+0]-f2[3*j+0],2)+pow(f1[3*i+1]-f2[3*j+1],2)+pow(f1[3*i+2]-f2[3*j+2],2));
		}
		di[i] = dj[min_element(dj.begin(), dj.end())-dj.begin()];
                dj.clear();
	}	

	if(type == 'C')
	{ 
		distance = accumulate(di.begin(), di.end(),0.0)/di.size();
	        di.clear(); 
	}
	
	if(type == 'H')
	{
		distance = std::max_element(di.begin(), di.end()) - di.begin();
        	di.clear();
	}

	return distance;
     }
}

//
// Local orientation measure
float Orientation_measure(const std::vector<float>& f1,const std::vector<float>& f2)
{    
  if(f1==f2) return 0; //same fiber --> Local orientation = 0
  else	
  {    
        int length, pas;     
        if(f1.size()<f2.size()) 
	{
	  length = f1.size()/3;
	  pas = (int)(f2.size()/f1.size());
	}
        else 
	{
	  length = f2.size()/3;
	  pas= (int)(f1.size()/f2.size());
	}
	
	float V1x, V1y, V1z, V2x, V2y, V2z;
	float orientation = 0;
	
	for(int i=0; i< length; i++)
	{
	  if(f1.size()<f2.size())
	  { 
	    V1x = f1[3*(i+1)+0]-f1[3*i+0]; 
	    V1y = f1[3*(i+1)+1]-f1[3*i+1]; 
	    V1z = f1[3*(i+1)+2]-f1[3*i+2];  
	    V1x = V1x/sqrt(pow(V1x,2) + pow(V1y,2) + pow(V1z,2));
	    V1y = V1y/sqrt(pow(V1x,2) + pow(V1y,2) + pow(V1z,2));
	    V1z = V1z/sqrt(pow(V1x,2) + pow(V1y,2) + pow(V1z,2));

	    V2x = f2[3*(pas*(i+1))+0]-f2[3*(i*pas)+0]; 
	    V2y = f2[3*(pas*(i+1))+1]-f2[3*(i*pas)+1]; 
	    V2z = f2[3*(pas*(i+1))+2]-f2[3*(i*pas)+2];
	    V2x = V2x/sqrt(pow(V2x,2) + pow(V2y,2) + pow(V2z,2));
	    V2y = V2y/sqrt(pow(V2x,2) + pow(V2y,2) + pow(V2z,2));
	    V2z = V2z/sqrt(pow(V2x,2) + pow(V2y,2) + pow(V2z,2));
	    
	    orientation += 1-abs((V1x*V2x) + (V1y*V2y) + (V1z*V2z));
	  }
	  else
	  {
	    V2x = f2[3*(i+1)+0]-f2[3*i+0]; 
	    V2y = f2[3*(i+1)+1]-f2[3*i+1]; 
	    V2z = f2[3*(i+1)+2]-f2[3*i+2];
	    V2x = V2x/sqrt(pow(V2x,2) + pow(V2y,2) + pow(V2z,2));
	    V2y = V2y/sqrt(pow(V2x,2) + pow(V2y,2) + pow(V2z,2));
	    V2z = V2z/sqrt(pow(V2x,2) + pow(V2y,2) + pow(V2z,2));
	    
	    V1x = f1[3*(pas*(i+1))+0]-f1[3*(i*pas)+0]; 
	    V1y = f1[3*(pas*(i+1))+1]-f1[3*(i*pas)+1]; 
	    V1z = f1[3*(pas*(i+1))+2]-f1[3*(i*pas)+2];
	    V1x = V1x/sqrt(pow(V1x,2) + pow(V1y,2) + pow(V1z,2));
	    V1y = V1y/sqrt(pow(V1x,2) + pow(V1y,2) + pow(V1z,2));
	    V1z = V1z/sqrt(pow(V1x,2) + pow(V1y,2) + pow(V1z,2));
	    
	    orientation += 1-abs((V1x*V2x) + (V1y*V2y) + (V1z*V2z));
	  }  
	}
    return orientation; 
  }
}

//
// Usage   : btkDistanceFibersKMeansApproximationClustering -b inputfile.vtk -d distance_measure -c number_of_clusters -a aplha -e eps -o outputfile.vtk
// Example : btkDistanceFibersKMeansApproximationClustering -b data.vtk -d 7 -c 5 -o clustering2-data.vtk

//
// main
int main ( int argc, char *argv[] )
{
  // Define command line parser
  TCLAP::CmdLine cmd("White-matter fiber tracts clustering based on a combination of the distance and orientation measures and the k-means approximation algorithm");
 
  // Define command line arguments
  TCLAP::ValueArg<std::string> bundleFileNameArg("b", "bundles", "Fiber tracts filename (vtk file)", true, "", "string", cmd);
  TCLAP::ValueArg<int> ChoiceDistanceArg("d", "distance", "Distance metric", false, 1, "int",cmd);
  TCLAP::ValueArg<int> NumberOfClustersArg("c", "clusters", "Number of clusters", false, 1, "int",cmd);
  TCLAP::ValueArg<float> AlphaArg("a", "alpha", "Parameter in [0,1]: used to control the balance between the distance and orientation similarities", false, 0.5, "float",cmd);
  TCLAP::ValueArg<float> EpsilonArg("e", "epsilon", "Parameter in ]0,1[: used to fix the sampling parameter", false, 0.5, "float",cmd);
  TCLAP::ValueArg<std::string> outputFileNameArg("o", "output", "Fibers clustering (vtk file)", false, "", "string", cmd);

  // Parse arguments
  cmd.parse(argc, argv);

  // Define command line variables
  // Get back arguments' values        
  std::string bundleFileName = bundleFileNameArg.getValue();
  int ChoiceDistance = ChoiceDistanceArg.getValue();
  int NumberOfClusters = NumberOfClustersArg.getValue();
  float alpha = AlphaArg.getValue();
  float eps = EpsilonArg.getValue();
  std::string outputFileName = outputFileNameArg.getValue();
 
  // Load bundle (vtk file)
  vtkSmartPointer<vtkPolyDataReader> bundleReader =  vtkSmartPointer<vtkPolyDataReader>::New();
  bundleReader->SetFileName(bundleFileName.c_str());
  bundleReader->Update();
  
  vtkSmartPointer<vtkPolyData> bundle = bundleReader->GetOutput();
  vtkSmartPointer<vtkCellArray> lines = bundle->GetLines();
 
  // For each fibers
  vtkIdType numberOfPoints, *pointIds;
  unsigned int i = 0;  

  float number_fiber = bundle->GetNumberOfLines();
  std::vector< std::vector<float> > tab_data(number_fiber,std::vector<float>(0));

  while(lines->GetNextCell(numberOfPoints, pointIds) != 0)
  {     
        tab_data[i].resize( 3*numberOfPoints , 0 );
  	for(int j=0; j < (int)numberOfPoints; j++) //For all the points of each fiber
  	{ // Get current point's coordinates
    	double worldCoordinates[3];
    	bundle->GetPoint(pointIds[j], worldCoordinates);
    	// Convert world coordinates to index coordinates (from worldcoordinate to indexIn3DImage)
    	Image3DType::PointType idx;
    	idx[0] = -worldCoordinates[0]; 
    	idx[1] = -worldCoordinates[1]; 
    	idx[2] =  worldCoordinates[2];
    	tab_data[i][3*j+0] = idx[0];
    	tab_data[i][3*j+1] = idx[1];
    	tab_data[i][3*j+2] = idx[2];
   	}
  i++; //next bundle 
  }

std::cout << "\nNumber of fibers: " << bundle->GetNumberOfLines();

//------------------------------------------------------------------------------
// Distance data matrix
//------------------------------------------------------------------------------
//
// OPTION: calculation of the threshold (--> this is applied with the case 10)
//------------------------------------------------------------------------------
vnl_matrix<float> dist( number_fiber, number_fiber );
float maxd = 0, threshold;
if (ChoiceDistance == 10)
{
  for( int c = 0; c < number_fiber; ++c )
  {
    for( int r = c; r < number_fiber; ++r )
    {
      if (c==r) dist(c,r) = 0;
      else dist(c,r) = Distance_measure(tab_data[c],tab_data[r],'C');
	
      if (maxd < dist(c,r)) maxd = dist(c,r);
      dist(r,c) = dist(c,r);
    }
      std::cout<<"\rComputation of the threshold: "
	   <<100*(float)(c+1)/(float)number_fiber
           <<"% "
           <<std::flush;
  }

std::cout << " --> done" << std::endl;

std::cout << "Maximum distance between fibers: " << maxd << std::endl;
threshold = maxd/NumberOfClusters;
std::cout << "Threshold: " << threshold;
}
//
//------------------------------------------------------------------------------
//
vnl_matrix<float> Distance_matrix( number_fiber, number_fiber );
std::cout << std::endl;

for( int c = 0; c < number_fiber; ++c )
{   
    for( int r = 0; r < number_fiber; ++r )
    { 
	switch(ChoiceDistance) 
	{
		case 1: Distance_matrix(r,c) = Distance_measure(tab_data[c],tab_data[r],'C');break;

		case 2: Distance_matrix(r,c) = (Distance_measure(tab_data[c],tab_data[r],'C')+Distance_measure(tab_data[r],tab_data[c],'C'))/2;break;
	
		case 3: Distance_matrix(r,c) = Distance_measure(tab_data[c],tab_data[r],'H');break;
		
		case 4: Distance_matrix(r,c) = (Distance_measure(tab_data[c],tab_data[r],'H')+Distance_measure(tab_data[r],tab_data[c],'H'))/2; break;

		case 5: Distance_matrix(r,c) = (alpha*Distance_measure(tab_data[c],tab_data[r],'H'))+((1-alpha)*(tab_data[c],tab_data[r],'C')); break;

		case 6: Distance_matrix(r,c) = (alpha*(Distance_measure(tab_data[c],tab_data[r],'H')+Distance_measure(tab_data[r],tab_data[c],'H'))/2)+((1-alpha)*(Distance_measure(tab_data[c],tab_data[r],'C')+Distance_measure(tab_data[r],tab_data[c],'C'))/2);break;

		case 7: Distance_matrix(r,c) = alpha*Orientation_measure(tab_data[c],tab_data[r])+((1-alpha)*(Distance_measure(tab_data[c],tab_data[r],'C')+Distance_measure(tab_data[r],tab_data[c],'C'))/2);break;

		case 8: Distance_matrix(r,c) = alpha*Orientation_measure(tab_data[c],tab_data[r])+((1-alpha)*(Distance_measure(tab_data[c],tab_data[r],'H')+Distance_measure(tab_data[r],tab_data[c],'H'))/2);break;
		
		case 9: Distance_matrix(r,c) = Orientation_measure(tab_data[c],tab_data[r]); break;
		
		case 10: 
		  if (dist(c,r) < threshold) Distance_matrix(r,c) = alpha*Orientation_measure(tab_data[c],tab_data[r])+((1-alpha)*dist(c,r));
		  else Distance_matrix(r,c) = dist(c,r);break;
		
		default: std::cout<< "Error: ''the value of d must between 1 and 10'' \n. Choose 				\n" 
				  << "(1)  Asymmetric Chamfer distance								\n"
				  << "(2)  Symmetric Chamfer distance								\n"
				  << "(3)  Asymmetric Hausdorff distance							\n"				  
				  << "(4)  Symmetric Hausdorff distance								\n"
				  << "(5)  Combined Chamfer-Hausdorff distances							\n"				  
				  << "(6)  Combined Symmetric Chamfer-Symmetric Hausdorff distances				\n"
				  << "(7)  Combined Local Orientation measure-Symmetric Chamfer distance			\n"
				  << "(8)  Combined Local Orientation measure-Symmetric Hausdorff distance			\n"
				  << "(9)  Loacl Orientation measure								\n"
				  << "(10) Combined Local Orientation measure(if distance<threshold)-Symmetric Chamfer distance	\n"
				  <<std::endl; exit (EXIT_FAILURE);
	}
    }
    std::cout<<"\rComputation of the distance matrix: "
	     <<100*(float)(c+1)/(float)number_fiber
             <<"% "
             <<std::flush;
}
std::cout << " --> done";

//------------------------------------------------------------------------------
// A randomized feature selection for the k-means clustering 
//------------------------------------------------------------------------------
vnl_svd<float> svd (Distance_matrix);
vnl_matrix<float> R_Distance_matrix(number_fiber,NumberOfClusters);

R_Distance_matrix=svd.V().extract(number_fiber,NumberOfClusters,0,0);
// leverage scores
vnl_vector<float> Score(number_fiber);
for(uint r=0; r<number_fiber; r++)
{
  float s=0;
  for(uint c=0; c<NumberOfClusters; c++)
    s+=pow(R_Distance_matrix(r,c),2);
  Score(r)=sqrt(s)/NumberOfClusters;
}

//sp is sampling parameter and eps is a parameter in ]0,1[
int sp=(int)( NumberOfClusters*log(NumberOfClusters/eps)/pow(eps,2) ); 

if(sp > number_fiber) sp = number_fiber;

std::vector<int> rand_perm;
for(uint i=1; i<number_fiber; ++i) rand_perm.push_back(i);
std::random_shuffle( rand_perm.begin(), rand_perm.end() );// using built-in random generator:

vnl_vector<float> idx(number_fiber);
uint t=0;
for (std::vector<int>::iterator it=rand_perm.begin(); it!=rand_perm.end(); ++it) 
{idx(t)=*it;t++;}

// get the points into the format needed for K-means
vtkSmartPointer<vtkTable> inputData = vtkSmartPointer<vtkTable>::New(); 
std::cout << std::endl;
for( int t = 0; t < sp; ++t )
{   
    std::stringstream colName;
    colName << "distance" << t;
    vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
    doubleArray->SetNumberOfComponents(1);
    doubleArray->SetName( colName.str().c_str() );
    doubleArray->SetNumberOfTuples(number_fiber);
    for( int r = 0; r < number_fiber; ++r )
    { 
      doubleArray->SetValue( r, Distance_matrix(r,idx(t))*pow(sp*Score(idx(t)),-0.5));
    }
    inputData->AddColumn( doubleArray );
        
    std::cout<<"\rComputation of the reduced distance matrix: "
	     <<100*(float)(t+1)/(float)sp
             <<"% "
             <<std::flush;
}

std::cout << " --> done"<< std::endl;
std::cout << "Size of the reduced distance matrix: [" << number_fiber << "," << sp << "]" << std::endl;

//------------------------------------------------------------------------------
// Apply the k-mean clustering
//------------------------------------------------------------------------------

vtkSmartPointer<vtkKMeansStatistics> kMeansStatistics = vtkSmartPointer<vtkKMeansStatistics>::New();

kMeansStatistics->SetInput( vtkStatisticsAlgorithm::INPUT_DATA, inputData );

// K-means initialization
for(int i=0 ; i < sp; ++i) 
kMeansStatistics->SetColumnStatus( inputData->GetColumnName( i ) , 1 );
// You can used these two columns choosen randomly to initialize the k-means
//kMeansStatistics->SetColumnStatus( inputData->GetColumnName( 1 ) , 1 );
//kMeansStatistics->SetColumnStatus( inputData->GetColumnName( sp-10 ) , 1 );

kMeansStatistics->RequestSelectedColumns();
kMeansStatistics->SetDefaultNumberOfClusters( NumberOfClusters );

// Test Learn and Derive options
kMeansStatistics->SetLearnOption( true );
kMeansStatistics->SetDeriveOption( true );
kMeansStatistics->SetTestOption( true );
kMeansStatistics->SetAssessOption( true );
kMeansStatistics->Update() ;

vtkSmartPointer<vtkDoubleArray> clusterArray = vtkSmartPointer<vtkDoubleArray>::New();
clusterArray->SetNumberOfComponents(1);
clusterArray->SetName( "ClusterId" );
 
for(unsigned int r = 0; r < kMeansStatistics->GetOutput()->GetNumberOfRows(); r++)
{
vtkVariant v = kMeansStatistics->GetOutput()->GetValue(r,kMeansStatistics->GetOutput()->GetNumberOfColumns() - 1);
//std::cout << "Fiber " << r << " is in cluster " << v.ToInt() << std::endl; 
clusterArray->InsertNextValue(v.ToDouble());
    
std::cout<<"\rFiber clustering process: "
         <<100*(float)(r+1)/(float)kMeansStatistics->GetOutput()->GetNumberOfRows()
         <<"% "
         <<std::flush;
}

std::cout << " --> done" << std::endl;

//------------------------------------------------------------------------------
// Create the output vtk file
//------------------------------------------------------------------------------

vtkSmartPointer<vtkPolyDataReader> bundleReader2 =  vtkSmartPointer<vtkPolyDataReader>::New();
bundleReader2->SetFileName(bundleFileName.c_str());
bundleReader2->Update();

vtkSmartPointer<vtkPolyData> bundle2 = bundleReader2->GetOutput();
vtkSmartPointer<vtkCellArray> lines2 = bundle2->GetLines();  
vtkIdType numberOfPoints2, *pointIds2;

vtkSmartPointer<vtkDoubleArray> label_data = vtkSmartPointer<vtkDoubleArray>::New();
label_data->SetNumberOfComponents(1);

int i = 0; 
while(lines2->GetNextCell(numberOfPoints2, pointIds2) != 0)
{
double l = clusterArray->GetValue(i);
for(int p=0; p < (int)numberOfPoints2; p++)
label_data->InsertNextTupleValue(&l);   
i++; //next bundle 
}

bundle2->GetPointData()->SetScalars(label_data);

// Save the output data into VTK file
vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
writer->SetInput(bundle2);
// If no filename is given for output, set it up with input name
    if(outputFileName.empty())
    {
        std::stringstream filename;
        filename << "clustering2-"+bundleFileName;
        writer->SetFileName(filename.str().c_str());
        writer->Write();
    }
    else
    {
        writer->SetFileName(outputFileName.c_str());
        writer->SetFileTypeToASCII();
        writer->Write();
    }
std::cout << "Save the clustering result --> done.\n" << std::endl;

return EXIT_SUCCESS;
}
