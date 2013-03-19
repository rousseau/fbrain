/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 14/01/2013
 Author(s): Larbi Boubchir, François Rousseau
 
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


//ITK STUFF
typedef float PixelType;
const unsigned int Dimension = 4;
typedef itk::Image< PixelType,Dimension > Image4DType;
typedef Image4DType::Pointer Image4DPointer;
typedef itk::ImageFileReader< Image4DType >  Reader4DType;
typedef itk::Image<PixelType,Dimension-1> Image3DType;
typedef Image3DType::Pointer Image3DPointer;

double meanDistanceBetweenFibers(std::vector<Image4DType::PointType> & f1, std::vector<Image4DType::PointType> & f2)
{
  double distance = 0;
  unsigned int length = f1.size();
  for(unsigned int p=0; p < length; p++)
  {
    double tmp0 = f1[p][0] - f2[p][0];
    double tmp1 = f1[p][1] - f2[p][1];
    double tmp2 = f1[p][2] - f2[p][2];
          
    distance += (tmp0*tmp0 + tmp1*tmp1 + tmp2*tmp2) / (1.0 * length );  
  }   
  return distance;
}


// main
int main ( int argc, char *argv[] )
{
  try{
    // Define command line parser
    TCLAP::CmdLine cmd("White-matter fibers labeling using atlas priors");

    // Define command line arguments
    TCLAP::ValueArg<std::string>   bundleFileNameArg("b", "bundles", "Fibers bundles filename (vtk file)", true, "", "string", cmd);
    TCLAP::ValueArg<std::string>   mapFileNameArg   ("p", "probability-map", "Probability Map 4D image (nifti file)", true, "", "string", cmd);
    TCLAP::ValueArg<std::string>   outputFileNameArg("o", "output", "fibers clustering (vtk file)", false, "", "string", cmd);
    TCLAP::ValueArg<unsigned int>  samplingArg      ("s", "sampling", "Number of points of fibers during clustering (10)", false, 10, "int", cmd);
    

    // Parse arguments
    cmd.parse(argc, argv);

    // Get back arguments' values
    std::string bundleFileName = bundleFileNameArg.getValue();
    std::string mapFileName    = mapFileNameArg.getValue();
    std::string outputFileName = outputFileNameArg.getValue();
    unsigned int samplingValue = samplingArg.getValue();

    std::cout<<"Loading bundles :"<<bundleFileName<<std::endl;
    vtkSmartPointer<vtkPolyDataReader> bundleReader =  vtkSmartPointer<vtkPolyDataReader>::New();
    bundleReader->SetFileName(bundleFileName.c_str());
    bundleReader->Update();

    vtkSmartPointer<vtkPolyData> bundle = bundleReader->GetOutput();
    vtkSmartPointer<vtkCellArray> lines = bundle->GetLines();
    
    unsigned int numberOfFibers = lines->GetNumberOfCells();
    std::cout<<"Number of fibers: "<< numberOfFibers <<std::endl;

    std::cout<<"Loading probability map (4D image) : "<<mapFileName<<std::endl;
    Reader4DType::Pointer reader = Reader4DType::New();
    reader->SetFileName(mapFileName);
    reader->Update();
    Image4DPointer mapImage = reader->GetOutput();

    Image4DType::RegionType input4DRegion = mapImage->GetLargestPossibleRegion();
    Image4DType::SizeType input4DSize = input4DRegion.GetSize();
    Image4DType::IndexType index = input4DRegion.GetIndex();

    unsigned int numberOfLabels = input4DSize[3];
    std::cout << "Number of Labels : " << numberOfLabels << std::endl;

    // For each fibers
    vtkIdType numberOfPoints, *pointIds;
    unsigned int fiberIndex = 0;

    //std::vector<float> label_final_point;

    vtkSmartPointer<vtkDoubleArray> label_data = vtkSmartPointer<vtkDoubleArray>::New();
    label_data->SetNumberOfComponents(1);
    vtkSmartPointer<vtkDoubleArray> fcm_label = vtkSmartPointer<vtkDoubleArray>::New();
    fcm_label->SetNumberOfComponents(1);



    std::vector< std::vector<Image4DType::PointType> > undersampledBundle(numberOfFibers);
    std::vector< unsigned int> fiberLengths(numberOfFibers);
    
    std::vector< std::vector< float > > labelWeightOfFibers(numberOfFibers);
    for(unsigned int i=0; i < numberOfFibers; i++)
    {
      labelWeightOfFibers[i].resize(numberOfLabels);
      for(int l=0;l < numberOfLabels;l++)
        labelWeightOfFibers[i][l] = 0;
    }
    
    std::vector<double> maximumScore(numberOfLabels);
    std::vector<int> indexWithMaximumScore(numberOfLabels,-1);
    unsigned int minimumLength = 100000000;
    unsigned int maximumLength = 0;
    std::vector<unsigned int> fiberLength(numberOfFibers,0);

    while(lines->GetNextCell(numberOfPoints, pointIds) != 0)
    {

      //for undersampling
      std::vector<Image4DType::PointType> tmpVector;
      fiberLengths[fiberIndex] = (int)numberOfPoints;

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

        for(int l=0;l < numberOfLabels;l++)
        {     
          //Continuous index
          worldPoint[3] = l;
          Image4DType::IndexType index;
          mapImage->TransformPhysicalPointToIndex(worldPoint, index);

          labelWeightOfFibers[fiberIndex][l] += mapImage->GetPixel(index);

        //Undersampling
        tmpVector.push_back(worldPoint);
        //std::cout<<worldPoint<<"\n";
        
        }        
      }
      
      if(minimumLength > (int)numberOfPoints) minimumLength = (int)numberOfPoints;
      if(maximumLength < (int)numberOfPoints) maximumLength = (int)numberOfPoints;
      fiberLength[fiberIndex] = (int)numberOfPoints;
      
      double maxlabel_point = std::max_element(labelWeightOfFibers[fiberIndex].begin(), labelWeightOfFibers[fiberIndex].end()) - labelWeightOfFibers[fiberIndex].begin();
      for(int p=0; p < (int)numberOfPoints; p++)
        label_data->InsertNextTupleValue(&maxlabel_point);

      double score = labelWeightOfFibers[fiberIndex][maxlabel_point];
      if(score > maximumScore[maxlabel_point])
      {
        maximumScore[maxlabel_point] = score;
        indexWithMaximumScore[maxlabel_point] = fiberIndex;  
      }

      //undersampling
      unsigned int fiberLength = tmpVector.size();
      undersampledBundle[fiberIndex].resize(samplingValue);
      for(unsigned int j=0; j<samplingValue; j++)
        undersampledBundle[fiberIndex][j] = tmpVector[(int)(j*(fiberLength-1)/(samplingValue-1))];

      fiberIndex++; //next bundle

    }


    //std::cout << "Number of tuples: " << label_data->GetNumberOfTuples() << std::endl;
    std::cout << "Number of points: " << bundle->GetPoints()->GetNumberOfPoints() << std::endl;

    std::cout<<"Lenghts of fibers: "<<maximumLength<<" "<<minimumLength<<std::endl;

    for(unsigned int l=0; l < numberOfLabels; l++)
    {
      std::cout<<"For label "<<l<<", max score : "<<maximumScore[l]<<", corresponding fiber : "<<indexWithMaximumScore[l];
      if( indexWithMaximumScore[l] >= 0 )
        std::cout<<" of size "<<fiberLengths[ indexWithMaximumScore[l] ]<<", "<<maximumScore[l]*100.0/fiberLengths[ indexWithMaximumScore[l] ]<<"%"<<std::endl;
      else
        std::cout<<std::endl;  
      
    }
    
    int myF = 13462;
    std::cout<<"Fiber 13462 \n";
    for(unsigned int l=0; l < numberOfLabels; l++)
      std::cout<<labelWeightOfFibers[myF][l]<<" ";
    std::cout<<"\n";  

    //FCM ------------------------------------------------------------------------------------------------
    int iterMax = 5;
    
    //Initialisation centroids
    std::vector< std::vector<Image4DType::PointType> > centroid(numberOfLabels);
    for(unsigned int l=0; l < numberOfLabels; l++)
    {
      centroid[l].resize(samplingValue);
      double sumOfWeights = 0;
      std::vector<Image4DType::PointType> vectorPoint(samplingValue);
      for(unsigned int p=0; p < samplingValue; p++)
      {
        Image4DType::PointType point;
        point[0] = 0;      point[1] = 0;      point[2] = 0;      point[3] = 0;
        vectorPoint[p] = point;        
      }
      //Weighted init -------------------------------------
      /*
      for(unsigned int f=0; f < numberOfFibers; f++)
      {  
        double weight = labelWeightOfFibers[f][l];
        sumOfWeights += weight;
        
        for(unsigned int p=0; p < samplingValue; p++)
        {
          vectorPoint[p][0] += weight *  undersampledBundle[f][p][0]; 
          vectorPoint[p][1] += weight *  undersampledBundle[f][p][1]; 
          vectorPoint[p][2] += weight *  undersampledBundle[f][p][2]; 
        }
      }
      
      //normalization of the centroid
      for(unsigned int p=0; p < samplingValue; p++)
      {
        centroid[l][p][0] = vectorPoint[p][0] / sumOfWeights;
        centroid[l][p][1] = vectorPoint[p][1] / sumOfWeights;
        centroid[l][p][2] = vectorPoint[p][2] / sumOfWeights;
      }
      */
      //Hard init
      for(unsigned int p=0; p < samplingValue; p++)
      {
        if( indexWithMaximumScore[l] >= 0 )
        {
          centroid[l][p][0] = undersampledBundle[indexWithMaximumScore[l]][p][0];
          centroid[l][p][1] = undersampledBundle[indexWithMaximumScore[l]][p][1];
          centroid[l][p][2] = undersampledBundle[indexWithMaximumScore[l]][p][2];
        }
        else
        {
          centroid[l][p][0] = 0;
          centroid[l][p][1] = 0;
          centroid[l][p][2] = 0;
        }
      }
      
    }

    //Cost Function
    double cost = 0;
    for(unsigned int f=0; f < numberOfFibers; f++)
      for(unsigned int l=0; l < numberOfLabels; l++)    
        cost += labelWeightOfFibers[f][l] * labelWeightOfFibers[f][l] * meanDistanceBetweenFibers( undersampledBundle[f], centroid[l] );  
    
    std::cout<<"Cost function : "<<cost<<std::endl;


    std::cout<<"Centroid du label 0 :\n";
    for(unsigned int p=0; p < samplingValue; p++)
      std::cout<<centroid[3][p];
    std::cout<<std::endl;

    for(unsigned int iter=0; iter < iterMax; iter++)
    {
      //Compute new weights
      for(unsigned int f=0; f < numberOfFibers; f++)
      {
        std::vector<double> distance(numberOfLabels,0);
        double sumOfDistances = 0;
        for(unsigned int l=0; l < numberOfLabels; l++)
        {
          distance[l] = meanDistanceBetweenFibers( undersampledBundle[f], centroid[l] ); 
                    
          if(distance[l] < 0.0001)
            distance[l] = 0.0001; 
                       
          sumOfDistances += (1.0 / distance[l]);
        }  
               
        for(unsigned int l=0; l < numberOfLabels; l++)
        {
          labelWeightOfFibers[f][l] = 1.0 / ( distance[l] * sumOfDistances );
          
          if(f==13462)
            std::cout<<"f 13462: "<<distance[0]<<" "<<labelWeightOfFibers[f][0]<<" "<<distance[1]<<" "<<labelWeightOfFibers[f][1]<<" "<<distance[2]<<" "<<labelWeightOfFibers[f][2]<<" "<<std::endl;
          
        }        
      }
      
      //Update centroids  
      for(unsigned int l=0; l < numberOfLabels; l++)
      {
        double sumOfWeights = 0;
        std::vector<Image4DType::PointType> vectorPoint(samplingValue);
        for(unsigned int p=0; p < samplingValue; p++)
        {
          Image4DType::PointType point;
          point[0] = 0;      point[1] = 0;      point[2] = 0;      point[3] = 0;
          vectorPoint[p] = point;        
        }
        for(unsigned int f=0; f < numberOfFibers; f++)
        {  
          double weight = labelWeightOfFibers[f][l];
          sumOfWeights += weight;
        
          for(unsigned int p=0; p < samplingValue; p++)
          {
            vectorPoint[p][0] += weight *  undersampledBundle[f][p][0]; 
            vectorPoint[p][1] += weight *  undersampledBundle[f][p][1]; 
            vectorPoint[p][2] += weight *  undersampledBundle[f][p][2]; 
          }
        }
      
        //normalization of the centroid
        for(unsigned int p=0; p < samplingValue; p++)
        {
          centroid[l][p][0] = vectorPoint[p][0] / sumOfWeights;
          centroid[l][p][1] = vectorPoint[p][1] / sumOfWeights;
          centroid[l][p][2] = vectorPoint[p][2] / sumOfWeights;
        }
      }
      
    //Cost Function
    double cost = 0;
    for(unsigned int f=0; f < numberOfFibers; f++)
      for(unsigned int l=0; l < numberOfLabels; l++)    
        cost += labelWeightOfFibers[f][l] * labelWeightOfFibers[f][l] * meanDistanceBetweenFibers( undersampledBundle[f], centroid[l] );  
    
    std::cout<<"Cost function : "<<cost<<std::endl;

      
      std::cout<<"Fiber 13462 \n";
      for(unsigned int l=0; l < numberOfLabels; l++)
        std::cout<<labelWeightOfFibers[myF][l]<<" ";
      std::cout<<"\n";  

      
      std::cout<<"Centroid du label 0 :\n";
      for(unsigned int p=0; p < samplingValue; p++)
        std::cout<<centroid[3][p];
      std::cout<<std::endl;
                 
    }
    
    
    for(unsigned int f=0; f < numberOfFibers; f++)
    {
      double maxlabel_point = std::max_element(labelWeightOfFibers[f].begin(), labelWeightOfFibers[f].end()) - labelWeightOfFibers[f].begin();        
      for(int p=0; p < fiberLength[f]; p++)
        fcm_label->InsertNextTupleValue(&maxlabel_point);
        
      std::cout<<"Fibre "<<f<<" : "<<maxlabel_point<<std::endl;      
    }
    std::cout<<"Nb fibres labellisees : "<<fiberIndex<<std::endl;
    // FIN DE FCM --------------------------------------------------------------------------------------------------
    
    //QB:
    //initialiser avec les fibres qui appartiennent au moins à un certain seuil
    //Ensuite, aapliquer QB
    //Relabelliser avec l'atlas eventuellement les clusters non labelisés
    
    
    //LABELING
    // create output vtk data
    //bundle->GetPointData()->SetScalars(label_data);
    bundle->GetPointData()->SetScalars(fcm_label);
    
    std::cout << "Number of tuples: " << label_data->GetNumberOfTuples() << std::endl;
    std::cout << "Number of tuples: " << fcm_label->GetNumberOfTuples() << std::endl;


    // Save output data into VTK file
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInput(bundle);
    // If no filename is given for output, set it up with input name
    if(outputFileName.empty())
    {
      std::stringstream filename;
      filename << "clustering-"+bundleFileName;
      writer->SetFileName(filename.str().c_str());
      writer->Write();
    }
    else // !outputFileName.empty()
    {
      writer->SetFileName(outputFileName.c_str());
      writer->SetFileTypeToASCII();
      writer->Write();
    }
  std::cout << "done." << std::endl;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
 
  return EXIT_SUCCESS;
}