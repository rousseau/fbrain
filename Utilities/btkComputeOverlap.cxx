/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 10/09/2010
  Author(s): François Rousseau (rousseau@unistra.fr)

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

/* Standard includes */
#include <string>
#include <iomanip>
#include <tclap/CmdLine.h>


/* Itk includes */
#include "itkImage.h"
#include "itkMinimumMaximumImageCalculator.h"

/*Btk includes*/
#include "btkImageHelper.h"

int main (int argc, char* argv[])
{
  
  try {  
  
  TCLAP::CmdLine cmd("It computes the overlap (for a set of labels) between two 3D images.", ' ', "1.0", true);
  TCLAP::ValueArg<std::string> groundTruthImageArg("g","ground_truth_file","input ground truth file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> inputImageArg("i","input_file","input segmentation file",true,"","string", cmd);
  TCLAP::ValueArg<int> minLabelArg("","min_label","minimum label (default : compute from the ground truth image)", false,-1,"int",cmd);
  TCLAP::ValueArg<int> maxLabelArg("","max_label","maximum label (default : compute from the ground truth image)", false,-1,"int",cmd);

  // Parse the args.
  cmd.parse( argc, argv );  

  std::string ground_truth = groundTruthImageArg.getValue();
  std::string input_file = inputImageArg.getValue();
  int minLabel = minLabelArg.getValue();
  int maxLabel = maxLabelArg.getValue();
    
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  
  //Reading the 2 input files
  ImageType::Pointer ground_truth_image = btk::ImageHelper<ImageType>::ReadImage(ground_truth);
  ImageType::Pointer input_image = btk::ImageHelper<ImageType>::ReadImage(input_file);

  itk::MinimumMaximumImageCalculator< ImageType >::Pointer minmax = itk::MinimumMaximumImageCalculator< ImageType >::New();
  minmax->SetImage(ground_truth_image);
  minmax->Compute();
    
  if(minLabel == -1) minLabel = minmax->GetMinimum();
  if(maxLabel == -1) maxLabel = minmax->GetMaximum();

  std::map<PixelType,uint> intersection;          //intersection between the 2 segmentations (i.e. true positive)
  std::map<PixelType,uint> cardinal_ground_truth;
  std::map<PixelType,uint> cardinal_input_image;

  typedef itk::ImageRegionIterator< ImageType > Iterator;

  Iterator itGroundTruth( ground_truth_image, ground_truth_image->GetLargestPossibleRegion() );
  Iterator itImage( input_image, input_image->GetLargestPossibleRegion() );

  for(itGroundTruth.GoToBegin(), itImage.GoToBegin(); !itGroundTruth.IsAtEnd(); ++itGroundTruth, ++itImage){

    PixelType ground_truth_label = itGroundTruth.Get();
    PixelType input_label = itImage.Get();

    if( ground_truth_label == input_label ) intersection[ground_truth_label] += 1;

    cardinal_input_image[input_label] += 1;
    cardinal_ground_truth[ground_truth_label] += 1;
  
  }

  std::map<PixelType, uint>::iterator cgtIt; //iterator over the map "cardinal_ground_truth"
  
  //std::cout<<"label | kappa index | overlap |\n";
  for(int i = minLabel; i < maxLabel + 1; i++){
    float overlap = 0;   //overlap or jaccard index = intersection / union = TP / (TP+FP+FN) = kappa / (2-kappa)
    float kappa = 0;     //kappa index or dice coefficient = 2* intersection / ( sum of cardinals) = 2*TP / (2*TP+FP+FN) = 2 * overlap / (overlap+1)
    
    PixelType label = i;

    kappa = 2.0 * intersection[label] / ( cardinal_ground_truth[label] + cardinal_input_image[label] );
    overlap = kappa / (2.0 - kappa);
    std::cout<<label<<" "<<kappa<<" "<<overlap<<"\n";
  }

    
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return 1;
}




