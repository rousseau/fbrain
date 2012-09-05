/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 05/09/2012
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
#include <vector>
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include <itkHausdorffDistanceImageFilter.h>
/*Btk includes*/
#include "btkImageHelper.h"


int main (int argc, char* argv[])
{
  
  //This is an ITK filter
  //Computes the Hausdorff distance between the set of non-zero pixels of two images.
  
  try {  

  TCLAP::CmdLine cmd("Computes the Hausdorff distance between the set of non-zero pixels of two images.", ' ', "1.0", true);
  TCLAP::MultiArg<std::string> inputImageArg("i","input_file","input anatomical image files (2 inputs required) ",true,"string", cmd);
  TCLAP::ValueArg< bool > spacingArg("s","spacing","use spacing (default is true)",false,true,"bool",cmd);


  // Parse the args.
  cmd.parse( argc, argv );
  
  std::vector<std::string> input_file    = inputImageArg.getValue();
  bool use_spacing                       = spacingArg.getValue();
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef ImageType::Pointer                    ImagePointer;
  typedef itk::HausdorffDistanceImageFilter<ImageType, ImageType> FilterType;

  std::cout<<"Read input images \n";
  std::vector<ImagePointer> inputImages;
  inputImages.resize(input_file.size());
  inputImages = btk::ImageHelper<ImageType>::ReadImage(input_file);

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput1(inputImages[0]);
  filter->SetInput2(inputImages[1]);
  filter->SetUseImageSpacing(use_spacing);
  filter->Update();

  std::cout<<"Hausdorff distance : "<<filter->GetHausdorffDistance()<<std::endl;
  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return 1;
}
