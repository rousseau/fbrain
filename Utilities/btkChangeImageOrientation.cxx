/*==========================================================================

  Date: 22/07/2015
  Author(s): François Rousseau

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
#include <string>
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include "itkOrientImageFilter.h"

/*Btk includes*/
#include "btkImageHelper.h"


int main (int argc, char* argv[])
{
  
  try {  

  TCLAP::CmdLine cmd("It changes image orientation using ITK OrientImageFilter.", ' ', "1.0", true);
  TCLAP::ValueArg<std::string> inputImageArg ("i","input_file","input image file ",true,"", "string", cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> orientationArg("","orientation","desired orientation",true,"","string", cmd);

  // Parse the args.
  cmd.parse( argc, argv );
  
  std::string input_file = inputImageArg.getValue();
  std::string output_file = outputImageArg.getValue();
  std::string orientation = orientationArg.getValue();

  //ITK declaration
  typedef float PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef ImageType::Pointer ImagePointer;
  typedef itk::SpatialOrientation::ValidCoordinateOrientationFlags SpatialOrientation;

  ImagePointer inputImage = btk::ImageHelper<ImageType>::ReadImage(input_file);
  SpatialOrientation itkOrientation;
  if (orientation=="RAS")
    itkOrientation = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
  else if (orientation=="LPS")
    itkOrientation = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS;
  else if (orientation=="LPI")
    itkOrientation = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI;
  else
  {
    std::cout<<"Default mode : RAS"<<std::endl;
    itkOrientation = itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS;
  }

  itk::OrientImageFilter<ImageType,ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType,ImageType>::New();
  orienter->UseImageDirectionOn();
  orienter->SetDesiredCoordinateOrientation(itkOrientation);
  orienter->SetInput(inputImage);
  orienter->Update();
  ImagePointer outputImage = orienter->GetOutput();

  btk::ImageHelper<ImageType>::WriteImage(outputImage, output_file);
    
  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return 1;
}
