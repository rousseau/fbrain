/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 03/06/2013 (last updated: 27/08/2013)
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

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// Local includes
#include "btkImageHelper.h"

// Type definitions
typedef float PixelType;
const unsigned int Dimension = 4;
typedef itk::Image< PixelType,Dimension >          ImageType;


// Usage:   btkCreate4DEmptyImageFromReferenceImage -i ref_image.nii.gz -o 4D-empty-image.nii.gz
// Example: btkCreate4DEmptyImageFromReferenceImage -i 01019-natbrain4D.nii.gz -o 4D-empty-image.nii.gz
//
// Test example: btkCreate4DEmptyImageFromReferenceImage -i 01019-natbrain4D.nii.gz -o 4D-empty-image.nii.gz
// The reference 4D image "01019-natbrain4D.nii.gz" is the probabilistic brain tissue segmentation MAP 4D image generated from the NAMIC data (case01019)
//

// main
//
int main ( int argc, char *argv[] )
{
  // Define command line parser
  TCLAP::CmdLine cmd("Create a 4D empty image from the probabilistic brain tissue segmentation MAP 4D image");
 
  // Define command line arguments
  TCLAP::ValueArg<std::string> inputFileNameArg("i", "input", "4D reference image (nifti file)", true, "", "string", cmd);  
  TCLAP::ValueArg<std::string> outputFileNameArg("o", "output", "4D empty image (nifti file)", true, "", "string", cmd);

  // Parse arguments
  cmd.parse(argc, argv);

  // Define command line variables
  // Get back arguments' values        
  std::string inputFileName      =  inputFileNameArg.getValue();
  std::string outputFileName     =  outputFileNameArg.getValue();
 
  
  // Reading images (4D image)--------------------------------------------------
  ImageType::Pointer inputImage = btk::ImageHelper<ImageType>::ReadImage(inputFileName);

  ImageType::Pointer outputImage= btk::ImageHelper<ImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());
  outputImage->FillBuffer(0.0);
 
// Write images-----------------------------------------------------------------
  btk::ImageHelper<ImageType>::WriteImage(outputImage, outputFileName);

return EXIT_SUCCESS;
}
