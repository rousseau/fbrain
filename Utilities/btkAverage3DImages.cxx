/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 09/08/2012
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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

/*Btk includes*/
#include "btkImageHelper.h"


int main (int argc, char* argv[])
{
  
  try {  

  TCLAP::CmdLine cmd("It computes the average and the variance at every voxel of a set of 3D images.", ' ', "1.0", true);
  TCLAP::MultiArg<std::string> anatomicalImageArg("a","anatomical_file","anatomical image files (possible multiple inputs) ",true,"string", cmd);
  TCLAP::ValueArg<std::string> meanImageArg("m","mean_output_file","output mean image file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> varianceImageArg("v","variance_output_file","output variance image file",true,"","string", cmd);

  // Parse the args.
  cmd.parse( argc, argv );
  
  std::vector<std::string> anatomical_file = anatomicalImageArg.getValue();
  std::string mean_file = meanImageArg.getValue();
  std::string variance_file = varianceImageArg.getValue();

  //ITK declaration
  typedef float PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef ImageType::Pointer ImagePointer;

  std::vector<ImagePointer> anatomicalImages;
  anatomicalImages.resize(anatomical_file.size());

  anatomicalImages = btk::ImageHelper<ImageType>::ReadImageArray(anatomical_file);
  
  ImageType::RegionType region;
  ImageType::SizeType size;
  ImageType::SpacingType spacing;

  //we assume that all images have the same size, spacing etc.
  region = anatomicalImages[0]->GetLargestPossibleRegion();
  size = region.GetSize();
  spacing = anatomicalImages[0]->GetSpacing();

  //Allocate output images
  ImagePointer meanImage = ImageType::New();
  meanImage->SetRegions( region );
  meanImage->SetSpacing( spacing );
  meanImage->SetOrigin( anatomicalImages[0]->GetOrigin() );
  meanImage->SetDirection( anatomicalImages[0]->GetDirection() );
  meanImage->Allocate();
  meanImage->FillBuffer(0); 

  ImagePointer varianceImage = ImageType::New();
  varianceImage->SetRegions( region );
  varianceImage->SetSpacing( spacing );
  varianceImage->SetOrigin( anatomicalImages[0]->GetOrigin() );
  varianceImage->SetDirection( anatomicalImages[0]->GetDirection() );
  varianceImage->Allocate();
  varianceImage->FillBuffer(0); 

  int x,y,z;
  uint n = anatomicalImages.size();

  #pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)size[2]; z++)
    for(y=0; y < (int)size[1]; y++)
      for(x=0; x < (int)size[0]; x++){

        ImageType::IndexType p;
        p[0] = x;
        p[1] = y;
        p[2] = z;

        float m = 0;
        float m2= 0;
        float val = 0;

        for(uint i=0; i < n; i++){
            val = anatomicalImages[i]->GetPixel(p);
            m += val;
            m2+= val*val;
        }
        float mean = m / n;
        float variance = (m2 / n) - (mean * mean) ;
            
        meanImage->SetPixel( p, mean );	  	  
        varianceImage->SetPixel( p, variance );
  }
  
  btk::ImageHelper<ImageType>::WriteImage(meanImage, mean_file);
  btk::ImageHelper<ImageType>::WriteImage(varianceImage, variance_file);
    
  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return 1;
}
