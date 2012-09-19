/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 20/10/2009
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
#include "itkImageDuplicator.h"
#include "itkApproximateSignedDistanceMapImageFilter.h"
#include "itkFastChamferDistanceImageFilter.h"

/*Btk includes*/
#include "btkImageHelper.h"

void chamfer_distance(itk::Image<short,3>::Pointer & inputImage, itk::Image<short,3>::Pointer & outputImage);

/* It would be nicer to use the ITK filter itkApproximateSignedDistanceMapImageFilter, but it crashes on real images sometimes

  typedef  itk::ApproximateSignedDistanceMapImageFilter< ImageType, FloatImageType  > ApproximateSignedDistanceMapImageFilterType;
  ApproximateSignedDistanceMapImageFilterType::Pointer filter = ApproximateSignedDistanceMapImageFilterType::New();
  filter->SetInput(image);
  filter->SetInsideValue(1);
  filter->SetOutsideValue(0);
  filter->Update();

*/

int main (int argc, char* argv[])
{
  
  try { 
    
  TCLAP::CmdLine cmd("It computes the chamfer distance to an input image.", ' ', "1.0", true);
  TCLAP::ValueArg<std::string> inputImageArg("i","input_file","input 3D image (cast to short)",true,"","string", cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output 3D image (corresponding to the chamfer distance)",true,"","string", cmd);

  // Parse the args.
  cmd.parse( argc, argv );  

  std::string inputFilename = inputImageArg.getValue();
  std::string outputFilename = outputImageArg.getValue();

    
  typedef short      	      PixelType;
  const   unsigned int      Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef itk::Image< float, Dimension >    FloatImageType;
  typedef itk::ImageDuplicator< ImageType > DuplicatorType;

  ImageType::Pointer inputImage = btk::ImageHelper<ImageType>::ReadImage(inputFilename);

  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage( inputImage );
  duplicator->Update();
  ImageType::Pointer outputImage = duplicator->GetOutput();

  chamfer_distance(inputImage, outputImage);

  btk::ImageHelper<ImageType>::WriteImage(outputImage, outputFilename);

  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return 1;
}


void chamfer_distance(itk::Image<short,3>::Pointer & inputImage, itk::Image<short,3>::Pointer & outputImage)
{
  std::cout<<"chamfer_distance\n";

  int radius = 1;
  itk::Image<short,3>::RegionType region = inputImage->GetLargestPossibleRegion();
  itk::Image<short,3>::SizeType size = region.GetSize();

  std::cout<<"Image size is : "<<size[0]<<" "<<size[1]<<" "<<size[2]<<"\n";
  itk::Image<short,3>::IndexType pixelIndex;

  for(uint k=radius; k<size[2]-radius; k++)
  for(uint j=radius; j<size[1]-radius; j++)
  for(uint i=radius; i<size[0]-radius; i++)
  {
    std::vector<short> vec;

    for(int z=-radius;z<radius+1;z++)
    for(int y=-radius;y<radius+1;y++)
    for(int x=-radius;x<radius+1;x++)
    {
      pixelIndex[0] = i+x;
      pixelIndex[1] = j+y;
      pixelIndex[2] = k+z;
      vec.push_back( inputImage->GetPixel( pixelIndex ) );
    }
    std::sort(vec.begin(), vec.end());
    pixelIndex[0] = i;
    pixelIndex[1] = j;
    pixelIndex[2] = k;
    outputImage->SetPixel( pixelIndex, vec[vec.size()/2.0]);
  }

  //--------------------------------------------------------------------
  int m;
  int dist[14],min;
  int d1=3,d2=4,d3=5;


  for(uint k=0; k<size[2]; k++)
  for(uint j=0; j<size[1]; j++)
  for(uint i=0; i<size[0]; i++)
  {
    pixelIndex[0] = i;
    pixelIndex[1] = j;
    pixelIndex[2] = k;
    if(inputImage->GetPixel(pixelIndex) != 0)
      outputImage->SetPixel(pixelIndex, 0);
    else
      outputImage->SetPixel(pixelIndex, 10000);
  }

  std::cout<<"------ Forward pass ------\n";

  for(uint k=1; k<size[2]-1; k++)
  for(uint j=1; j<size[1]-1; j++)
  for(uint i=1; i<size[0]-1; i++)
  {
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k;
    dist[0]=outputImage->GetPixel(pixelIndex);
    pixelIndex[0] = i;    pixelIndex[1] = j-1;  pixelIndex[2] = k;
    dist[1]=outputImage->GetPixel(pixelIndex)+d1;
    pixelIndex[0] = i-1;  pixelIndex[1] = j-1;  pixelIndex[2] = k;
    dist[2]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i-1;  pixelIndex[1] = j;    pixelIndex[2] = k;
    dist[3]=outputImage->GetPixel(pixelIndex)+d1;
    pixelIndex[0] = i-1;  pixelIndex[1] = j+1;  pixelIndex[2] = k;
    dist[4]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i-1;  pixelIndex[1] = j-1;  pixelIndex[2] = k-1;
    dist[5]=outputImage->GetPixel(pixelIndex)+d3;
    pixelIndex[0] = i-1;  pixelIndex[1] = j;    pixelIndex[2] = k-1;
    dist[6]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i-1;  pixelIndex[1] = j+1;  pixelIndex[2] = k-1;
    dist[7]=outputImage->GetPixel(pixelIndex)+d3;
    pixelIndex[0] = i;    pixelIndex[1] = j-1;  pixelIndex[2] = k-1;
    dist[8]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k-1;
    dist[9]=outputImage->GetPixel(pixelIndex)+d1;
    pixelIndex[0] = i;    pixelIndex[1] = j+1;  pixelIndex[2] = k-1;
    dist[10]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i+1;  pixelIndex[1] = j-1;  pixelIndex[2] = k-1;
    dist[11]=outputImage->GetPixel(pixelIndex)+d3;
    pixelIndex[0] = i+1;  pixelIndex[1] = j;    pixelIndex[2] = k-1;
    dist[12]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i+1;  pixelIndex[1] = j+1;  pixelIndex[2] = k-1;
    dist[13]=outputImage->GetPixel(pixelIndex)+d3;

    min=dist[0];
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k;

    for(m=1;m<14;m++)
      if(dist[m]<min)
	{
	  min=dist[m];
          outputImage->SetPixel(pixelIndex, dist[m]);
	}


  }

  std::cout<<"------ Backward pass ------\n";

  for(uint k=size[2]-2; k>0; k--)
  for(uint j=size[1]-2; j>0; j--)
  for(uint i=size[0]-2; i>0; i--)
  {
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k;
    dist[0]=outputImage->GetPixel(pixelIndex);
    pixelIndex[0] = i;    pixelIndex[1] = j+1;  pixelIndex[2] = k;
    dist[1]=outputImage->GetPixel(pixelIndex)+d1;
    pixelIndex[0] = i+1;  pixelIndex[1] = j+1;  pixelIndex[2] = k;
    dist[2]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i+1;  pixelIndex[1] = j;    pixelIndex[2] = k;
    dist[3]=outputImage->GetPixel(pixelIndex)+d1;
    pixelIndex[0] = i+1;  pixelIndex[1] = j-1;  pixelIndex[2] = k;
    dist[4]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i-1;  pixelIndex[1] = j-1;  pixelIndex[2] = k+1;
    dist[5]=outputImage->GetPixel(pixelIndex)+d3;
    pixelIndex[0] = i-1;  pixelIndex[1] = j;    pixelIndex[2] = k+1;
    dist[6]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i-1;  pixelIndex[1] = j+1;  pixelIndex[2] = k+1;
    dist[7]=outputImage->GetPixel(pixelIndex)+d3;
    pixelIndex[0] = i;    pixelIndex[1] = j-1;  pixelIndex[2] = k+1;
    dist[8]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k+1;
    dist[9]=outputImage->GetPixel(pixelIndex)+d1;
    pixelIndex[0] = i;    pixelIndex[1] = j+1;  pixelIndex[2] = k+1;
    dist[10]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i+1;  pixelIndex[1] = j-1;  pixelIndex[2] = k+1;
    dist[11]=outputImage->GetPixel(pixelIndex)+d3;
    pixelIndex[0] = i+1;  pixelIndex[1] = j;    pixelIndex[2] = k+1;
    dist[12]=outputImage->GetPixel(pixelIndex)+d2;
    pixelIndex[0] = i+1;  pixelIndex[1] = j+1;  pixelIndex[2] = k+1;
    dist[13]=outputImage->GetPixel(pixelIndex)+d3;

    min=dist[0];
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k;
    outputImage->SetPixel(pixelIndex, dist[0]);

    for(m=1;m<14;m++)
      if(dist[m]<min)
	{
	  min=dist[m];
          outputImage->SetPixel(pixelIndex, dist[m]);
	}


  }
  for(uint j=0; j<size[1]; j++)
  for(uint i=0; i<size[0]; i++)
  {
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = 1;
    short tmp = outputImage->GetPixel(pixelIndex);
    pixelIndex[2] = 0;
    outputImage->SetPixel(pixelIndex, tmp);
    pixelIndex[2] = size[2]-2;
    tmp = outputImage->GetPixel(pixelIndex);
    pixelIndex[2] = size[2]-1;
    outputImage->SetPixel(pixelIndex, tmp);
  }
  for(uint k=0; k<size[2]; k++)
  for(uint i=0; i<size[0]; i++)
  {
    pixelIndex[0] = i;    pixelIndex[1] = 1;    pixelIndex[2] = k;
    short tmp = outputImage->GetPixel(pixelIndex);
    pixelIndex[1] = 0;
    outputImage->SetPixel(pixelIndex, tmp);
    pixelIndex[1] = size[1]-2;
    tmp = outputImage->GetPixel(pixelIndex);
    pixelIndex[1] = size[1]-1;
    outputImage->SetPixel(pixelIndex, tmp);
  }
  for(uint k=0; k<size[2]; k++)
  for(uint j=0; j<size[1]; j++)
  {
    pixelIndex[0] = 1;    pixelIndex[1] = j;    pixelIndex[2] = k;
    short tmp = outputImage->GetPixel(pixelIndex);
    pixelIndex[0] = 0;
    outputImage->SetPixel(pixelIndex, tmp);
    pixelIndex[0] = size[0]-2;
    tmp = outputImage->GetPixel(pixelIndex);
    pixelIndex[0] = size[0]-1;
    outputImage->SetPixel(pixelIndex, tmp);
  }
  for(uint k=0; k<size[2]; k++)
  for(uint j=0; j<size[1]; j++)
  for(uint i=0; i<size[0]; i++)
  {
    pixelIndex[0] = i;    pixelIndex[1] = j;    pixelIndex[2] = k;
    short tmp = outputImage->GetPixel(pixelIndex);
    outputImage->SetPixel(pixelIndex, tmp/3);
  }


}



