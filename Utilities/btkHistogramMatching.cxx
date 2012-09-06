/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/08/2012
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
#include <itkHistogramMatchingImageFilter.h>

/*Btk includes*/
#include "btkImageHelper.h"


int main (int argc, char* argv[])
{
  
  //This is an ITK filter
  //Purpose: Normalize the grayscale values between two images by histogram matching.
  //HistogramMatchingImageFilter normalizes the grayscale values of a source image based on the grayscale values of a reference image. This filter uses a histogram matching technique where the histograms of the two images are matched only at a specified number of quantile values.
  //This filter was originally designed to normalize MR images of the same MR protocol and same body part. The algorithm works best if background pixels are excluded from both the source and reference histograms. A simple background exclusion method is to exclude all pixels whose grayscale values are smaller than the mean grayscale value. ThresholdAtMeanIntensityOn() switches on this simple background exclusion method.
  //The source image can be set via either SetInput() or SetSourceImage(). The reference image can be set via SetReferenceImage().
  //SetNumberOfHistogramLevels() sets the number of bins used when creating histograms of the source and reference images. SetNumberOfMatchPoints() governs the number of quantile values to be matched.
  //This filter assumes that both the source and reference are of the same type and that the input and output image type have the same number of dimension and have scalar pixel types.
  //REFERENCE
  //Laszlo G. Nyul, Jayaram K. Udupa, and Xuan Zhang, "New Variants of a Method of MRI Scale Standardization", IEEE Transactions on Medical Imaging, 19(2):143-150, 2000.
  
  
  try {  

  TCLAP::CmdLine cmd("Normalize the grayscale values between two images by histogram matching.", ' ', "1.0", true);
  TCLAP::ValueArg<std::string> inputImageArg("i","input_file","input anatomical image file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> refImageArg("r","reference_file","reference anatomical image file",true,"","string", cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output anatomical image file",true,"","string", cmd);
  TCLAP::ValueArg< unsigned int > binArg("b","bin","number of histogram bins (default is 1024)",false,1024,"unsigned int",cmd);
  TCLAP::ValueArg< unsigned int > matchArg("p","points","number of match points (default is 7)",false,7,"unsigned int",cmd);
  TCLAP::ValueArg< bool > thresholdArg("t","threshold","threshold at mean intensity (default is true, meaning excluding voxels of the background)",false,true,"bool",cmd);


  // Parse the args.
  cmd.parse( argc, argv );
  
  std::string input_file                 = inputImageArg.getValue();
  std::string ref_file                   = refImageArg.getValue();
  std::string output_file                = outputImageArg.getValue();
  unsigned int numberOfBins              = binArg.getValue();
  unsigned int numberOfMatchPoints       = matchArg.getValue();
  bool threshold                         = thresholdArg.getValue();
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType;
  typedef ImageType::Pointer                    ImagePointer;
  typedef itk::HistogramMatchingImageFilter<ImageType, ImageType, PixelType> FilterType;

  std::cout<<"Read input images \n";
  ImagePointer inputImage = btk::ImageHelper<ImageType>::ReadImage(input_file);
  ImagePointer refImage   = btk::ImageHelper<ImageType>::ReadImage(ref_file);

  std::cout<<"Performing histogram matching\n";
  std::cout<<"number of bins : "<<numberOfBins<<std::endl;
  std::cout<<"number of match points : "<<numberOfMatchPoints<<std::endl;
  std::cout<<"threshold at mean intensity : "<<threshold<<std::endl;

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(inputImage);
  filter->SetReferenceImage(refImage);
  filter->SetNumberOfHistogramLevels(numberOfBins);
  filter->SetNumberOfMatchPoints(numberOfMatchPoints);
  filter->SetThresholdAtMeanIntensity(threshold);
  filter->Update();

  ImagePointer outputImage = filter->GetOutput();

  
  btk::ImageHelper<ImageType>::WriteImage(outputImage, output_file);
    
  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return 1;
}
