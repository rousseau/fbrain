/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 27/09/2011
  Author(s): Estanislao Oubel (oubel@unistra.fr)

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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

/* Standard includes */
#include <tclap/CmdLine.h>
#include "stdio.h"

/* Itk includes */
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageFileReader.h"
#include "itkEuler3DTransform.h"

/* BTK includes */
#include "btkImageHelper.h"

int main( int argc, char *argv[] )
{
  try {

  std::string imageAFile;
  std::string imageBFile;
  std::string maskFile;

  // Parse arguments

  TCLAP::CmdLine cmd("Calculates MSE, MI, and NC of two images A and B. The use "
      "of a mask is possible.", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> imageAArg("a","imageA","Image A",true,"","string",cmd);
  TCLAP::ValueArg<std::string> imageBArg("b","imageB","Image B",true,"","string",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","mask","Mask file on image A",false,"","string",cmd);

  TCLAP::SwitchArg ncSwitch("","nc","Calculates the normalized correlation.", false);
  TCLAP::SwitchArg miSwitch("","mi","Calculates the mutual information.", false);
  TCLAP::SwitchArg mseSwitch("","mse","Calculates the mean squared error.", false);
  TCLAP::SwitchArg nmiSwitch("","nmi","Calculates the normalized mutual information.", false);

  std::vector<TCLAP::Arg*>  xorlist;
  xorlist.push_back(&ncSwitch);
  xorlist.push_back(&miSwitch);
  xorlist.push_back(&mseSwitch);
  xorlist.push_back(&nmiSwitch);

  cmd.xorAdd( xorlist );

  // Parse the argv array.
  cmd.parse( argc, argv );

  imageAFile = imageAArg.getValue();
  imageBFile = imageBArg.getValue();
  maskFile 	 = maskArg.getValue();

  bool nc_switch  = ncSwitch.getValue();
  bool mi_switch  = miSwitch.getValue();
  bool mse_switch = mseSwitch.getValue();
  bool nmi_switch = nmiSwitch.getValue();

  // typedefs

  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;

  typedef ImageType::RegionType               RegionType;
  typedef std::vector< RegionType >           RegionArrayType;

  typedef itk::ImageFileReader< ImageType  >  ImageReaderType;

  typedef itk::Image< unsigned char, Dimension >    ImageMaskType;
  typedef ImageMaskType::Pointer                    ImageMaskPointer;

  typedef itk::ImageFileReader< ImageMaskType >     MaskReaderType;

  typedef itk::ImageMaskSpatialObject< Dimension >  MaskType;
  typedef MaskType::Pointer  MaskPointer;

  typedef itk::NormalizedCorrelationImageToImageMetric< ImageType,
                                                        ImageType > NCMetricType;

  typedef itk::MattesMutualInformationImageToImageMetric< ImageType,
                                                         ImageType > MIMetricType;

  typedef itk::MeanSquaresImageToImageMetric< ImageType,
                                              ImageType > MSEMetricType;

  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric< ImageType,
                                                         ImageType > NMIMetricType;
  // Read images

  ImagePointer imageA = btk::ImageHelper<ImageType>::ReadImage(imageAFile);
  ImagePointer imageB = btk::ImageHelper<ImageType>::ReadImage(imageBFile);


  // Read mask
  MaskPointer mask = MaskType::New();

  if (strcmp(maskFile.c_str(),"") != 0)
  {
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader -> SetFileName( maskFile.c_str() );
    maskReader -> Update();

    mask -> SetImage( maskReader -> GetOutput() );
  }

  // Compute similarity

  typedef itk::Euler3DTransform< double > EulerTransformType;
  EulerTransformType::Pointer identity = EulerTransformType::New();
  identity -> SetIdentity();

  typedef itk::LinearInterpolateImageFunction<
                                    ImageType,
                                    double>     InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  if (nc_switch)
  {
    NCMetricType::Pointer nc = NCMetricType::New();
    nc -> SetFixedImage(  imageA );
    nc -> SetMovingImage( imageB );

    if (strcmp(maskFile.c_str(),"") != 0)
    {
      nc -> SetFixedImageRegion( mask -> GetAxisAlignedBoundingBoxRegion() );
      nc -> SetFixedImageMask( mask );
    } else
      {
        nc -> SetFixedImageRegion( imageA -> GetLargestPossibleRegion() );
      }

    nc -> SetTransform( identity );
    nc -> SetInterpolator( interpolator );
    nc -> Initialize();

    float nc_value = -nc -> GetValue( identity -> GetParameters() );
    std::cout << "NC = " << nc_value << std::endl;

  }
  else if (mi_switch)
  {
    MIMetricType::Pointer mi = MIMetricType::New();
    mi -> SetFixedImage(  imageA );
    mi -> SetMovingImage( imageB );

    if (strcmp(maskFile.c_str(),"") != 0)
    {
      mi -> SetFixedImageRegion( mask -> GetAxisAlignedBoundingBoxRegion() );
      mi -> SetFixedImageMask( mask );
    } else
      {
        mi -> SetFixedImageRegion( imageA -> GetLargestPossibleRegion() );
      }

    mi -> SetTransform( identity );
    mi -> SetInterpolator( interpolator );
    mi -> Initialize();

    float mi_value = -mi -> GetValue( identity -> GetParameters() );
    std::cout << "MI = " << mi_value << std::endl;

   }
   else if (mse_switch)
    {
      MSEMetricType::Pointer mse = MSEMetricType::New();
      mse -> SetFixedImage(  imageA );
      mse -> SetMovingImage( imageB );

      if (strcmp(maskFile.c_str(),"") != 0)
      {
        mse -> SetFixedImageRegion( mask -> GetAxisAlignedBoundingBoxRegion() );
        mse -> SetFixedImageMask( mask );
      }
      else
      {
        mse -> SetFixedImageRegion( imageA -> GetLargestPossibleRegion() );
      }

      mse -> SetTransform( identity );
      mse -> SetInterpolator( interpolator );
      mse -> Initialize();

      float mse_value = mse -> GetValue( identity -> GetParameters() );
      std::cout << "MSE = " << mse_value << std::endl;
    }
  else if (nmi_switch)
   {
     NMIMetricType::Pointer nmi = NMIMetricType::New();
     nmi -> SetFixedImage(  imageA );
     nmi -> SetMovingImage( imageB );

     if (strcmp(maskFile.c_str(),"") != 0)
     {
       nmi -> SetFixedImageRegion( mask -> GetAxisAlignedBoundingBoxRegion() );
       nmi -> SetFixedImageMask( mask );
     }
     else
     {
       nmi -> SetFixedImageRegion( imageA -> GetLargestPossibleRegion() );
     }

     unsigned int nBins = 256;
     NMIMetricType::HistogramType::SizeType histSize;
     histSize.SetSize(2);
     histSize[0] = nBins;
     histSize[1] = nBins;
     nmi->SetHistogramSize(histSize);

     nmi -> SetTransform( identity );
     nmi -> SetInterpolator( interpolator );
     nmi -> Initialize();

     float nmi_value = nmi -> GetValue( identity -> GetParameters() );
     std::cout << "NMI = " << nmi_value << std::endl;
   }

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

