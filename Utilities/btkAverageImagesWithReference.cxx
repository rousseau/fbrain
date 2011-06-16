/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 21/12/2010
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

#include "itkImage.h"
#include <tclap/CmdLine.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


int main( int argc, char *argv[] )
{

  try {

  std::vector< std::string > input;
  const char *outImage = NULL, *refImageName = NULL;

  TCLAP::CmdLine cmd("Creates a high resolution image from a set of low resolution images", ' ', "Unversioned");

  TCLAP::MultiArg<std::string> inputArg("i","input","image file",true,"string",cmd);
  TCLAP::ValueArg<std::string> outArg("o","output","High resolution image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> refArg("r","reference","Reference image",true,"","string",cmd);
  TCLAP::SwitchArg nnSwitch("","nn","Nearest Neighbor interpolation", cmd, false);

  // Parse the argv array.
  cmd.parse( argc, argv );

  input = inputArg.getValue();
  outImage = outArg.getValue().c_str();
  refImageName = refArg.getValue().c_str();

  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;
  typedef itk::ImageFileReader< ImageType  >  ImageReaderType;

  unsigned int numberOfImages = input.size();

  std::vector< ImagePointer >	imageArray;

   // Read images

  for (unsigned int i = 0; i<numberOfImages; i++)
  {
    ImageReaderType::Pointer  imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input[i].c_str() );
    imageReader -> Update();
    imageArray.push_back( imageReader -> GetOutput() );

  }

  ImageReaderType::Pointer refImageReader = ImageReaderType::New();
  refImageReader -> SetFileName( refImageName );
  refImageReader->Update();
  ImageType::Pointer refImage = refImageReader -> GetOutput();

  // Compute average

  typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
  IteratorType refIt( refImage, refImage -> GetLargestPossibleRegion() );

  typedef itk::LinearInterpolateImageFunction<ImageType,double> InterpolatorType;
  typedef itk::NearestNeighborInterpolateImageFunction< ImageType >  NNInterpolatorType;

  NNInterpolatorType::Pointer nn_interpolator = NNInterpolatorType::New();
  InterpolatorType::Pointer interpolator = InterpolatorType::New();

  ImageType::IndexType index;
  ImageType::PointType point;
  double value;
  unsigned int counter;

  for (refIt.GoToBegin(); !refIt.IsAtEnd(); ++refIt )
  {
    index = refIt.GetIndex();
    refImage -> TransformIndexToPhysicalPoint( index, point );
    value = 0;
    counter = 0;

    for (unsigned int i=0; i < numberOfImages; i++)
    {
      if (nnSwitch.isSet())
      {
        nn_interpolator -> SetInputImage( imageArray[i] );
        if ( nn_interpolator -> IsInsideBuffer( point ) )
        {
          value += nn_interpolator -> Evaluate( point );
          counter++;
        }
      }
      else
      {
        interpolator -> SetInputImage( imageArray[i] );
        if ( interpolator -> IsInsideBuffer( point ) )
        {
          value += interpolator -> Evaluate( point );
          counter++;
        }
      }
    }

    if ( counter>0 ) refIt.Set(value/counter);
    else refIt.Set(0.0);
  }

  // Write average

  typedef itk::ImageFileWriter< ImageType >  WriterType;

  WriterType::Pointer writer =  WriterType::New();
  writer->SetFileName( outImage );
  writer->SetInput( refImage );
  writer->Update();


  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

