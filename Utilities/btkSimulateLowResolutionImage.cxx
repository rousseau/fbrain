/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 24/05/2011
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

#include "itkEuler3DTransform.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageMaskSpatialObject.h"

#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "btkSliceBySliceTransform.h"

#include "btkSuperResolutionImageFilter.h"

#include <tclap/CmdLine.h>

int main( int argc, char *argv[] )
{

  try {

  const char *inputFile = NULL;
  const char *maskFile  = NULL;
  const char *transformFile = NULL;
  const char *simFile  = NULL;
  const char *refFile  = NULL;

  // Parse arguments

  TCLAP::CmdLine cmd("Simulates a low resolution image from a high resolution "
      "image (reconstructed, super-resolution, or acquired) and a known "
      "transformation between both images. The use of either a boxcar or a Gaussian "
      "kernel (default) is possible.", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Low-resolution image file.",
      true,"","string",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","mask","Low-resolution image mask file.",
      false,"","string",cmd);
  TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image. "
      "Typically the output of btkImageReconstruction is used." ,true,"","string",cmd);
  TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",
      true,"","string",cmd);
  TCLAP::SwitchArg  boxcarSwitchArg("","boxcar","A boxcar-shaped PSF is assumed "
      "as imaging model (by default a Gaussian-shaped PSF is employed.).",false);
  TCLAP::ValueArg<std::string> transArg("t","transform","transform file",false,
      "","string",cmd);

  // Parse the argv array.
  cmd.parse( argc, argv );

  inputFile = inputArg.getValue().c_str();
  maskFile = maskArg.getValue().c_str();
  refFile = refArg.getValue().c_str();
  simFile = outArg.getValue().c_str();
  transformFile = transArg.getValue().c_str();

  // typedefs
  const   unsigned int    Dimension = 3;
  typedef btk::SliceBySliceTransform< double, Dimension > TransformType;
  itk::TransformFactory<TransformType>::RegisterTransform();

  typedef float  PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;
  typedef std::vector<ImagePointer>           ImagePointerArray;

  typedef itk::Image< unsigned char, Dimension >  ImageMaskType;
  typedef itk::ImageFileReader< ImageMaskType > MaskReaderType;
  typedef itk::ImageMaskSpatialObject< Dimension > MaskType;

  typedef ImageType::RegionType               RegionType;
  typedef std::vector< RegionType >           RegionArrayType;

  typedef itk::ImageFileReader< ImageType >   ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >   WriterType;

  typedef itk::TransformFileReader     TransformReaderType;
  typedef TransformReaderType::TransformListType * TransformListType;

  typedef btk::SuperResolutionImageFilter< ImageType, ImageType >  ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();

  // Filter setup

  ImageType::IndexType  roiIndex;
  ImageType::SizeType   roiSize;

  // add image
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader -> SetFileName( inputFile );
  imageReader -> Update();
  resampler   -> AddInput( imageReader -> GetOutput() );

  // add region
  if ( strcmp(maskFile,"") != 0 )
  {
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader -> SetFileName( maskFile );
    maskReader -> Update();

    MaskType::Pointer mask = MaskType::New();
    mask -> SetImage( maskReader -> GetOutput() );

    RegionType roi = mask -> GetAxisAlignedBoundingBoxRegion();
    roiIndex = roi.GetIndex();
    roiSize  = roi.GetSize();

  } else
      {
        roiSize  = imageReader -> GetOutput() -> GetLargestPossibleRegion().GetSize();
        roiIndex = imageReader -> GetOutput() -> GetLargestPossibleRegion().GetIndex();
      }

  RegionType imageRegion;
  imageRegion.SetIndex(roiIndex);
  imageRegion.SetSize (roiSize);
  resampler -> AddRegion( imageRegion );


  if ( strcmp(transformFile,"") != 0 )
  {
    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader -> SetFileName( transformFile );
    transformReader -> Update();

    TransformListType transforms = transformReader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator titr = transforms->begin();
    TransformType::Pointer trans = dynamic_cast< TransformType * >( titr->GetPointer() );

    for(unsigned int j=0; j< trans -> GetNumberOfSlices(); j++)
      resampler -> SetTransform(0, j, trans -> GetSliceTransform(j) ) ;
  }

  // Set reference image

  ImageReaderType::Pointer refReader = ImageReaderType::New();
  refReader -> SetFileName( refFile );
  refReader -> Update();

  resampler -> UseReferenceImageOn();
  resampler -> SetReferenceImage( refReader -> GetOutput() );

  if ( boxcarSwitchArg.isSet() )
    resampler -> SetPSF( ResamplerType::BOXCAR );
  resampler -> CreateH();

  // Write simulated image

  WriterType::Pointer writer =  WriterType::New();
  writer -> SetFileName( simFile );
  writer -> SetInput( resampler -> GetSimulatedImage(0) );

  if ( strcmp(simFile,"") != 0)
  {
    std::cout << "Writing " << simFile << " ... ";
    writer -> Update();
    std::cout << "done." << std::endl;
  }

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

