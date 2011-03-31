/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/03/2010
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

#include "CmdLine.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "btkGroupwiseS2SDistortionCorrection.h"
#include "itkImageMaskSpatialObject.h"


int main( int argc, char *argv[] )
{

  try {

  const char *input = NULL, *output = NULL, *mgrad = NULL, *folder = NULL, *maskFile = NULL;
  const char *bvec = NULL, *cvec = NULL;

  TCLAP::CmdLine cmd("Correct distortions caused by eddy currents in dwi sequences", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Original sequence",true,"homer","string",cmd);
  TCLAP::ValueArg<std::string> bvecArg("g","bvec","Gradient directions",true,"","string",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","mask","Mask in the B0 image",true,"","string",cmd);

  TCLAP::ValueArg<std::string> outputArg("o","output","Corrected sequence",true,"homer","string",cmd);
  TCLAP::ValueArg<std::string> cvecArg("c","cvec","Corrected directions",true,"","string",cmd);

  TCLAP::ValueArg<std::string> folderArg("f","folder","Folder for transformatios",true,"homer","string",cmd);
  TCLAP::ValueArg<std::string> mgradArg("","meanGradient","Mean gradient",false,"none","string",cmd);

  // Parse the argv array.
  cmd.parse( argc, argv );

  input  = inputArg.getValue().c_str();
  output = outputArg.getValue().c_str();
  mgrad = mgradArg.getValue().c_str();
  folder = folderArg.getValue().c_str();
  maskFile = maskArg.getValue().c_str();

  bvec = bvecArg.getValue().c_str();
  cvec = cvecArg.getValue().c_str();

  // Read sequence

  typedef short         PixelType;
  const   unsigned int  Dimension = 4;

  typedef itk::Image< PixelType, Dimension >   SequenceType;

  typedef itk::ImageFileReader< SequenceType >  SequenceReaderType;
  SequenceReaderType::Pointer sequenceReader = SequenceReaderType::New();
  sequenceReader -> SetFileName( input );
  sequenceReader -> Update();

  SequenceType::Pointer sequence = sequenceReader -> GetOutput();


  // Read image mask and create spatial object

  typedef itk::Image< unsigned char, 3 >   ImageMaskType;
  typedef itk::ImageFileReader< ImageMaskType >  ImageMaskReaderType;

  ImageMaskReaderType::Pointer imageMaskReader = ImageMaskReaderType::New();
  imageMaskReader -> SetFileName( maskFile );
  imageMaskReader -> Update();

  ImageMaskType::Pointer imageMask = imageMaskReader -> GetOutput();

  typedef itk::ImageMaskSpatialObject< 3 >   MaskType;
  MaskType::Pointer mask = MaskType::New();

  mask -> SetImage( imageMask );


  // Distortion correction

  typedef btk::GroupwiseS2SDistortionCorrection< SequenceType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter -> SetInput( sequence );
  filter -> SetFixedImageRegion( mask -> GetAxisAlignedBoundingBoxRegion() );
  filter -> SetGradientTable( bvec );
  filter -> StartRegistration( );

  SequenceType::Pointer correctedSequence = filter -> GetOutput();

  SequenceType::SizeType sequenceSize = sequence -> GetLargestPossibleRegion().GetSize();

  // Write corrected sequence

  typedef itk::ImageFileWriter< SequenceType >  WriterType;

  WriterType::Pointer writer =  WriterType::New();
  writer->SetFileName( output );
  writer->SetInput( correctedSequence );
  writer->Update();

  // Write transformations

  filter -> WriteTransforms( folder );

  typedef itk::Image< PixelType, 3 >   ImageType;
  typedef itk::ImageFileWriter< ImageType >  ImageWriterType;
  ImageWriterType::Pointer imageWriter =  ImageWriterType::New();

  if ( strcmp(mgrad,"none") != 0 )
  {
    imageWriter->SetFileName( mgrad );
    imageWriter->SetInput( filter -> GetMeanGradient() );
    imageWriter->Update();
  }


  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

