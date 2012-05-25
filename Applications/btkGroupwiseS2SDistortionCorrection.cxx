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

/* Standard includes */
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"

/* Btk includes */
#include "btkFileNameTools.h"
#include "btkGroupwiseS2SDistortionCorrection.h"


int main( int argc, char *argv[] )
{

try {

  std::string input, bvec, bval, maskFile;
  std::string output, bvec_out, bval_out, mgrad;
  std::string folder, inRadix, outRadix;

  TCLAP::CmdLine cmd("Correct distortions caused by eddy currents in dwi "
      "sequences", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Original sequence",true,"",
      "string",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","mask","Mask in the B0 image",true,"",
      "string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Corrected sequence (spatial resampling, not in the Q-space)",true,"",
      "string",cmd);
  TCLAP::ValueArg<std::string> folderArg("f","transformation_folder","Folder for"
      " transformations (main output of the algorithm)",true,"","string",cmd);

  TCLAP::ValueArg<std::string> mgradArg("","mean-gradient","Mean gradient (output of the algorithm)",false,
      "","string",cmd);

  // Parse the argv array.
  cmd.parse( argc, argv );

  input = inputArg.getValue();
  output = outputArg.getValue();
  mgrad = mgradArg.getValue();
  folder = folderArg.getValue();
  maskFile = maskArg.getValue();

  std::cout << folder.c_str() << std::endl;

  inRadix = btk::GetRadixOf(input);
  bvec = inRadix + ".bvec";
  bval = inRadix + ".bval";

  outRadix = btk::GetRadixOf(output);
  bvec_out = outRadix + ".bvec";
  bval_out = outRadix + ".bval";

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

  // Distortion correction -- (main part of the program)

  typedef btk::GroupwiseS2SDistortionCorrection< SequenceType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter -> SetInput( sequence );
  filter -> SetFixedImageRegion( mask -> GetAxisAlignedBoundingBoxRegion() );
  filter -> SetGradientTable( bvec.c_str() );
  filter -> StartRegistration( );

  SequenceType::Pointer correctedSequence = filter -> GetOutput();

  SequenceType::SizeType sequenceSize = sequence -> GetLargestPossibleRegion().GetSize();

  // Write corrected sequence

  typedef itk::ImageFileWriter< SequenceType >  WriterType;

  WriterType::Pointer writer =  WriterType::New();
  writer->SetFileName( output );
  writer->SetInput( correctedSequence );
  writer->Update();

  // Write gradient table
  char clcopybvec[255];
  sprintf(clcopybvec,"cp %s %s",bvec.c_str(),bvec_out.c_str());
  system(clcopybvec);

  // Write b-values
  char clcopybval[255];
  sprintf(clcopybval,"cp %s %s",bval.c_str(),bval_out.c_str());
  system(clcopybval);

  // Write transformations

  filter -> WriteTransforms( folder.c_str() );

  typedef itk::Image< PixelType, 3 >   ImageType;
  typedef itk::ImageFileWriter< ImageType >  ImageWriterType;
  ImageWriterType::Pointer imageWriter =  ImageWriterType::New();

  if ( strcmp(mgrad.c_str(),"") != 0 )
  {
    imageWriter->SetFileName( mgrad );
    imageWriter->SetInput( filter -> GetMeanGradient() );
    imageWriter->Update();
  }

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

