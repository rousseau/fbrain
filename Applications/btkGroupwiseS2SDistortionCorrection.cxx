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


int main( int argc, char *argv[] )
{

  try {

  const char *input = NULL, *output = NULL, *mgrad = NULL, *folder = NULL;
  const char *bvec = NULL, *cvec = NULL;

  int x1, y1, z1, x2, y2, z2;

  TCLAP::CmdLine cmd("Correct distortions caused by eddy currents in dwi sequences", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Original sequence",true,"homer","string");
  TCLAP::ValueArg<std::string> outputArg("o","output","Corrected sequence",true,"homer","string");
  TCLAP::ValueArg<std::string> folderArg("f","folder","Folder for transformatios",true,"homer","string");
  TCLAP::ValueArg<std::string> mgradArg("m","mgrad","Mean gradient",false,"none","string");

  TCLAP::ValueArg<std::string> bvecArg("g","bvec","Gradient directions",true,"","string");
  TCLAP::ValueArg<std::string> cvecArg("c","cvec","Corrected directions",true,"","string");

  TCLAP::ValueArg<int> x1Arg("","x1","ROI's min x value",false,0,"int");
  TCLAP::ValueArg<int> y1Arg("","y1","ROI's min y value",false,0,"int");
  TCLAP::ValueArg<int> z1Arg("","z1","ROI's min z value",false,0,"int");
  TCLAP::ValueArg<int> x2Arg("","x2","ROI's max x value",false,0,"int");
  TCLAP::ValueArg<int> y2Arg("","y2","ROI's max x value",false,0,"int");
  TCLAP::ValueArg<int> z2Arg("","z2","ROI's max x value",false,0,"int");


  // Add the argument nameArg to the CmdLine object. The CmdLine object
  // uses this Arg to parse the command line.
  cmd.add( z2Arg );
  cmd.add( y2Arg );
  cmd.add( x2Arg );
  cmd.add( z1Arg );
  cmd.add( y1Arg );
  cmd.add( x1Arg );

  cmd.add( outputArg );
  cmd.add( inputArg );
  cmd.add( folderArg );
  cmd.add( mgradArg );

  cmd.add( bvecArg );
  cmd.add( cvecArg );

  // Parse the argv array.
  cmd.parse( argc, argv );

  input  = inputArg.getValue().c_str();
  output = outputArg.getValue().c_str();
  mgrad = mgradArg.getValue().c_str();
  folder = folderArg.getValue().c_str();

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

  // ROI

  x1 = x1Arg.getValue();
  y1 = y1Arg.getValue();
  z1 = z1Arg.getValue();

  x2 = x2Arg.getValue();
  y2 = y2Arg.getValue();
  z2 = z2Arg.getValue();

  // Create fixed image region

  typedef itk::Image< PixelType, 3 >   ImageType;

  ImageType::RegionType fixedImageRegion;
  ImageType::IndexType  fixedImageRegionIndex;
  ImageType::SizeType   fixedImageRegionSize;

  fixedImageRegionIndex[0] = x1; fixedImageRegionIndex[1]= y1; fixedImageRegionIndex[2]= z1;
  fixedImageRegionSize[0]  = x2 - x1 + 1; fixedImageRegionSize[1] = y2 - y1 + 1; fixedImageRegionSize[2] = z2 - z1 + 1;

  fixedImageRegion.SetIndex(fixedImageRegionIndex);
  fixedImageRegion.SetSize(fixedImageRegionSize);

  // Distortion correction

  typedef btk::GroupwiseS2SDistortionCorrection< SequenceType >  FilterType;
  FilterType::Pointer filter = FilterType::New();

  filter -> SetInput( sequence );
  filter -> SetFixedImageRegion( fixedImageRegion );
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

  // Write corrected gradient table
//  The following line is commented since it has no sence for coupe à coupe corrections
//  filter -> WriteGradientTable( cvec );

  // Write transformations

  filter -> WriteTransforms( folder );

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

