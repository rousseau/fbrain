/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 15/11/2010
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

#include "tclap/CmdLine.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageDuplicator.h"

int main( int argc, char *argv[] )
{

  try {

  const char *inputName = NULL, *outputName = NULL;
  int label;

  TCLAP::CmdLine cmd("Splits a label image into binary components", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"none","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Output folder",true,"none","string",cmd);
  TCLAP::ValueArg<int> labelArg("l","label","Label value",true,0,"int",cmd);

  // Parse the argv array.
  cmd.parse( argc, argv );

  inputName = inputArg.getValue().c_str();
  outputName = outputArg.getValue().c_str();
  label = labelArg.getValue();

  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;

  // Read input

  typedef itk::ImageFileReader< ImageType  >  ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader -> SetFileName( inputName );
  imageReader -> Update();
  ImagePointer input = imageReader -> GetOutput();

  // Duplicate input

  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator -> SetInputImage (input);
  duplicator -> Update();
  ImageType::Pointer labelImage = duplicator -> GetOutput();

  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  IteratorType labelIt( labelImage, labelImage -> GetLargestPossibleRegion() );

  for (labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt)
  {
    if ( labelIt.Get() == label )
      labelIt.Set(1);
    else
      labelIt.Set(0);
  }

  typedef itk::ImageFileWriter< ImageType >  ImageWriterType;
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( outputName );
  imageWriter->SetInput( labelImage );
  imageWriter->Update();

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

