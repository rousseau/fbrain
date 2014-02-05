/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 13/12/2013
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


/*  Standard includes */
#include "vector"
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkMultiLabelSTAPLEImageFilter.h"


//ITK Implementation of : 
//T. Rohlfing, D. B. Russakoff, and C. R. Maurer, Jr., "Performance-based classifier combination in atlas-based image segmentation using expectation-maximization parameter estimation," IEEE Transactions on Medical Imaging, vol. 23, pp. 983-994, Aug. 2004.


int main( int argc, char *argv[] )
{
  try {
    
    TCLAP::CmdLine cmd("Apply the multi-label STAPLE algorithm implemented in ITK", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg ("i","input","Label image file.", true,"string",cmd);
    TCLAP::ValueArg<std::string> outArg   ("o","output","Output file of the estimated label image.", true,"","string",cmd);
    
    // Parse the argv array.
    cmd.parse( argc, argv );
    
    // Get the value parsed by each arg. 
    std::vector<std::string> label_file  = inputArg.getValue();
    std::string output_file              = outArg.getValue();

  //ITK declaration
  const   unsigned int        Dimension = 3;
  typedef itk::Image< short, Dimension >    ShortImageType;

  typedef itk::MultiLabelSTAPLEImageFilter< ShortImageType >  STAPLEFilterType;
  typedef STAPLEFilterType::Pointer STAPLEFilterPointer;

  typedef ShortImageType::Pointer ShortImagePointer;

  typedef itk::ImageFileReader< ShortImageType >  ShortReaderType;
  typedef itk::ImageFileWriter< ShortImageType >  ShortWriterType;


  //STAPLEFilterPointer filter = STAPLEFilterType::New();

  /*
  filter->SetMaximumNumberOfIterations( 1 );

  std::vector<ShortReaderType::Pointer> readers;
  for ( int i = 0; i < label_file.size(); ++i )
    {
    ShortReaderType::Pointer reader = ShortReaderType::New();
    reader->SetFileName( label_file[i] );
    filter->SetInput( i, reader->GetOutput() );
    readers.push_back( reader );
    }

  ShortWriterType::Pointer writer = ShortWriterType::New();
  writer->SetFileName( output_file );
  writer->SetInput( filter->GetOutput() );

  writer->Update();
*/

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}
