/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/11/2010
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
#include "itkImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


int main( int argc, char *argv[] )
{

  try {

  const char *inputName = NULL, *outputName = NULL;
  unsigned int dim;

  TCLAP::CmdLine cmd("Sets the direction to the identity, and the origin to the center of the image",
  ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"none","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Output folder",true,"none","string",cmd);
  TCLAP::ValueArg<unsigned int> dimArg("d","dimension","Image dimension (3 / 4)",true,3,"string",cmd);

  // Parse the argv array.
  cmd.parse( argc, argv );

  inputName  = inputArg.getValue().c_str();
  outputName = outputArg.getValue().c_str();
  dim        = dimArg.getValue();

  typedef  short  PixelType;

  if (dim==4)
  {
    typedef itk::Image< PixelType, 4 >  ImageType;
    typedef ImageType::PointType  OriginType;

    typedef itk::ImageFileReader< ImageType  > ImageReaderType;
    ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
    imageReader -> SetFileName( inputName );
    imageReader -> Update();

    ImageType::Pointer image = imageReader -> GetOutput();

    ImageType::DirectionType imageDirection = image -> GetDirection();
    vnl_matrix<double> imageDirectionVnl(4,4);
    imageDirectionVnl = imageDirection.GetVnlMatrix();

    vnl_matrix<double> imageSpatialDirectionVnl(3,3);
    imageSpatialDirectionVnl = imageDirectionVnl.extract(3,3);

    imageDirection.SetIdentity();

    if ( vnl_determinant( imageSpatialDirectionVnl) < 0 )
      imageDirection(2,2) = -1;

    image -> SetDirection (imageDirection);

    ImageType::SpacingType spacing = image -> GetSpacing();
    ImageType::SizeType    size    = image -> GetLargestPossibleRegion().GetSize();

    ImageType::PointType origin;

    for (unsigned int i=0; i<3; i++)
    {
      origin[i] = (1.0-size[i])*0.5*spacing[i];
    }

    if ( vnl_determinant( imageSpatialDirectionVnl) < 0 )
      origin[2] = -origin[2];


    image -> SetOrigin( origin );

    typedef itk::ImageFileWriter< ImageType >  WriterType;

    WriterType::Pointer writer =  WriterType::New();
    writer->SetFileName( outputName );
    writer->SetInput( image );
    writer->Update();
  } else if (dim == 3)
  {
    typedef itk::Image< PixelType, 3 >  ImageType;
    typedef ImageType::PointType  OriginType;

    typedef itk::ImageFileReader< ImageType  > ImageReaderType;
    ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
    imageReader -> SetFileName( inputName );
    imageReader -> Update();

    ImageType::Pointer image = imageReader -> GetOutput();

    ImageType::DirectionType imageDirection = image -> GetDirection();
    vnl_matrix<double> imageSpatialDirectionVnl(3,3);
    imageSpatialDirectionVnl = imageDirection.GetVnlMatrix();

    imageDirection.SetIdentity();

    // FIXME To do things in a more general (correct) way, we should detect the
    // orientation and then do a voxel reorganization. Or just force an image
    // orientation with det = 1. This can be done with itk classes.
    // (itkOrientImageFilter)

    if ( vnl_determinant( imageSpatialDirectionVnl) < 0 )
      imageDirection(2,2) = -1;

    image -> SetDirection (imageDirection);

    ImageType::SpacingType spacing = image -> GetSpacing();
    ImageType::SizeType    size    = image -> GetLargestPossibleRegion().GetSize();

    ImageType::PointType origin;

    for (unsigned int i=0; i<3; i++)
    {
      origin[i] = (1.0-size[i])*0.5*spacing[i];
    }

    if ( vnl_determinant( imageSpatialDirectionVnl) < 0 )
      origin[2] = -origin[2];


    image -> SetOrigin( origin );

    typedef itk::ImageFileWriter< ImageType >  WriterType;

    WriterType::Pointer writer =  WriterType::New();
    writer->SetFileName( outputName );
    writer->SetInput( image );
    writer->Update();

  } else
  {
    std::cout << "ERROR: Image dimension incorrect." << std::endl;
    return EXIT_FAILURE;

  }


  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}

