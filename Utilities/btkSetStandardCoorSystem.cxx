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

#include "itkOrientImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkJoinSeriesImageFilter.h"

#include "btkDiffusionGradientTable.h"

#include "itkEuler3DTransform.h"

int main( int argc, char *argv[] )
{

  try {

  const char *inputName = NULL, *outputName = NULL, *gTableFile = NULL, *cTableFile = NULL;
  unsigned int dim;

  TCLAP::CmdLine cmd("Sets the direction to the identity, and the origin to the center of the image",
  ' ', "Unversioned");

  TCLAP::ValueArg<std::string>  inputArg("i","input","Input image",true,"","string",cmd);
  TCLAP::ValueArg<std::string>  outputArg("o","output","Output folder",true,"","string",cmd);
  TCLAP::ValueArg<unsigned int> dimArg("d","dimension","Image dimension (3 / 4)",true,3,"unsigned int",cmd);
  TCLAP::ValueArg<std::string>  gTableArg("g","bvec","Gradient table",false,"","string",cmd);
  TCLAP::ValueArg<std::string>  cTableArg("c","ctable","Corrected table",false,"","string",cmd);

  // Parse the argv array.
  cmd.parse( argc, argv );

  inputName  = inputArg.getValue().c_str();
  outputName = outputArg.getValue().c_str();
  gTableFile = gTableArg.getValue().c_str();
  cTableFile = cTableArg.getValue().c_str();
  dim        = dimArg.getValue();

  typedef  short  PixelType;
  typedef itk::Image< PixelType, 3 >  ImageType;
  typedef itk::OrientImageFilter< ImageType, ImageType >  FilterType;


  if (dim==4)
  {
    typedef itk::Image< PixelType, 4 >  SequenceType;
    typedef SequenceType::PointType  OriginType;

    typedef itk::ImageFileReader< SequenceType  > ImageReaderType;
    ImageReaderType::Pointer  sequenceReader  = ImageReaderType::New();
    sequenceReader -> SetFileName( inputName );
    sequenceReader -> Update();

    SequenceType::Pointer     sequence = sequenceReader -> GetOutput();
    SequenceType::SpacingType spacing  = sequence -> GetSpacing();
    SequenceType::SizeType    size     = sequence -> GetLargestPossibleRegion().GetSize();

    SequenceType::DirectionType sequenceDirection = sequence -> GetDirection();
    vnl_matrix<double> originalSpatialDirectionVnl(3,3);
    originalSpatialDirectionVnl = sequenceDirection.GetVnlMatrix().extract(3,3);

    typedef itk::ExtractImageFilter< SequenceType, ImageType > ExtractorType;
    ExtractorType::Pointer extractor  = ExtractorType::New();
    extractor -> SetInput( sequence );

    SequenceType::RegionType inputRegion = sequence -> GetLargestPossibleRegion();

    unsigned int numberOfFrames = size[3];
    size[3] = 0;

    SequenceType::IndexType start = inputRegion.GetIndex();

    typedef itk::JoinSeriesImageFilter< ImageType, SequenceType > JoinerType;
    JoinerType::Pointer joiner = JoinerType::New();
    joiner -> SetOrigin( 0.0 );
    joiner -> SetSpacing( spacing[3] );

    for (unsigned int i = 0; i < numberOfFrames; i++)
    {
      start[3] = i;

      SequenceType::RegionType desiredRegion;
      desiredRegion.SetSize(  size  );
      desiredRegion.SetIndex( start );

      extractor -> SetExtractionRegion( desiredRegion );

      FilterType::Pointer filter = FilterType::New();
      filter -> SetInput( extractor -> GetOutput()  );
      filter -> UseImageDirectionOn();
      filter -> SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
      filter -> Update();

/*      std::cout << "Extractor direction = " << std::endl;
      std::cout << extractor -> GetOutput() -> GetDirection() << std::endl;

      std::cout << "Reoriented direction = " << std::endl;
      std::cout << filter -> GetOutput() -> GetDirection() << std::endl; */

      joiner -> SetInput( i, filter -> GetOutput() );
    }

    joiner -> Update();

    sequence = joiner -> GetOutput();
    spacing  = sequence -> GetSpacing();
    size     = sequence -> GetLargestPossibleRegion().GetSize();
    sequenceDirection = sequence -> GetDirection();

    // Change gradient table

    vnl_matrix<double> reorientedSpatialDirectionVnl(3,3);
    reorientedSpatialDirectionVnl = sequenceDirection.GetVnlMatrix().extract(3,3);

    vnl_matrix<double> rotationMatrix;
    rotationMatrix = reorientedSpatialDirectionVnl.transpose()*originalSpatialDirectionVnl;

    typedef btk::DiffusionGradientTable< SequenceType > GradientTableType;
    GradientTableType::Pointer gradientTable = GradientTableType::New();

    typedef itk::Euler3DTransform<double> TransformType;
    TransformType::Pointer transform = TransformType::New();

    gradientTable -> SetNumberOfGradients(numberOfFrames);
    gradientTable -> SetImage( sequence );
//    gradientTable -> SetTransform( transform );
    gradientTable -> SetRotationMatrix( rotationMatrix );
    gradientTable -> LoadFromFile( gTableFile);
    gradientTable -> RotateGradients();
    gradientTable -> SaveToFile( cTableFile);

    // Change sequence

    sequenceDirection.SetIdentity();
    sequence -> SetDirection (sequenceDirection);

    SequenceType::PointType origin;

    for (unsigned int i=0; i<3; i++)
    {
      origin[i] = (1.0-size[i])*0.5*spacing[i];
    }

    sequence -> SetOrigin( origin );

    typedef itk::ImageFileWriter< SequenceType >  WriterType;

    WriterType::Pointer writer =  WriterType::New();
    writer->SetFileName( outputName );
    writer->SetInput( sequence );
    writer->Update();


  } else if (dim == 3)
  {
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

    FilterType::Pointer filter = FilterType::New();

    filter -> SetInput( image );
    filter -> UseImageDirectionOn();
    filter -> SetDesiredCoordinateOrientation( itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI );
    filter -> Update();

    image = filter -> GetOutput();
    image -> SetDirection (imageDirection);

    ImageType::SpacingType spacing = image -> GetSpacing();
    ImageType::SizeType    size    = image -> GetLargestPossibleRegion().GetSize();

    ImageType::PointType origin;

    for (unsigned int i=0; i<3; i++)
    {
      origin[i] = (1.0-size[i])*0.5*spacing[i];
    }

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

