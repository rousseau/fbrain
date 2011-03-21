/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 26/01/2010
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
#include "itkImageRegionIteratorWithIndex.h"

#include "itkImageFileReader.h"
#include "itkSpatialOrientationAdapter.h"
#include "itkExtractImageFilter.h"
#include "itkSpatialOrientationAdapter.h"

#include <iostream>

#include "CmdLine.h"

typedef itk::SpatialOrientation::ValidCoordinateOrientationFlags
SO_OrientationType;
std::string SO_OrientationToString(SO_OrientationType in)
{
  switch(in)
    {
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP:
      return std::string("RIP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP:
      return std::string("LIP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP:
      return std::string("RSP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP:
      return std::string("LSP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA:
      return std::string("RIA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA:
      return std::string("LIA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA:
      return std::string("RSA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA:
      return std::string("LSA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP:
      return std::string("IRP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP:
      return std::string("ILP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP:
      return std::string("SRP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP:
      return std::string("SLP");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA:
      return std::string("IRA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA:
      return std::string("ILA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA:
      return std::string("SRA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA:
      return std::string("SLA");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI:
      return std::string("RPI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI:
      return std::string("LPI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI:
      return std::string("RAI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI:
      return std::string("LAI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS:
      return std::string("RPS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS:
      return std::string("LPS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS:
      return std::string("RAS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS:
      return std::string("LAS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI:
      return std::string("PRI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI:
      return std::string("PLI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI:
      return std::string("ARI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI:
      return std::string("ALI");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS:
      return std::string("PRS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS:
      return std::string("PLS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS:
      return std::string("ARS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS:
      return std::string("ALS");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR:
      return std::string("IPR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR:
      return std::string("SPR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR:
      return std::string("IAR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR:
      return std::string("SAR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL:
      return std::string("IPL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL:
      return std::string("SPL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL:
      return std::string("IAL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL:
      return std::string("SAL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR:
      return std::string("PIR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR:
      return std::string("PSR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR:
      return std::string("AIR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR:
      return std::string("ASR");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL:
      return std::string("PIL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL:
      return std::string("PSL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL:
      return std::string("AIL");
    case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL:
      return "ASL";
    default:
      {
      std::stringstream x;
      x << (in & 0xff) << ", " << ((in >> 8) & 0xff) << ", " << ((in >> 16) && 0xff);
      return x.str();
      }
    }
}

int main( int argc, char *argv[] )
{

  const char *inputFile = NULL;
  unsigned int dim;

  TCLAP::CmdLine cmd("Prints some image information.", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Input image",true,"","string",cmd);
  TCLAP::ValueArg<unsigned int> dimArg("d","dimension","Image dimension (3 / 4)",true,3,"unsigned int",cmd);

  cmd.parse( argc, argv );

  inputFile = inputArg.getValue().c_str();
  dim       = dimArg.getValue();

  typedef itk::Image< short, 3 >  ImageType;

  if (dim==4)
  {
    typedef itk::Image< short, 4 >  SequenceType;

    typedef itk::ImageFileReader< SequenceType > SequenceReaderType;
    SequenceReaderType::Pointer  sequenceReader  = SequenceReaderType::New();
    sequenceReader->SetFileName(  inputFile );
    sequenceReader->Update();

    SequenceType::Pointer sequence = sequenceReader->GetOutput();
    sequence->Print(std::cout);

    typedef itk::ExtractImageFilter< SequenceType, ImageType > ExtractorType;
    ExtractorType::Pointer extractor  = ExtractorType::New();
    extractor -> SetInput( sequence );

    SequenceType::RegionType region = sequence -> GetLargestPossibleRegion();
    SequenceType::IndexType  start  = region.GetIndex();
    SequenceType::SizeType   size   = sequence -> GetLargestPossibleRegion().GetSize();

    size[3] = 0;
    region.SetSize( size );
    extractor -> SetExtractionRegion( region );
    extractor -> Update();

    std::cout << std::endl << "Anatomical orientation = " << SO_OrientationToString(
        itk::SpatialOrientationAdapter().FromDirectionCosines( extractor -> GetOutput() -> GetDirection() )
        ) << std::endl;

  } else if (dim==3)
    {
      typedef itk::ImageFileReader< ImageType > ImageReaderType;
      ImageReaderType::Pointer  imageReader  = ImageReaderType::New();
      imageReader->SetFileName(  inputFile );
      imageReader->Update();
      ImageType::Pointer image = imageReader->GetOutput();
      image->Print(std::cout);

      std::cout << std::endl << "Anatomical orientation = " << SO_OrientationToString(
        itk::SpatialOrientationAdapter().FromDirectionCosines( image -> GetDirection() )
        ) << std::endl;
    } else
      {
        std::cout << "ERROR: Image dimension unsupported." << std::endl;
        return EXIT_FAILURE;

      }

  return EXIT_SUCCESS;
}

