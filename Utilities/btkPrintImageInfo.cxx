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


// TCLAP includes
#include <tclap/CmdLine.h>

// STL includes
#include "iostream"
#include "string"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkSpatialOrientationAdapter.h"
#include "vnl/vnl_matrix.h"


// Spatial orientation (RAS, LPS, etc.)
typedef itk::SpatialOrientation::ValidCoordinateOrientationFlags SpatialOrientation;

/**
  * @brief Get the spatial orientation token (RAS, LPS, etc.) from the ITK spatial orientation.
  * @param orientation ITK spatial orientation
  * @return The corresponding spatial orientation token (RAS, LPS, etc.)
  */
std::string SpatialOrientationToString(SpatialOrientation orientation)
{
    std::string orientationString;

    switch(orientation)
    {
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIP:
            orientationString = "RIP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIP:
            orientationString = "LIP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSP:
            orientationString = "RSP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSP:
            orientationString = "LSP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RIA:
            orientationString = "RIA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LIA:
            orientationString = "LIA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RSA:
            orientationString = "RSA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LSA:
            orientationString = "LSA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRP:
            orientationString = "IRP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILP:
            orientationString = "ILP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRP:
            orientationString = "SRP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLP:
            orientationString = "SLP";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IRA:
            orientationString = "IRA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ILA:
            orientationString = "ILA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SRA:
            orientationString = "SRA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SLA:
            orientationString = "SLA";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPI:
            orientationString = "RPI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPI:
            orientationString = "LPI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAI:
            orientationString = "RAI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAI:
            orientationString = "LAI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RPS:
            orientationString = "RPS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LPS:
            orientationString = "LPS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_RAS:
            orientationString = "RAS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_LAS:
            orientationString = "LAS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRI:
            orientationString = "PRI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLI:
            orientationString = "PLI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARI:
            orientationString = "ARI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALI:
            orientationString = "ALI";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PRS:
            orientationString = "PRS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PLS:
            orientationString = "PLS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ARS:
            orientationString = "ARS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ALS:
            orientationString = "ALS";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPR:
            orientationString = "IPR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPR:
            orientationString = "SPR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAR:
            orientationString = "IAR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAR:
            orientationString = "SAR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IPL:
            orientationString = "IPL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SPL:
            orientationString = "SPL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_IAL:
            orientationString = "IAL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_SAL:
            orientationString = "SAL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIR:
            orientationString = "PIR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSR:
            orientationString = "PSR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR:
            orientationString = "AIR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASR:
            orientationString = "ASR";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PIL:
            orientationString = "PIL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_PSL:
            orientationString = "PSL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL:
            orientationString = "AIL";
            break;
        case itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_ASL:
            orientationString = "ASL";
            break;
        default:
            std::stringstream x;
            x << (orientation & 0xff) << ", " << ((orientation >> 8) & 0xff) << ", " << ((orientation >> 16) && 0xff);
            orientationString = x.str();
    }

    return orientationString;
}


// Input sequence (4D)
typedef itk::Image<short,4>            Sequence;
typedef itk::ImageFileReader<Sequence> SequenceReader;

// Input volume (3D)
typedef itk::Image<short,3>         Image;
typedef itk::ImageFileReader<Image> ImageReader;


int main( int argc, char *argv[] )
{
    //
    // Parse program's arguments
    //

    // Define command line parameters
    std::string inputFileName;
    unsigned int dimension;

    // Define command line parser
    TCLAP::CmdLine cmd("Prints some image information.", ' ', "Unversioned");

    // Define command line arguments
    TCLAP::ValueArg<std::string> inputFileNameArg("i", "input", "Input image", true, "", "string", cmd);
    TCLAP::ValueArg<unsigned int> dimensionArg("d", "dimension", "Image dimension (default: 3; values: 3-4)", true, 3, "unsigned int", cmd);

    // Parse arguments
    cmd.parse( argc, argv );

    // Get back arguments' values
    inputFileName = inputFileNameArg.getValue();
    dimension     = dimensionArg.getValue();


    //
    // Processing
    //

    try
    {
        if(dimension == 4)
        {
            // Read sequence
            SequenceReader::Pointer  sequenceReader  = SequenceReader::New();
            sequenceReader->SetFileName(inputFileName);
            sequenceReader->Update();

            Sequence::Pointer sequence = sequenceReader->GetOutput();

            // Display informations
            std::cout << "Dimension: " << Sequence::GetImageDimension() << std::endl << std::endl;

            std::cout << "Size:    " << sequence->GetLargestPossibleRegion().GetSize() << std::endl;
            std::cout << "Origin:  " << sequence->GetOrigin() << std::endl;
            std::cout << "Spacing: " << sequence->GetSpacing() << std::endl << std::endl;

            std::cout << "Direction:" << std::endl << sequence->GetDirection() << std::endl;

            // Display spatial orientation
            Image::DirectionType direction(sequence->GetDirection().GetVnlMatrix().extract(3,3));
            std::cout << "Anatomical orientation: " << SpatialOrientationToString(itk::SpatialOrientationAdapter().FromDirectionCosines(direction)) << std::endl;
        }
        else if (dimension == 3)
        {
            // Read image
            ImageReader::Pointer imageReader = ImageReader::New();
            imageReader->SetFileName(inputFileName);
            imageReader->Update();

            Image::Pointer image = imageReader->GetOutput();

            // Display informations
            std::cout << "Dimension: " << Image::GetImageDimension() << std::endl << std::endl;

            std::cout << "Size:    " << image->GetLargestPossibleRegion().GetSize() << std::endl;
            std::cout << "Origin:  " << image->GetOrigin() << std::endl;
            std::cout << "Spacing: " << image->GetSpacing() << std::endl << std::endl;

            std::cout << "Direction:" << std::endl << image->GetDirection() << std::endl;

            // Display spatial orientation
            std::cout << std::endl << "Anatomical orientation: " << SpatialOrientationToString(itk::SpatialOrientationAdapter().FromDirectionCosines(image->GetDirection())) << std::endl;
        }
        else // dimension != 3 && dimension != 4
        {
            std::stringstream message;
            message << "Image dimension (" << dimension << ") is unsupported ! You should use either dimension 3 or 4.";
            throw std::string(message.str());
        }
    }
    catch(itk::ExceptionObject &error)
    {
        std::cout << "ITK error: " << error << std::endl;
    }
    catch(std::string &message)
    {
        std::cout << "Error: " << message << std::endl;
    }

    return EXIT_SUCCESS;
}

