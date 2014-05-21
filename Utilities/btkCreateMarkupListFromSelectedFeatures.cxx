/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 05/12/2013
  Author(s): Julien Pontabry (pontabry@unistra.fr)

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


/**
 * @file btkCreateMarkupListFromSelectedFeatures.cxx
 * @author Julien Pontabry
 * @date 18/03/2013
 * @ingroup FeatureSelection
 * @brief Create a markup list from selected features that can be opened within Slicer.
 */

// TCLAP includes
#include <tclap/CmdLine.h>

// STL includes
#include "fstream"

// ITK includes
#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegionIterator.h"

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"

// Definitions
typedef itk::Image< unsigned char,3 >         ImageMask;
typedef itk::ImageMaskSpatialObject< 3 >      Mask;
typedef itk::ImageRegionIterator< ImageMask > ImageMaskIterator;


/**
 * @brief Main function of the program.
 */
int main(int argc, char *argv[])
{
    try
    {
        //
        // Program's command line
        //

        // Defines the command line parser
        TCLAP::CmdLine cmd("BTK markup list creation from selected features", ' ', "1.0");

        // Defines arguments
        TCLAP::ValueArg< std::string > selectedFeaturesFileNamesArg("r", "reduced_mask", "Selected features filename", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output markup list filename (Default: \"SelectedFeatures.fcsv\")", false, "SelectedFeatures.fcsv", "string", cmd);

        // Options
        TCLAP::SwitchArg vtkFileArg("", "vtk", "Produce VTK files in output", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string selectedFeaturesFileName = selectedFeaturesFileNamesArg.getValue();
        std::string           outputFileName = outputFileNameArg.getValue();

        bool vtkFile = vtkFileArg.getValue();


        //
        // Read input file
        //

        ImageMask::Pointer selectedFeaturesImage = btk::ImageHelper< ImageMask >::ReadImage(selectedFeaturesFileName);


        //
        // Processing
        //

        // Compute the object mask
        Mask::Pointer mask = Mask::New();
        mask->SetImage(selectedFeaturesImage);
        mask->Update();

        // Get the masked region
        Mask::RegionType maskedRegion = mask->GetAxisAlignedBoundingBoxRegion();

        // Open output file
        std::ofstream outputFileStream(outputFileName.c_str());

        if(!outputFileStream.is_open())
        {
            outputFileStream.close();
            btkException("Error when writing in file !");
        }

        if(vtkFile)
        {
            outputFileStream << "# vtk DataFile Version 3.0" << std::endl;
            outputFileStream << "vtk output" << std::endl;
            outputFileStream << "ASCII" << std::endl;
            outputFileStream << "DATASET POLYDATA" << std::endl;

            // Run through selected features to add to markup list
            ImageMaskIterator selectedFeaturesIt(selectedFeaturesImage, maskedRegion);

            // Get the number of features
            unsigned int featureCount = 0;

            for(selectedFeaturesIt.GoToBegin(); !selectedFeaturesIt.IsAtEnd(); ++selectedFeaturesIt)
            {
                if(selectedFeaturesIt.Get() > 0)
                {
                    featureCount++;
                }
            }

            outputFileStream << "POINTS " << featureCount << " float" << std::endl;

            // Write features
            for(selectedFeaturesIt.GoToBegin(); !selectedFeaturesIt.IsAtEnd(); ++selectedFeaturesIt)
            {
                if(selectedFeaturesIt.Get() > 0)
                {
                    // Convert image coordinate to physical coordinates
                    ImageMask::IndexType index = selectedFeaturesIt.GetIndex();
                    ImageMask::PointType point;
                    selectedFeaturesImage->TransformIndexToPhysicalPoint(index,point);

                    outputFileStream << -point[0] << " " << -point[1] << " " << point[2] << std::endl;
                }
            }
        }
        else // fcsv file
        {
            outputFileStream << "# Markups fiducial file version = 4.3" << std::endl;
            outputFileStream << "# CoordinateSystem = 0" << std::endl;
            outputFileStream << "# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID" << std::endl;

            // Get origin of the image
            ImageMask::PointType origin = selectedFeaturesImage->GetOrigin();

            // Run through selected features to add to markup list
            ImageMaskIterator selectedFeaturesIt(selectedFeaturesImage, maskedRegion);
            unsigned int featureCount = 1;

            for(selectedFeaturesIt.GoToBegin(); !selectedFeaturesIt.IsAtEnd(); ++selectedFeaturesIt)
            {
                if(selectedFeaturesIt.Get() > 0)
                {
                    // Convert image coordinate to physical coordinates
                    ImageMask::IndexType index = selectedFeaturesIt.GetIndex();
                    ImageMask::PointType point;
                    selectedFeaturesImage->TransformIndexToPhysicalPoint(index,point);

                    outputFileStream << "vtkMRMLMarkupsFiducialNode_0," << -point[0] << "," << -point[1] << "," << point[2] << "," << origin[0] << "," << origin[1] << "," << origin[2] << ",1,1,1,1,F-" << featureCount << ",Feature " << featureCount << "," << std::endl;

                    featureCount++;
                }
            }
        }

        // Close output stream
        outputFileStream.close();
    }
    catch(TCLAP::ArgException &exception)
    {
        btkCoutMacro("Command line error: " << exception.error() << " for argument " << exception.argId());
        exit(EXIT_FAILURE);
    }
    catch(itk::ExceptionObject &e)
    {
        btkCoutMacro("Error (" << e.GetLocation() << "): " << e.GetDescription() << std::endl);
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
