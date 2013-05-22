/*
 * Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 *
 * 31/01/2012
 * < pontabry at unistra dot fr >
 *
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-B license and that you accept its terms.
 */


 // TCLAP : Templatized C++ Command Line Parser includes
#include <tclap/CmdLine.h>

// STL includes
#include "cstdlib"
#include "string"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkPoint.h"
#include "itkContinuousIndex.h"


const unsigned int Dimension = 3;
typedef float PixelType;

typedef itk::Image<PixelType,Dimension> Image;
typedef itk::ImageFileReader<Image>     ImageFileReader;
typedef itk::ImageFileWriter<Image>     ImageFileWriter;
typedef itk::ImageRegionIterator<Image> ImageIterator;


int main(int argc, char *argv[])
{
    try
    {
        //
        // Parse program's arguments
        //

        // Define command line variables
        std::vector<std::string> bundleFileName;
        std::string referenceFileName;
        std::string outputFileName;

        // Define command line parser
        TCLAP::CmdLine cmd("Convert a fibers bundle to its associated connectivity map", ' ', "0.1");

        // Define command line arguments
        TCLAP::MultiArg<std::string> bundleFileNameArg("b", "bundle", "Fibers bundle filename (possible multiple inputs)", true, "string", cmd);
        TCLAP::ValueArg<std::string> referenceFileNameArg("r", "reference", "Reference image (space)", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputFileNameArg("o", "output", "Output map filename", false, "map.nii.gz", "string", cmd);

        // Parse arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        bundleFileName    = bundleFileNameArg.getValue();
        referenceFileName = referenceFileNameArg.getValue();
        outputFileName    = outputFileNameArg.getValue();


        //
        // Load files
        //

        std::cout << "Loading reference image..." << std::flush;

        // Load reference
        ImageFileReader::Pointer referenceReader = ImageFileReader::New();
        referenceReader->SetFileName(referenceFileName.c_str());
        referenceReader->Update();
        Image::Pointer reference = referenceReader->GetOutput();

        std::cout << "done." << std::endl;


        //
        // Create output image
        //

        std::cout << "Creating empty output image..." << std::flush;

        Image::Pointer output = Image::New();
        output->SetRegions(reference->GetLargestPossibleRegion());
        output->SetOrigin(reference->GetOrigin());
        output->SetSpacing(reference->GetSpacing());
        output->SetDirection(reference->GetDirection());
        output->Allocate();
        output->FillBuffer(0);

        std::cout << "done." << std::endl;


        //
        // Conversion
        //

        std::cout << "Converting fibers bundle to connectivity map..." << std::flush;


        for(int bundleIndex=0; bundleIndex<bundleFileName.size(); bundleIndex++)
        {
          
        // Load bundle
        vtkSmartPointer<vtkPolyDataReader> bundleReader = vtkSmartPointer<vtkPolyDataReader>::New();
        bundleReader->SetFileName(bundleFileName[bundleIndex].c_str());
        bundleReader->Update();

        vtkSmartPointer<vtkPolyData> bundle = bundleReader->GetOutput();
        vtkSmartPointer<vtkCellArray> lines = bundle->GetLines();

        vtkIdType numberOfPoints, *pointIds;

        // For each fibers
        while(lines->GetNextCell(numberOfPoints, pointIds) != 0)
        {
            // For each point of the current fiber
            for(unsigned int i = 0; i < numberOfPoints; i++)
            {
                // Get current point's coordinates
                double worldCoordinates[3];
                bundle->GetPoint(pointIds[i], worldCoordinates);

                // Convert world coordinates to index coordinates
                itk::Point<PixelType,Dimension> worldPoint;
                worldPoint[0] = -worldCoordinates[0]; worldPoint[1] = -worldCoordinates[1]; worldPoint[2] = worldCoordinates[2];

                itk::ContinuousIndex<PixelType,Dimension> continuousIndexCoordinates;
                output->TransformPhysicalPointToContinuousIndex(worldPoint, continuousIndexCoordinates);

                //
                // Tri-linear interpolation
                //

                // Continuous index
                float cx = continuousIndexCoordinates[0];
                float cy = continuousIndexCoordinates[1];
                float cz = continuousIndexCoordinates[2];

                // Discrete index
                short ix = (short)std::floor(cx);
                short iy = (short)std::floor(cy);
                short iz = (short)std::floor(cz);

                // Weights
                float wx = (float)(cx - ix);
                float wy = (float)(cy - iy);
                float wz = (float)(cz - iz);

                // Add computed value to each pixels around the center point
                Image::IndexType index;

                index[0] = ix; index[1] = iy; index[2] = iz;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + (1-wz)*(1-wx)*(1-wy));

                index[0] = ix; index[1] = iy+1; index[2] = iz;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + (1-wz)*(1-wx)*wy);

                index[0] = ix+1; index[1] = iy; index[2] = iz;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + (1-wz)*wx*(1-wy));

                index[0] = ix+1; index[1] = iy+1; index[2] = iz;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + (1-wz)*wx*wy);

                index[0] = ix; index[1] = iy; index[2] = iz+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + wz*(1-wx)*(1-wy));

                index[0] = ix; index[1] = iy+1; index[2] = iz+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + wz*(1-wx)*wy);

                index[0] = ix+1; index[1] = iy; index[2] = iz+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + wz*wx*(1-wy));

                index[0] = ix+1; index[1] = iy+1; index[2] = iz+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < output->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < output->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < output->GetLargestPossibleRegion().GetSize(2))
                    output->SetPixel(index, output->GetPixel(index) + wz*wx*wy);
            }

        }
          
        
        }//End of loop of bundles

        ImageIterator it(output, output->GetLargestPossibleRegion());

        float maxIntensity = 0;

        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            if(it.Get() > maxIntensity)
                maxIntensity = it.Get();
        }

        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
            it.Set(it.Get()/maxIntensity);

        std::cout << "done." << std::endl;


        //
        // Write output
        //

        std::cout << "Writing output connectivity map..." << std::flush;

        ImageFileWriter::Pointer outputWriter = ImageFileWriter::New();
        outputWriter->SetInput(output);
        outputWriter->SetFileName(outputFileName.c_str());
        outputWriter->Update();

        std::cout << "done." << std::endl;
    }
    catch(TCLAP::ArgException &e)
    {
        std::cerr << "Error (arguments): " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch(itk::ImageFileReaderException &e)
    {
        std::cerr << "Error (input images): an image filename is wrong or does not exist ; program exit !" << std::endl;
        std::cerr << "ITK error trace:" << e;
        exit(EXIT_FAILURE);
    }
    catch(std::string &e)
    {
        std::cerr << std::endl << e << " ; program exit !" << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
