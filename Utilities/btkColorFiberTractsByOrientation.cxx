/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/01/2013
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
 * @file btkColorFiberTractsByOrientation.cxx
 * @author Julien Pontabry
 * @date 04/01/2013
 * @brief Color the fibers (vtk polydata) by their orientation (local or global).
 * @ingroup Diffusion
 */

// TCLAP includes
#include "tclap/CmdLine.h"

// STL includes
#include "cstdlib"
#include "string"
#include "iostream"

// ITK includes
#include "itkMacro.h"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkCallbackCommand.h"

// Local includes
#include "btkPolyDataColorLinesByOrientation.h"


/**
 * @brief Progress bar display.
 * @param caller Object which called the progress event.
 * @param eventId Id of the event.
 * @param clientData Client data.
 * @param callData Call data.
 */
void progress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
    btk::PolyDataColorLinesByOrientation *filter = static_cast< btk::PolyDataColorLinesByOrientation * >(caller);

    std::cout << "\tProgress: " << static_cast< unsigned int >(filter->GetProgress()*100) << "%\r" << std::flush;
}

/**
 * @brief Main function of the program.
 */
int main(int argc, char *argv[])
{
    try
    {
        //
        // Command line parser
        //

        // Command line definition
        TCLAP::CmdLine cmd("Color the fibers (polydata) by their mean orientation", ' ', "1.0", true);

        // Arguments
        TCLAP::ValueArg< std::string >  inputFileNameArg("i", "input", "Input polydata filename (VTK)", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output polydata filename (VTK)", true, "", "string", cmd);

        TCLAP::SwitchArg colorByLocalOrientationArg("", "local_orientation_color", "Color the output fibers by local orientation instead of mean orientation", cmd);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get arguments' values
        std::string  inputFileName = inputFileNameArg.getValue();
        std::string outputFileName = outputFileNameArg.getValue();

        bool colorByLocalOrientation = colorByLocalOrientationArg.getValue();


        //
        // Processing
        //

        // Load fiber bundle
        std::cout << "Reading file... " << std::flush;

        vtkSmartPointer< vtkPolyDataReader > reader = vtkSmartPointer< vtkPolyDataReader >::New();
        reader->SetFileName(inputFileName.c_str());
        reader->Update();

        std::cout << "done." << std::endl;

        // Process
        std::cout << "Processing... " << std::endl;

        vtkSmartPointer< vtkCallbackCommand > observer = vtkSmartPointer< vtkCallbackCommand >::New();
        observer->SetCallback(progress);

        vtkSmartPointer< btk::PolyDataColorLinesByOrientation > filter = vtkSmartPointer< btk::PolyDataColorLinesByOrientation >::New();
        filter->SetInput(reader->GetOutput());

        if(colorByLocalOrientation)
        {
            filter->SetColorOrientation(btk::PolyDataColorLinesByOrientation::COLOR_LOCAL_ORIENTATION);
        }

        filter->AddObserver(vtkCommand::ProgressEvent, observer);
        filter->Update();

        std::cout << std::endl << "done." << std::endl;

        // Write filtered fiber bundle
        std::cout << "Writing file... " << std::flush;

        vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
        writer->SetFileName(outputFileName.c_str());
        writer->SetInput(filter->GetOutput());
        writer->Update();

        std::cout << "done." << std::endl;
    }
    catch(TCLAP::ArgException &e)
    {
        std::cerr << "Error: " << std::endl << e.error() << " for arg " << e.argId() << std::endl;
        return EXIT_FAILURE;
    }
    catch(std::string &e)
    {
        std::cerr << "Error: " << std::endl << e << std::endl;
        return EXIT_FAILURE;
    }
    catch(itk::ExceptionObject &e)
    {
        std::cerr << "Error: " << std::endl << e.GetDescription() << std::endl;
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;
}
