/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 29/06/2012
  Author(s): Julien Pontabry (pontabry at unistra dot fr)
  
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
#include "cstdlib"
#include "iostream"
#include "string"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataWriter.h"

// Local includes
#include "btkMacro.h"
#include "btkFileHelper.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkImageHelper.h"
#include "btkDiffusionTensorReconstructionFilter.h"
#include "btkSphericalHarmonicsDiffusionDecompositionFilter.h"
#include "btkDiffusionModel.h"
#include "btkTensorModel.h"
#include "btkOrientationDiffusionFunctionModel.h"
#include "btkTractographyAlgorithm.h"
#include "btkStreamlineTractographyAlgorithm.h"
#include "btkCommandProgressUpdate.h"
#include "btkPolyDataColorLinesByOrientation.h"


int main(int argc, char *argv[])
{
    try
    {
        //
        // Program's command line
        //

        // Defines the command line parser
        TCLAP::CmdLine cmd("BTK Tractography program", ' ', "2.0");

        // Defines arguments
        TCLAP::ValueArg< std::string > dwiSequenceFileNameArg("d", "dwi_sequence", "DWI Sequence", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >        maskFileNameArg("m", "mask", "Propagation ROI", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >         roiFileNameArg("r", "roi_image", "ROI image", true, "", "string", cmd);
        TCLAP::MultiArg< unsigned short >           labelsArg("l", "label", "Label selection in ROI image", false, "int", cmd);

        TCLAP::SwitchArg     streamlineAlgorithmArg("", "streamline", "Use streamline algorithm for tractography", cmd, false);
        TCLAP::SwitchArg          tensorModelingArg("", "tensor", "Use tensor modeling for diffusion", cmd, false);

        TCLAP::ValueArg< std::string >        modelFileNameArg("", "model", "Model image", false, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNamePrefixArg("o", "output", "Prefix of the filenames of the outputs", false, "tractography", "string", cmd);

        TCLAP::ValueArg< float >       stepSizeArg("", "step_size", "Step size between two points of the solution", false, 0.5, "positive real", cmd);
        TCLAP::ValueArg< float >    seedSpacingArg("", "seed_spacing", "Spacing between two seeds (in mm)", false, 1, "positive real", cmd);
        TCLAP::SwitchArg                 useRK4Arg("", "RK4", "Use the Runge-Kutta method at order 4 instead of Euler's methods for path computing in streamline algorithm", cmd, false);
        TCLAP::ValueArg< float > thresholdAngleArg("", "threshold_angle", "Threshold angle (in degrees) between two successive displacements", false, 75, "positive real", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string      dwiSequenceFileName = dwiSequenceFileNameArg.getValue();
        std::string             maskFileName = maskFileNameArg.getValue();
        std::string              roiFileName = roiFileNameArg.getValue();
        std::vector< unsigned short > labels = labelsArg.getValue();

        bool     streamlineAlgorithm = streamlineAlgorithmArg.getValue();
        bool          tensorModeling = tensorModelingArg.getValue();

        std::string        modelFileName = modelFileNameArg.getValue();
        std::string outputFileNamePrefix = outputFileNamePrefixArg.getValue();

        float       stepSize = stepSizeArg.getValue();
        bool          useRK4 = useRK4Arg.getValue();
        float thresholdAngle = 0.017453292519943 * thresholdAngleArg.getValue();
        float    seedSpacing = seedSpacingArg.getValue();


        //
        // Read input files
        //

        btk::DiffusionSequence::Pointer         dwiSequence = btk::DiffusionSequenceHelper::ReadSequence(dwiSequenceFileName);
        btk::TractographyAlgorithm::MaskImage::Pointer mask = btk::ImageHelper< btk::TractographyAlgorithm::MaskImage >::ReadImage(maskFileName);
        btk::TractographyAlgorithm::LabelImage::Pointer roi = btk::ImageHelper< btk::TractographyAlgorithm::LabelImage >::ReadImage(roiFileName);


        //
        // Estimate diffusion modeling
        //

        btk::DiffusionModel::Pointer model = NULL;

        if(tensorModeling)
        {
            btk::DiffusionTensorReconstructionFilter::OutputImageType::Pointer tensorImage = NULL;

            if(modelFileName.empty() || !btk::FileHelper::FileExist(modelFileName))
            {
                btkCoutMacro("Estimating tensor modeling...");

                btk::DiffusionTensorReconstructionFilter::Pointer reconstructionFilter = btk::DiffusionTensorReconstructionFilter::New();
                reconstructionFilter->SetInput(dwiSequence);
                reconstructionFilter->SetThreshold(0);
                reconstructionFilter->Update();

                tensorImage = reconstructionFilter->GetOutput();

                if(!modelFileName.empty())
                {
                    btk::ImageHelper< btk::DiffusionTensorReconstructionFilter::OutputImageType >::WriteImage(tensorImage, modelFileName);
                }
            }
            else // There is an existant model image and we load it.
            {
                btkCoutMacro("Loading tensor image...");

                tensorImage = btk::ImageHelper< btk::DiffusionTensorReconstructionFilter::OutputImageType >::ReadImage(modelFileName);
            }

            btk::TensorModel::Pointer tensorModel = btk::TensorModel::New();
            tensorModel->SetInputModelImage(tensorImage);
            tensorModel->SetBValue(dwiSequence->GetBValues()[1]);
            tensorModel->SetSphericalResolution(300);
            tensorModel->Update();

            model = tensorModel;

            btkCoutMacro("done.");
        }
        else // ODF modeling
        {
            btk::SphericalHarmonicsDiffusionDecompositionFilter::OutputImageType::Pointer shCoefficientsImage = NULL;

            if(modelFileName.empty() || !btk::FileHelper::FileExist(modelFileName))
            {
                btkCoutMacro("Estimating spherical harmonics Q-Ball modeling...");

                btk::SphericalHarmonicsDiffusionDecompositionFilter::Pointer decompositionFilter = btk::SphericalHarmonicsDiffusionDecompositionFilter::New();
                decompositionFilter->SetInput(dwiSequence);
                decompositionFilter->Update();

                shCoefficientsImage = decompositionFilter->GetOutput();

                if(!modelFileName.empty())
                {
                    btk::ImageHelper< btk::SphericalHarmonicsDiffusionDecompositionFilter::OutputImageType >::WriteImage(shCoefficientsImage, modelFileName);
                }
            }
            else // There is an existant model image and we load it.
            {
                btkCoutMacro("Loading spherical harmonics coefficients image...");

                shCoefficientsImage = btk::ImageHelper< btk::SphericalHarmonicsDiffusionDecompositionFilter::OutputImageType >::ReadImage(modelFileName);
            }

            btk::OrientationDiffusionFunctionModel::Pointer odfModel = btk::OrientationDiffusionFunctionModel::New();
            odfModel->SetInputModelImage(shCoefficientsImage);
            odfModel->SetBValue(dwiSequence->GetBValues()[1]);
            odfModel->SetSphericalResolution(300);
            odfModel->Update();

            model = odfModel;

            btkCoutMacro("done.");
        }


        //
        // Tractography
        //

        btk::TractographyAlgorithm::Pointer algorithm = NULL;

//        if(streamlineAlgorithm)
        {
            btkCoutMacro("Setting up streamline algorithm...");

            btk::StreamlineTractographyAlgorithm::Pointer strAlgorithm = btk::StreamlineTractographyAlgorithm::New();
            strAlgorithm->SetStepSize(stepSize);
            strAlgorithm->UseRungeKuttaOrder4(useRK4);
            strAlgorithm->SetThresholdAngle(thresholdAngle);
            strAlgorithm->SetSeedSpacing(seedSpacing);

            algorithm = strAlgorithm;

            btkCoutMacro("done.");
        }
//        else // Particle filtering algorithm
//        {
//            btkCoutMacro("Setting up particle filtering algorithm...");

//            btkCoutMacro("done.");
//        }

        btkCoutMacro("Processing...");

        // Add an observer for progress bar
        btk::CommandProgressUpdate::Pointer observer = btk::CommandProgressUpdate::New();
        algorithm->AddObserver(itk::ProgressEvent(), observer);

        // Process
        algorithm->SetDiffusionSequence(dwiSequence);
        algorithm->SetDiffusionModel(model);
        algorithm->SetMask(mask);
        algorithm->SetRegionsOfInterest(roi);
        algorithm->SetSeedLabels(labels);
        algorithm->Update();

        btkCoutMacro("done.");


        //
        // Write output
        //

        btkCoutMacro("Writing outputs...");

        // Write each fibers corresponding to labels
        std::vector< vtkSmartPointer< vtkAppendPolyData > > fibers = algorithm->GetOutputFiber();

        for(unsigned int i = 0; i < fibers.size(); i++)
        {
            if(fibers[i] != NULL && fibers[i]->GetNumberOfInputPorts() > 0)
            {
                std::stringstream filename;
                filename << outputFileNamePrefix << "-" << i+1 << ".vtk";

                fibers[i]->Update();

                // Color fibers by mean orientation
                vtkSmartPointer< btk::PolyDataColorLinesByOrientation > colorFilter = vtkSmartPointer< btk::PolyDataColorLinesByOrientation >::New();
                colorFilter->SetInput(fibers[i]->GetOutput());
                colorFilter->Update();

                vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
                writer->SetInput(colorFilter->GetOutput());
                writer->SetFileName(filename.str().c_str());
                writer->Write();
            }
        }

        btkCoutMacro("done.");
    }
    catch(TCLAP::ArgException &exception)
    {
        btkCoutMacro("Command line error: " << exception.error() << " for argument " << exception.argId());
        exit(EXIT_FAILURE);
    }
    catch(std::string &message)
    {
        btkCoutMacro("Exception raised: " << message << std::endl);
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
