/*
 * Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 *
 * 24 januar 2011
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
#include "string"
#include "iostream"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataWriter.h"

// Local includes
#include "btkTypes.h"
#include "btkPoint.h"
#include "btkDirection.h"
#include "btkDTFPSignal.h"
#include "btkDTFPSignalExtractor.h"
#include "btkDTFPImportanceDensity.h"
#include "btkDTFPAPrioriDensity.h"
#include "btkDTFPLikelihoodDensity.h"
#include "btkDTFPParticleFilter.h"


using namespace btk;


int main(int argc, char *argv[])
{
    // Command line variables
    std::string dwiFileName;
    std::string vecFileName;
    unsigned int valFileName;
    std::string maskFileName;
    std::string labelFilename;

    std::string outMapFileName;
    std::string outFibersFileName;

    bool verboseMode;
    bool quietMode;
    bool saveTmpFiles;
    bool lps;

    unsigned int nbOfParticles;
    Real stepSize;
    Real epsilon;
    Real Kappa;
    Float seedSpacing;


    try
    {
        //
        // Program's command line parser definition
        //

            // Defines command line parser
            TCLAP::CmdLine cmd("BTK Tractography", ' ', "0.2");

            // Defines arguments
            TCLAP::ValueArg<std::string>   dwiArg("d", "dwi", "Dwi sequence", true, "", "string", cmd);
            TCLAP::ValueArg<std::string>   vecArg("g", "gradient_vectors", "Gradient vectors", true, "", "string", cmd);
            TCLAP::ValueArg<unsigned int>   valArg("b", "b_values", "B-values", true, 1500, "unsigned int", cmd);
            TCLAP::ValueArg<std::string>  maskArg("m", "mask", "White matter mask", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> labelArg("l", "label", "Label volume of seeds", true, "", "string", cmd);

            TCLAP::ValueArg<std::string>    outMapArg("", "map", "Output connection map file", false, "map", "string", cmd);
            TCLAP::ValueArg<std::string> outFibersArg("", "fibers", "Output fibers file", false, "fibers", "string", cmd);
            TCLAP::SwitchArg lpsSwitchArg("", "lps", "Word coordinates expressed in LPS (Left-Posterior-Superior). By default RAS (Right-Anterior-Superior) is used.", cmd, false);

            TCLAP::SwitchArg verboseSwitchArg("", "verbose", "Display more informations on standard output", cmd, false);
            TCLAP::SwitchArg quietSwitchArg("", "quiet", "Display no information on either standard and error outputs", cmd, false);
            TCLAP::SwitchArg saveTmpSwitchArg("", "save_temporary_files", "Save diffusion signal, model coefficients, variance and spherical coordinates of gradient directions in files", cmd, false);


            TCLAP::ValueArg<unsigned int> particlesArg("", "number_of_particles", "Number of particles", false, 1000, "unsigned int", cmd);
            TCLAP::ValueArg<Real>    epsilonArg("", "resampling_threshold", "Resampling treshold", false, 0.01, "Real", cmd);
            TCLAP::ValueArg<Real>    stepSizeArg("", "step_size", "Step size of particles displacement", false, 0.5, "Real", cmd);
            TCLAP::ValueArg<Real>    KappaArg("", "curve_constraint", "Curve constraint of a particle's trajectory", false, 30.0, "Real", cmd);
            TCLAP::ValueArg<Real>    seedSpacingArg("", "seed_spacing", "Spacing in mm between seeds", false, 1.0, "float", cmd);

            // Parsing arguments
            cmd.parse(argc, argv);

            // Get back arguments' values
            dwiFileName    = dwiArg.getValue();
            vecFileName    = vecArg.getValue();
            valFileName    = valArg.getValue();
            maskFileName   = maskArg.getValue();
            labelFilename  = labelArg.getValue();

            outMapFileName    = outMapArg.getValue();
            outFibersFileName = outFibersArg.getValue();

            verboseMode  = verboseSwitchArg.getValue();
            quietMode    = quietSwitchArg.getValue();
            saveTmpFiles = saveTmpSwitchArg.getValue();
            lps          = lpsSwitchArg.getValue();

            nbOfParticles  = particlesArg.getValue();
            epsilon        = epsilonArg.getValue();
            stepSize       = stepSizeArg.getValue();
            Kappa          = KappaArg.getValue();
            seedSpacing    = seedSpacingArg.getValue();
    }
    catch(TCLAP::ArgException &e)
    {
        std::cout << "TCLAP error: " << e.error() << " for argument " << e.argId() << std::endl;
        std::exit(EXIT_FAILURE);
    }


    char displayMode = quietMode ? 0 : (verboseMode ? 2 : 1);


    //
    // Diffusion signal extraction
    //

        DTFPSignalExtractor *extractor = new DTFPSignalExtractor(vecFileName, dwiFileName, maskFileName, displayMode);
        extractor->extract();

        std::vector<Vector> *directions  = extractor->GetDirections();
        Sequence::Pointer       signal   = extractor->GetSignal();
        std::vector<Real>      *sigmas   = extractor->GetSigmas();
        Mask::Pointer           mask     = extractor->GetMask();
        Image::Pointer          baseline = extractor->GetBaseline();

        if(saveTmpFiles)
            extractor->save();

        delete extractor;



    //
    // Tractography (using particles filter)
    //

        // Get signal function
        DTFPSignal *signalFun = new DTFPSignal(signal, sigmas, directions, valFileName, baseline, displayMode);

        // Read label image if any and verify sizes
        LabelMapReader::Pointer labelReader = LabelMapReader::New();
        labelReader->SetFileName(labelFilename);
        labelReader->Update();
        LabelMap::Pointer labelVolume = labelReader->GetOutput();


        // Initialize output structures
        std::vector<int> labels;
        std::vector<vtkSmartPointer<vtkAppendPolyData> > fibers;
        std::vector<Image::Pointer> connectMaps;


        // Resample label map
        LabelResampler::Pointer labelResampler = LabelResampler::New();
        labelResampler->SetInput(labelVolume);

        LabelInterpolator::Pointer labelInterpolator = LabelInterpolator::New();
        labelResampler->SetInterpolator(labelInterpolator);

        LabelMap::SpacingType spacing = labelVolume->GetSpacing();
        LabelMap::SizeType size = labelVolume->GetLargestPossibleRegion().GetSize();
        size[0] *= spacing[0]/seedSpacing; size[1] *= spacing[1]/seedSpacing; size[2] *= spacing[2]/seedSpacing;
        labelResampler->SetSize(size);

        LabelMap::DirectionType direction = labelVolume->GetDirection();
        labelResampler->SetOutputDirection(direction);

        LabelMap::PointType origin = labelVolume->GetOrigin();
        labelResampler->SetOutputOrigin(origin);

        spacing[0] = seedSpacing; spacing[1] = seedSpacing; spacing[2] = seedSpacing;
        labelResampler->SetOutputSpacing(spacing);

        labelResampler->Update();
        LabelMap::Pointer resampledLabelVolume = labelResampler->GetOutput();

        // Apply filter on each labeled voxels
        LabelMapIterator labelIt(resampledLabelVolume, resampledLabelVolume->GetLargestPossibleRegion());


        // Apply filter on each labeled voxels
        LabelMapIterator it(labelVolume, labelVolume->GetLargestPossibleRegion());

        for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        {
            int label = (int)it.Get();

            if(label != 0)
            {
                Image::IndexType index = it.GetIndex();

                itk::Point<Real,3> worldPoint;
                labelVolume->TransformIndexToPhysicalPoint(index, worldPoint);

                Mask::IndexType maskIndex;
                mask->TransformPhysicalPointToIndex(worldPoint, maskIndex);

                // If the seed is not in the mask, there is no need to continue this one
                if(mask->GetPixel(maskIndex) != 0)
                {
                    Point begin(worldPoint[0], worldPoint[1], worldPoint[2]);

                    // Set up filter's densities
                    DTFPImportanceDensity importance(signalFun);
                    DTFPAPrioriDensity    apriori(Kappa);
                    DTFPLikelihoodDensity likelihood(signalFun);


                    // Let's start filtering
                    Display1(displayMode, std::cout << "Filtering label " << label << "..." << std::endl);
                    Display2(displayMode, std::cout << "\tSeed's world coordinates: (" << worldPoint[0] << "," << worldPoint[1] << "," << worldPoint[2] << ")" << std::endl);
                    Display2(displayMode, std::cout << "\tSeed's image coordinates: (" << maskIndex[0] << "," << maskIndex[1] << "," << maskIndex[2] << ")" << std::endl);

                    DTFPParticleFilter filter(signalFun, apriori, likelihood, importance, mask, signalFun->getSize(), signalFun->getOrigin(), signalFun->getSpacing(), nbOfParticles, begin, epsilon, stepSize, displayMode);

                    if(lps)
                        filter.SetLPSOn();

                    filter.run(label);


                    // Save data

                    assert(fibers.size() == connectMaps.size());

                    for(int l=fibers.size(); l<label; l++)
                        labels.push_back(-1);

                    if(labels[label-1] == -1)
                    {
                        fibers.push_back(vtkSmartPointer<vtkAppendPolyData>::New());

                        connectMaps.push_back(Image::New());
                        connectMaps.back()->SetOrigin(signalFun->getOrigin());
                        connectMaps.back()->SetSpacing(signalFun->getSpacing());
                        connectMaps.back()->SetRegions(signalFun->getSize());
                        connectMaps.back()->SetDirection(signalFun->GetDirection());
                        connectMaps.back()->Allocate();
                        connectMaps.back()->FillBuffer(0);

                        labels[label-1] = fibers.size()-1;
                    }


                    fibers[labels[label-1]]->AddInput(filter.GetFiber());

                    ImageIterator out(connectMaps[labels[label-1]], connectMaps[labels[label-1]]->GetLargestPossibleRegion());
                    ImageIterator  in(filter.GetConnectionMap(), filter.GetConnectionMap()->GetLargestPossibleRegion());

                    for(in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
                        out.Set(out.Get() + in.Get());

                    Display1(displayMode, std::cout << "done." << std::endl);
                }
            }
        } // for each labeled voxels

        delete signalFun;
        delete directions;
        delete sigmas;


    //
    // Writing files
    //

        // Normalize and write connection map image
        for(unsigned int label = 0; label < labels.size(); label++)
        {
            Real max = 0;
            ImageIterator it(connectMaps[labels[label]], connectMaps[labels[label]]->GetLargestPossibleRegion());

            for(it.GoToBegin(); !it.IsAtEnd(); ++it)
            {
                if(max < it.Get())
                    max = it.Get();
            }

            if(max > 0)
            {
                for(it.GoToBegin(); !it.IsAtEnd(); ++it)
                    it.Set(it.Get() / max);
            }

            try
            {
                std::stringstream filename;
                filename << outMapFileName << "-" << label+1 << ".nii.gz";

                ImageWriter::Pointer writer = ImageWriter::New();
                writer->SetFileName(filename.str().c_str());
                writer->SetInput(connectMaps[labels[label]]);
                writer->Update();
            }
            catch(itk::ImageFileWriterException &err)
            {
                std::cout << "Error: " << std::endl;
                std::cout << err << std::endl;
            }
        }

        // Write fiber polydata
        for(unsigned int label = 0; label < labels.size(); label++)
        {
            unsigned int nbOfInputs = fibers[labels[label]]->GetNumberOfInputPorts();
            bool saveFibers = true;
            unsigned int i = 0;

            do
            {
                if(fibers[labels[label]]->GetNumberOfInputConnections(i++) < 1)
                    saveFibers = false;
            } while(saveFibers && i<nbOfInputs);

            if(saveFibers)
            {
                std::stringstream filename;
                filename << outFibersFileName << "-" << label+1 << ".vtk";

                fibers[label]->Update();
                vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
                writer->SetInput(fibers[labels[label]]->GetOutput());
                writer->SetFileName(filename.str().c_str());
                writer->SetFileTypeToBinary();
                writer->Write();
            }
        }


    return EXIT_SUCCESS;
}
