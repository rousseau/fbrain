/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

31 march 2010
< pontabry at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


// TCLAP : Templatized C++ Command Line Parser includes
#include "CmdLine.h"

// STL includes
#include "string"
#include "iostream"

// Local includes
#include "btkTypes.h"
#include "btkPoint.h"
#include "btkSignal.h"
#include "btkSHModel.h"
#include "btkSHModelDensity.h"
#include "btkNormalDensity.h"
#include "btkVonMisesFisherDensity.h"
#include "btkImportanceDensity.h"
#include "btkInitialDensity.h"
#include "btkAPrioriDensity.h"
#include "btkLikelihoodDensity.h"
#include "btkParticleFilter.h"


int main(int argc, char *argv[])
{
    try
    {
        //
        // Program's command line parser definition
        //

            // Command line variables
            std::string modelFileName;
            std::string signalFileName;
            std::string dirFileName;
            std::string sigmasFileName;
            std::string maskFileName;
            unsigned int nbOfParticles;
            btk::Real stepSize;
            btk::Real epsilon;
            btk::Real xBegin;
            btk::Real yBegin;
            btk::Real zBegin;
            std::string labelFilename;
            btk::Real Kappa;
            unsigned int maxLength;
            btk::Real angleThreshold;

            // Defines command line parser
            TCLAP::CmdLine cmd("BTK Particle filtering", ' ', "0.6");

            // Defines arguments
            TCLAP::ValueArg<std::string>  modelArg("m", "model", "Estimated model (.nii)", true, "", "string");
            TCLAP::ValueArg<std::string>  signalArg("s", "signal", "Diffusion Signal (.nii)", true, "", "string");
            TCLAP::ValueArg<std::string>  dirArg("d", "directions", "Gradient directions", true, "", "string");
            TCLAP::ValueArg<std::string>  sigArg("a", "sigmas", "Standard deviation of signal images", true, "", "string");
            TCLAP::ValueArg<unsigned int> particlesArg("p", "particles", "Number of particles", false, 1000, "unsigned int");
            TCLAP::ValueArg<btk::Real>    epsilonArg("e", "epsilon", "Resampling treshold", false, 0.6, "btk::Real");
            TCLAP::ValueArg<btk::Real>    stepSizeArg("t", "stepSize", "Step size", false, 0.5, "btk::Real");
            TCLAP::ValueArg<btk::Real>    xBeginArg("x", "x-coord", "Begin's x-coord", false, -1, "btk::Real");
            TCLAP::ValueArg<btk::Real>    yBeginArg("y", "y-coord", "Begin's y-coord", false, -1, "btk::Real");
            TCLAP::ValueArg<btk::Real>    zBeginArg("z", "z-coord", "Begin's z-coord", false, -1, "btk::Real");
            TCLAP::ValueArg<std::string>  maskArg("b", "mask", "Image mask", true, "", "string");
            TCLAP::ValueArg<std::string>  labelArg("l", "label", "Labels image", false, "", "string");
            TCLAP::ValueArg<btk::Real>    KappaArg("c", "concentration", "Vmf concentration", false, 30.0, "btk::Real");
            TCLAP::ValueArg<unsigned int> maxLengthArg("k", "max-length", "Maximal length of particles", true, 0, "unsigned int");
            TCLAP::ValueArg<btk::Real>    angleThreshArg("f", "angle-threshold", "Angle threshold", false, M_PI/3., "btk::Real");

            cmd.add(modelArg);
            cmd.add(signalArg);
            cmd.add(dirArg);
            cmd.add(sigArg);
            cmd.add(particlesArg);
            cmd.add(epsilonArg);
            cmd.add(stepSizeArg);
            cmd.add(xBeginArg);
            cmd.add(yBeginArg);
            cmd.add(zBeginArg);
            cmd.add(maskArg);
            cmd.add(labelArg);
            cmd.add(KappaArg);
            cmd.add(maxLengthArg);
            cmd.add(angleThreshArg);

            // Parsing arguments
            cmd.parse(argc, argv);

            // Get back arguments' values
            modelFileName  = modelArg.getValue();
            signalFileName = signalArg.getValue();
            dirFileName    = dirArg.getValue();
            sigmasFileName = sigArg.getValue();
            nbOfParticles  = particlesArg.getValue();
            epsilon        = epsilonArg.getValue();
            stepSize       = stepSizeArg.getValue();
            xBegin         = xBeginArg.getValue();
            yBegin         = yBeginArg.getValue();
            zBegin         = zBeginArg.getValue();
            maskFileName   = maskArg.getValue();
            labelFilename  = labelArg.getValue();
            Kappa          = KappaArg.getValue();
            maxLength      = maxLengthArg.getValue();
            angleThreshold = angleThreshArg.getValue();


            // verify entries
            if(labelFilename.empty())
            {
                if(xBegin == -1 && yBegin == -1 && zBegin == -1)
                {
                    std::cout << "Error: You must specify at least a label image or cartesian coordinates for starting point!" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }


        //
        // Tractography (using particles filter)
        //


            using namespace btk;


            // Get spherical harmonics model
            SHModel *model = new SHModel(modelFileName);

            // Get signal
            Signal *signal = new Signal(signalFileName, sigmasFileName, dirFileName);

            // Read label image if any and verify sizes
            Image::Pointer labels;
            if(!labelFilename.empty())
            {
                ImageReader::Pointer reader = ImageReader::New();
                reader->SetFileName(labelFilename);
                reader->Update();
                labels = reader->GetOutput();
            }

            // These densities will be used next
            SHModelDensity modelDensity(model);
            NormalDensity normalDensity;
            VonMisesFisherDensity vmf(Kappa);
            VonMisesFisherDensity testvmfimp(100);


            // Apply filter on each labeled voxels
            if(!labelFilename.empty())
            {
                ImageIterator it(labels, labels->GetLargestPossibleRegion());

                for(it.GoToBegin(); !it.IsAtEnd(); ++it)
                {
                    int label = (int)it.Get();

                    if(label != 0)
                    {
                        Image::IndexType index = it.GetIndex();
                        itk::Point<Real,3> worldPoint;
                        labels->TransformIndexToPhysicalPoint(index,worldPoint);
                        Point begin(worldPoint[0], worldPoint[1], worldPoint[2]);

                        // Set up filter's densities
                        std::cout << "Setting up needed densities for label " << label << begin << " ..." << std::endl;
                        ImportanceDensity importance(testvmfimp, model, angleThreshold);
                        InitialDensity    initial(modelDensity, begin);
                        APrioriDensity    apriori(vmf);
                        LikelihoodDensity likelihood(normalDensity, signal, model);
                        std::cout << "done." << std::endl;


                        // Let's start filtering
                        std::cout << "Filtering label " << label << "..." << std::endl;

                        ParticleFilter filter(model, initial, apriori, likelihood, importance, maskFileName,
                            signal->getSize(), signal->getOrigin(), signal->getSpacing(),
                            nbOfParticles, begin, epsilon, stepSize, maxLength);

                        filter.run(label);

                        std::cout << "done." << std::endl;
                    }
                } // for each labeled voxels
            }
            else
            {
                int label = 0;
                Point begin(xBegin, yBegin, zBegin);

                // Set up filter's densities
                std::cout << "Setting up needed densities for label " << label << " " << begin << " ..." << std::endl;
                ImportanceDensity importance(testvmfimp, model, angleThreshold);
                InitialDensity    initial(modelDensity, begin);
                APrioriDensity    apriori(vmf);
                LikelihoodDensity likelihood(normalDensity, signal, model);
                std::cout << "done." << std::endl;


                // Let's start filtering
                std::cout << "Filtering label " << label << "..." << std::endl;

                ParticleFilter filter(model, initial, apriori, likelihood, importance, maskFileName,
                        signal->getSize(), signal->getOrigin(), signal->getSpacing(),
                        nbOfParticles, begin, epsilon, stepSize, maxLength);

                filter.run(label);
                filter.saveConnectionMap(label,begin);
                filter.saveFiber(label, 0, begin);

                std::cout << "done." << std::endl;
            }


        //
        // Clean
        //

            delete model;
            delete signal;
    }
    catch(TCLAP::ArgException &e)
    {
        std::cout << "TCLAP error: " << e.error() << " for argument " << e.argId() << std::endl;
    }

    return EXIT_SUCCESS;
}
