/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 februar 2010
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

// Local includes
#include "btkSphericalHarmonics.h"
#include "btkSHModelEstimator.h"


int main(int argc, char *argv[])
{
    try
    {
        //
        // Program's command line parser definition
        //

            // Command line variables
            std::string signalFileName;
            std::string directionsFileName;
            std::string maskFileName;
            unsigned int order;

            // Defines command line parser
            TCLAP::CmdLine cmd("LSIIT-MIV Model estimator", ' ', "0.5");

            // Defines arguments
            TCLAP::ValueArg<std::string>  signalArg("s", "signal", "Signal image", true, "", "string");
            TCLAP::ValueArg<std::string>  directionsArg("d", "directions", "Gradient directions", true, "", "string");
            TCLAP::ValueArg<unsigned int> orderArg("o", "order", "Order of computation", false, 4, "unsigned int");
            TCLAP::ValueArg<std::string> maskArg("b", "mask", "Image mask", true, "", "string");

            cmd.add(signalArg);
            cmd.add(directionsArg);
            cmd.add(orderArg);
            cmd.add(maskArg);

            // Parsing arguments
            cmd.parse(argc, argv);

            // Get back arguments' values
            signalFileName     = signalArg.getValue();
            directionsFileName = directionsArg.getValue();
            order              = orderArg.getValue();
            maskFileName       = maskArg.getValue();


        //
        // Local model estimation
        //

            btk::SHModelEstimator *estimator =
                    new btk::SHModelEstimator(signalFileName, directionsFileName, maskFileName, order);
            estimator->estimate();
            estimator->save();

            delete estimator;
    }
    catch(TCLAP::ArgException &e)
    {
        std::cout << "Error: " << e.error() << " for argument " << e.argId() << std::endl;
    }

    return EXIT_SUCCESS;
}
