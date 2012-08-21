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

// Local includes
#include "btkMacro.h"


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
        TCLAP::MultiArg< unsigned int >             labelsArg("l", "label", "Label selection in ROI image", false, "int", cmd);

        TCLAP::SwitchArg particleFilterAlgorithmArg("", "particle_filter", "Use particle filtering algorithm for tractography", true);
        TCLAP::SwitchArg     streamlineAlgorithmArg("", "streamline", "Use streamline algorithm for tractography", false);
        TCLAP::SwitchArg             odfModelingArg("", "odf", "Use ODF modeling for diffusion", true);
        TCLAP::SwitchArg          tensorModelingArg("", "tensor", "Use tensor modeling for diffusion", false);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string    dwiSequenceFileName = dwiSequenceFileNameArg.getValue();
        std::string           maskFileName = maskFileNameArg.getValue();
        std::string            roiFileName = roiFileNameArg.getValue();
        std::vector< unsigned int > labels = labelsArg.getValue();

        bool particleFilterAlgorithm = particleFilterAlgorithmArg.getValue();
        bool     streamlineAlgorithm = streamlineAlgorithmArg.getValue();
        bool             odfModeling = odfModelingArg.getValue();
        bool          tensorModeling = tensorModelingArg.getValue();

        if(particleFilterAlgorithm)
            throw std::string("Particle filter algorithm not yet implemented !");

        if(odfModeling)
            throw std::string("ODF modeling not yet implemented !");


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
