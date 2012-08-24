/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 05/07/2012
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


// STL includes
#include "cstdlib"
#include "vector"
#include "string"

// TCLAP includes
#include "tclap/CmdLine.h"

// Local includes
#include "btkDiffusionSequenceHelper.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceFileReader.h"
#include "btkDiffusionSequenceFileWriter.h"


int main(int argc, char *argv[])
{
    try
    {
        // Define command line parser
        TCLAP::CmdLine cmd("Diffusion sequence test program.", ' ');

        // Define arguments
        TCLAP::UnlabeledMultiArg< std::string > filenamesArg("filenames", "Input and comparison filenames.", true, "string", cmd);

        // Parse the command line
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::vector< std::string > filenames = filenamesArg.getValue();

        // Check the number of arguments
        if(filenames.size() != 2)
            throw std::string("this test should have two filenames in command line !");


        //
        // Testing
        //

        // Test 1
        std::vector< std::string > test_filenames;
        test_filenames.push_back(filenames[0]);
        std::vector< btk::DiffusionSequence::Pointer > test_reader = btk::DiffusionSequenceHelper::ReadSequenceArray(test_filenames);

        if(test_reader.size() != 1)
            throw std::string("Test 1 failed ! (wrong number of images red)");

        // Test 2
        std::vector< btk::DiffusionSequence::Pointer > test_copy = btk::DiffusionSequenceHelper::CreateNewSequenceFromPhysicalSpaceOf(test_reader);

        if(test_reader[0]->GetLargestPossibleRegion() != test_copy[0]->GetLargestPossibleRegion())
            throw std::string("Test 2 failed ! (regions are different)");

        if(test_reader[0]->GetOrigin() != test_copy[0]->GetOrigin())
            throw std::string("Test 2 failed ! (origins are different)");

        if(test_reader[0]->GetSpacing() != test_copy[0]->GetSpacing())
            throw std::string("Test 2 failed ! (spacings are different)");

        if(test_reader[0]->GetDirection() != test_copy[0]->GetDirection())
            throw std::string("Test 2 failed ! (directions are different)");

        if(test_reader[0]->GetBValues() != test_copy[0]->GetBValues())
            throw std::string("Test 2 failed ! (b-values are different)");

        if(test_reader[0]->GetGradientTable() != test_copy[0]->GetGradientTable())
            throw std::string("Test 2 failed ! (gradient tables are different)");

        // Test 3
        test_filenames[0] = filenames[1];
        btk::DiffusionSequenceHelper::WriteSequenceArray(test_reader, test_filenames);
    }
    catch(itk::ExceptionObject &exception)
    {
        std::cerr << "Exception: " << exception << std::endl;
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
