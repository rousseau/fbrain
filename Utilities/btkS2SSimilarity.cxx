/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date:
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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
#include "tclap/CmdLine.h"

// ITK includes
#include "itkImage.h"
// Local includes
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkS2SSimilarityFilter.h"


// Image and sequence definitions
typedef short                                   PixelType;

typedef itk::Image< PixelType,3 >               TImage;

typedef btk::DiffusionSequence                  TSequence;
typedef btk::DiffusionSequence::Pointer         SequencePointer;

typedef std::vector< float  >                   S2SVectorType;
typedef std::vector<std::vector< float  > >     SimilarityS2SType;

// Filters definitions
typedef btk::S2SSimilarityFilter< TImage >                      S2SSimilarityType;



//----------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    try
    {

        //////////////////////////////////////////////////////////////////////////
        //
        // Parse program's arguments
        //

        // Define command line object for program
        TCLAP::CmdLine cmd("DWI volumes similarity.", ' ', "1.0");

        // Define arguments
        TCLAP::ValueArg< std::string >  inputSequenceFileNameArg("i", "sequence", "Input diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >  outputFileNameArg("f", "similarities", "Similarities file", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >  methodArg("m", "method", "method to calculate similarities: NonNormalized"
                                                  " or Normalized1 ", false, "Normalized", "string", cmd);
        TCLAP::ValueArg< std::string >  delimiterArg("d", "delimiter", "delimiter", false, ";", "string", cmd);

        TCLAP::ValueArg<int> sliceNumberArg("","sliceNumber", "Slice for which the similarity will be calculated", false, 0, "int", cmd);
        TCLAP::ValueArg<int> volumeNumberArg("","volumeNumber","Volume for which the similarity will be calculated", false,0,"int",cmd);

        TCLAP::SwitchArg verboseModArg("v", "verbose", "Verbose mode", cmd, false);

        // Parse command line
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string inputSequenceFileName = inputSequenceFileNameArg.getValue();
        std::string outputFileName        = outputFileNameArg.getValue();
        std::string method                = methodArg.getValue();
        std::string delimiter             = delimiterArg.getValue();

        unsigned int sliceNumber          = sliceNumberArg.getValue();
        unsigned int volumeNumber         = volumeNumberArg.getValue();

        bool verboseMod = verboseModArg.getValue();

        //////////////////////////////////////////////////////////////////////////
        //
        // Read images
        //
        // Test if mandatory files exists
        if(!btk::FileHelper::FileExist(inputSequenceFileName))
        {
            btkException("Input sequence file are missing or path or name is wrong (name.nii.gz is mandatory) !");

        }

        // Read input images
        SequencePointer inputSequence = btk::DiffusionSequenceHelper::ReadSequence(inputSequenceFileName);

        //////////////////////////////////////////////////////////////////////////
        //
        // Calculate similarities
        //
        S2SSimilarityType::Pointer SimilarityFilter = S2SSimilarityType::New();
        SimilarityFilter -> SetInputSequence(inputSequence.GetPointer());
        SimilarityFilter -> SetVerboseMod(verboseMod);
        SimilarityFilter -> SetDelimiter(delimiter);
        if(sliceNumberArg.isSet())
        {
            SimilarityFilter -> SetSliceNumber(sliceNumber);
        }
        SimilarityFilter -> SetVolumeNumber(volumeNumber);       
        SimilarityFilter -> SetMethod(method);
        SimilarityFilter -> Update();


        //////////////////////////////////////////////////////////////////////////
        //
        // Write data file
        //
        SimilarityFilter -> WriteData(outputFileName);

    }
    catch(itk::ExceptionObject &exception)
    {
        std::cerr << "ITK error:" << std::endl;
        std::cerr << exception << std::endl;
        exit(EXIT_FAILURE);
    }
}


