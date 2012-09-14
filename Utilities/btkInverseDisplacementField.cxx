/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 13/09/2012
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
#include "string"
#include "vector"

// TCLAP includes
#include <tclap/CmdLine.h>

// ITK includes
#include "itkDisplacementFieldTransform.h"
#include "itkInverseDisplacementFieldImageFilter.h"

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"


// ITK definitions
typedef itk::DisplacementFieldTransform< float,3 >::DisplacementFieldType DeformationField;


// Main function
int main(int argc, char *argv[])
{
    try
    {

        //
        // Command line parser
        //

        // Command line
        TCLAP::CmdLine cmd("Compute the inverse of a displacement field", ' ', "1.0", true);

        // Arguments
        TCLAP::ValueArg< std::string > inputFileNameArg("i", "input", "Input field filenames", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output field filename (default: \"out.nii.gz\")", false, "out.nii.gz", "string", cmd);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        std::string  inputFileName = inputFileNameArg.getValue();
        std::string outputFileName = outputFileNameArg.getValue();


        //
        // Reading
        //

        DeformationField::Pointer input = btk::ImageHelper< DeformationField >::ReadImage(inputFileName);


        //
        // Processing
        //

        typedef itk::InverseDisplacementFieldImageFilter< DeformationField, DeformationField > InverseDeformationFieldFilter;

        btkCoutMacro("Processing...");

        InverseDeformationFieldFilter::Pointer filter = InverseDeformationFieldFilter::New();
        filter->SetInput(input);
        filter->SetOutputSpacing(input->GetSpacing());
        filter->SetOutputOrigin(input->GetOrigin());
        filter->SetSize(input->GetLargestPossibleRegion().GetSize());
        filter->SetSubsamplingFactor(16);
        filter->Update();

        btkCoutMacro("done.");


        //
        // Writing
        //

        btk::ImageHelper< DeformationField >::WriteImage(filter->GetOutput(), outputFileName);

    }
    catch(TCLAP::ArgException &e)
    {
        btkCoutMacro("Exception: " << e.error() << " for arg " << e.argId());
    }
    catch(std::string &message)
    {
        btkCoutMacro("Exception: " << message);
    }

  return EXIT_SUCCESS;
}
