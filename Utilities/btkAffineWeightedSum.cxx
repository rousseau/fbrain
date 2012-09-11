/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 01/12/2011
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
#include "itkAffineTransform.h"

// Local includes
#include "btkMacro.h"
#include "btkIOTransformHelper.h"
#include "btkMatrixOperations.h"


// ITK definitions
typedef itk::AffineTransform< double,3 > AffineTransform;


int main(int argc, char *argv[])
{
    try
    {

        //
        // Command line parser
        //

        // Command line
        TCLAP::CmdLine cmd("Compute a weighted sum of affine transform", ' ', "2.0", true);

        // Arguments
        TCLAP::MultiArg< std::string > inputFileNamesArg("i", "input", "Input affine transform filenames", true, "string", cmd);
        TCLAP::MultiArg< float >              weightsArg("w", "weight", "Image weight (default: 1/N)", false, "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output affine transform filename (default: \"out.txt\")", false, "out.txt", "string", cmd);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();
        std::vector< float >              weights = weightsArg.getValue();
        std::string                outputFileName = outputFileNameArg.getValue();


        //
        // Read & check affine transforms
        //

        itk::TransformFactory< itk::MatrixOffsetTransformBase< double,3,3 > >::RegisterTransform();

        // Verify number of images and number of weights
        if(weights.size() > 0 && inputFileNames.size() != weights.size())
            throw std::string("The number of input images is different than the number of weights !");

        // Read affine transforms
        std::vector< AffineTransform::Pointer > inputTransforms = btk::IOTransformHelper< AffineTransform >::ReadTransform(inputFileNames);

        // Verify if affine transforms have the same center
        AffineTransform::CenterType center = inputTransforms[0]->GetCenter();

        for(unsigned int i = 1; i < inputTransforms.size(); i++)
        {
            if(center != inputTransforms[i]->GetCenter())
                throw std::string("The centers of affine transforms are different !");
        }

        // Set the weight if no weights were given
        float weight = 1.0f / static_cast< float >(inputFileNames.size());

        for(unsigned int i = 0; i < inputFileNames.size(); i++)
        {
            weights.push_back(weight);
        }


        //
        // Processing
        //

        btk::MatrixOperations::Matrix outputTransformMatrix(4,4);
        outputTransformMatrix.Fill(0);

        for(unsigned int i = 0; i < inputTransforms.size(); i++)
        {
            AffineTransform::MatrixType           matrix = inputTransforms[i]->GetMatrix();
            AffineTransform::TranslationType translation = inputTransforms[i]->GetTranslation();

            btk::MatrixOperations::Matrix inputTransformMatrix(4,4);
            inputTransformMatrix(0,0) = matrix(0,0); inputTransformMatrix(0,1) = matrix(0,1); inputTransformMatrix(0,2) = matrix(0,2); inputTransformMatrix(0,3) = translation[0];
            inputTransformMatrix(1,0) = matrix(1,0); inputTransformMatrix(1,1) = matrix(1,1); inputTransformMatrix(1,2) = matrix(1,2); inputTransformMatrix(1,3) = translation[1];
            inputTransformMatrix(2,0) = matrix(2,0); inputTransformMatrix(2,1) = matrix(2,1); inputTransformMatrix(2,2) = matrix(2,2); inputTransformMatrix(2,3) = translation[2];
            inputTransformMatrix(3,0) =           0; inputTransformMatrix(3,1) =           0; inputTransformMatrix(3,2) =           0; inputTransformMatrix(3,3) =              1;

            outputTransformMatrix += btk::MatrixOperations::Logarithm(inputTransformMatrix) * weights[i];
        }

        outputTransformMatrix = btk::MatrixOperations::Exponential(outputTransformMatrix);


        //
        // Write output affine transform
        //

        AffineTransform::Pointer outputTransform = AffineTransform::New();
        outputTransform->SetIdentity();

        AffineTransform::MatrixType matrix;
        matrix(0,0) = outputTransformMatrix(0,0); matrix(0,1) = outputTransformMatrix(0,1); matrix(0,2) = outputTransformMatrix(0,2);
        matrix(1,0) = outputTransformMatrix(1,0); matrix(1,1) = outputTransformMatrix(1,1); matrix(1,2) = outputTransformMatrix(1,2);
        matrix(2,0) = outputTransformMatrix(2,0); matrix(2,1) = outputTransformMatrix(2,1); matrix(2,2) = outputTransformMatrix(2,2);

        AffineTransform::TranslationType translation;
        translation[0] = outputTransformMatrix(0,3); translation[1] = outputTransformMatrix(1,3); translation[2] = outputTransformMatrix(2,3);

        outputTransform->SetCenter(inputTransforms[0]->GetCenter());
        outputTransform->SetMatrix(matrix);
        outputTransform->SetTranslation(translation);

        btk::IOTransformHelper< AffineTransform >::WriteTransform(outputTransform, outputFileName);
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





