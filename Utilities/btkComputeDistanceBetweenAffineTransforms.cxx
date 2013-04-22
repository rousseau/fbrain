/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 17/04/2013
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
        TCLAP::CmdLine cmd("Compute the distance between two affine transforms", ' ', "1.0", true);

        // Arguments
        TCLAP::UnlabeledMultiArg< std::string > inputFileNamesArg("input", "Input affine transform filenames", true, "string", cmd);

        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();


        //
        // Read & check affine transforms
        //

        itk::TransformFactory< itk::MatrixOffsetTransformBase< double,3,3 > >::RegisterTransform();

        // Verify number of images and number of weights
        if(inputFileNames.size() != 2)
            throw std::string("Main: There should be 2 affine transform !");

        // Read affine transforms
        std::vector< AffineTransform::Pointer > inputTransforms = btk::IOTransformHelper< AffineTransform >::ReadTransform(inputFileNames);

        // Verify if affine transforms have the same center
        AffineTransform::CenterType center = inputTransforms[0]->GetCenter();

        for(unsigned int i = 1; i < inputTransforms.size(); i++)
        {
            if(center != inputTransforms[i]->GetCenter())
                throw std::string("The centers of affine transforms are different !");
        }


        //
        // Processing
        //

        // Compute euclidian distance in log-euclidian space
        btk::MatrixOperations::Matrix M1(4,4);
        M1.Fill(0);

        AffineTransform::MatrixType           matrix1 = inputTransforms[0]->GetMatrix();
        AffineTransform::TranslationType translation1 = inputTransforms[0]->GetTranslation();

        M1(0,0) = matrix1(0,0); M1(0,1) = matrix1(0,1); M1(0,2) = matrix1(0,2); M1(0,3) = translation1[0];
        M1(1,0) = matrix1(1,0); M1(1,1) = matrix1(1,1); M1(1,2) = matrix1(1,2); M1(1,3) = translation1[1];
        M1(2,0) = matrix1(2,0); M1(2,1) = matrix1(2,1); M1(2,2) = matrix1(2,2); M1(2,3) = translation1[2];
        M1(3,0) =            0; M1(3,1) =            0; M1(3,2) =            0; M1(3,3) =               1;

        btk::MatrixOperations::Matrix M2(4,4);
        M2.Fill(0);

        AffineTransform::MatrixType           matrix2 = inputTransforms[1]->GetMatrix();
        AffineTransform::TranslationType translation2 = inputTransforms[1]->GetTranslation();

        M2(0,0) = matrix2(0,0); M2(0,1) = matrix2(0,1); M2(0,2) = matrix2(0,2); M2(0,3) = translation2[0];
        M2(1,0) = matrix2(1,0); M2(1,1) = matrix2(1,1); M2(1,2) = matrix2(1,2); M2(1,3) = translation2[1];
        M2(2,0) = matrix2(2,0); M2(2,1) = matrix2(2,1); M2(2,2) = matrix2(2,2); M2(2,3) = translation2[2];
        M2(3,0) =            0; M2(3,1) =            0; M2(3,2) =            0; M2(3,3) =               1;

        btk::MatrixOperations::Matrix M(4,4);
        M.Fill(0);
        M = btk::MatrixOperations::Logarithm(M1) - btk::MatrixOperations::Logarithm(M2);
//        M = btk::MatrixOperations::Exponential(M);

        std::cerr << btk::MatrixOperations::Norm(M) << std::endl;
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
