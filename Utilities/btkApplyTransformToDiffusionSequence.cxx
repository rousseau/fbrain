/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/08/2013
  Author(s): Marc Schweitzer (marc.schweitzer (at) unistra.fr)

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
#include "iostream"
#include "string"
#include "sstream"

// TCLAP includes
#include "tclap/CmdLine.h"

// ITK includes
#include "itkImage.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkTransform.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"

// Local includes
#include "btkFileHelper.h"
#include "btkImageHelper.h"
#include "btkIOTransformHelper.h"
#include "btkCommandIterationUpdate.h"
#include "btkMacro.h"
#include "btkDiffusionGradientTable.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"


// Image and sequence definitions
typedef short                                   PixelType;
typedef double                                  PrecisionType;
typedef itk::Image< PixelType,3 >               Image;
typedef btk::DiffusionSequence                  Sequence;
typedef btk::DiffusionGradientTable< Sequence > GradientTable;
typedef itk::Transform< PrecisionType, 3 >      itkTransform;
typedef itk::MatrixOffsetTransformBase<PrecisionType,3,3> MatrixTransformType;
typedef itk::Euler3DTransform< PrecisionType >  EulerTransform;
// Filters definitions
typedef itk::ExtractImageFilter< Sequence,Image >                      SequenceExtractor;
typedef itk::JoinSeriesImageFilter< Image,Sequence >                   ImageSequenceJoiner;
typedef itk::ResampleImageFilter< Image,Image >                        ImageResampler;
typedef itk::BSplineInterpolateImageFunction< Image,PrecisionType,PrecisionType >   BSplineInterpolator;

int main( int argc, char *argv[])
{
    try
    {
        //
        // Parse program's arguments
        //

        // Define command line object for program
        TCLAP::CmdLine cmd("Registers diffusion sequence to anatomical data.", ' ', "2.0");

        // Define arguments
        TCLAP::ValueArg< std::string > inputSequenceFileNameArg("d", "diffusion_sequence", "Input diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputSequenceFileNameArg("o", "output", "Registred output diffusion sequence", true, "", "string", cmd);
        TCLAP::SwitchArg invArg  ("","inv","inverse the transformation", cmd, false);
        TCLAP::ValueArg< std::string > transformFileNameArg("t", "transformation", "Transformation to apply", true, "", "string", cmd);

        //TCLAP::SwitchArg verboseModeArg("v", "verbose", "Verbose mode", cmd, false);

        // Parse command line
        cmd.parse(argc, argv);
        // Get back arguments' values
        std::string  inputSequenceFileName = inputSequenceFileNameArg.getValue();
        std::string outputSequenceFileName = outputSequenceFileNameArg.getValue();
        std::string transformationFileName = transformFileNameArg.getValue();

        bool inverseTheTransform = invArg.getValue();

        // Set B-values and gradients' filenames according to diffusion sequences' filenames
        std::string     inputSequenceFileNameRadix = btk::FileHelper::GetRadixOf(inputSequenceFileName);
        std::string inputSequenceGradientsFileName = inputSequenceFileNameRadix + ".bvec";
        std::string inputSequenceBValuesFileName = inputSequenceFileNameRadix + ".bval";

        std::string     outputSequenceFileNameRadix = btk::FileHelper::GetRadixOf(outputSequenceFileName);
        std::string outputSequenceGradientsFileName = outputSequenceFileNameRadix + ".bvec";

        Sequence::Pointer inputSequence = btk::DiffusionSequenceHelper::ReadSequence(inputSequenceFileName);


        itk::TransformFactory< MatrixTransformType >::RegisterTransform();
        MatrixTransformType::Pointer transform = btk::IOTransformHelper< MatrixTransformType >::ReadTransform(transformationFileName);

        MatrixTransformType::Pointer inverseTransform = MatrixTransformType::New();

        transform->GetInverse(inverseTransform);

        if(inverseTheTransform)
        {
           transform = inverseTransform ;
        }

        // Extract baseline image from diffusion sequence
        Sequence::RegionType extractionRegion = inputSequence->GetLargestPossibleRegion();
        extractionRegion.SetSize(3, 0);

        SequenceExtractor::Pointer extractor = SequenceExtractor::New();
        extractor->SetInput(inputSequence);
        extractor->SetExtractionRegion(extractionRegion);
        extractor->SetDirectionCollapseToSubmatrix();
        extractor->Update();


        std::cout<<"Warping diffusion sequence ..."<<std::endl;

        // Define image sequence joiner
        ImageSequenceJoiner::Pointer joiner = ImageSequenceJoiner::New();
        joiner->SetOrigin(0);
        joiner->SetSpacing(1);

        // Each image of diffusion sequence is warped to anatomical reference and resampled
        for(unsigned int i = 0; i < inputSequence->GetLargestPossibleRegion().GetSize(3); i++)
        {
            // Define the current region of interest (1 image is sequence) for extractor
            Sequence::RegionType currentRegion = inputSequence->GetLargestPossibleRegion();
            currentRegion.SetSize(3,0);
            currentRegion.SetIndex(3,i);

            extractor->SetExtractionRegion(currentRegion);
            extractor->Update();

            // Define resampler
            ImageResampler::Pointer resampler = ImageResampler::New();
            resampler->SetTransform(transform);
            resampler->SetInput(extractor->GetOutput());
            resampler->SetUseReferenceImage(true);
            resampler->SetReferenceImage(extractor->GetOutput());
            resampler->SetDefaultPixelValue(0);

            // Define interpolation
            BSplineInterpolator::Pointer interpolator = BSplineInterpolator::New();
            interpolator->SetSplineOrder(3);
            resampler->SetInterpolator(interpolator);

            // Resample image and add to output sequence
            resampler->Update();
            joiner->SetInput(i, resampler->GetOutput());
        }

        joiner->Update();

        Sequence::Pointer outputSequence = joiner->GetOutput();

        // Set diffusion information
        outputSequence->SetBValues(inputSequence->GetBValues());
        outputSequence->SetGradientTable(inputSequence->GetGradientTable());

        std::cout << "done." << std::endl;


        //
        // Update gradient table
        //

        std::cout << "Updating gradient table..." << std::endl;

        // Decompose affine transform into rigid and affine part
        // (only rigid part is necessary for gradient table reorientation)
        vnl_matrix< PrecisionType > R(3,3);
        R.set_identity();
        R = transform->GetMatrix().GetVnlMatrix();

        vnl_matrix< PrecisionType > PQ = R;
        vnl_matrix< PrecisionType > NQ = R;
        vnl_matrix< PrecisionType > PQNQDiff;

        for(unsigned int ni = 0; ni < 100; ni++ )
        {
            // Average current Qi with its inverse transpose
            NQ = ( PQ + vnl_inverse_transpose( PQ ) ) / 2.0;
            PQNQDiff = NQ - PQ;

            if( PQNQDiff.frobenius_norm() < 1e-7 )
            {
                break;
            }
            else
            {
                PQ = NQ;
            }
        }

        // Define rigid transformation for gradient table reorientation
        EulerTransform::Pointer gradientTableTransform = EulerTransform::New();
        gradientTableTransform->SetMatrix(NQ);

        // Load input gradient table, rotate and save to new gradient table file
        GradientTable::Pointer inputGradientTable = GradientTable::New();
        inputGradientTable->SetNumberOfGradients(inputSequence->GetLargestPossibleRegion().GetSize(3));
        inputGradientTable->SetImage(inputSequence);
        inputGradientTable->SetTransform(gradientTableTransform);
        inputGradientTable->LoadFromFile(inputSequenceGradientsFileName.c_str());
        inputGradientTable->RotateGradientsInWorldCoordinates();

        std::cout << "done." << std::endl;

        //
        // Write output sequence
        //

        // Save diffusion sequence
        btk::DiffusionSequenceHelper::WriteSequence(outputSequence, outputSequenceFileName);

        // Save gradient table
        inputGradientTable->SaveToFile(outputSequenceGradientsFileName.c_str());

    }
    catch(itk::ExceptionObject &exception)
    {
        std::cerr << "ITK error:" << std::endl;
        std::cerr << exception << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
