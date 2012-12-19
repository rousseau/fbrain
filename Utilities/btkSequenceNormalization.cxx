/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 03/12/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr), Frederic Champ (champ@unistra.fr)

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

/* Standard includes */
#include "string"
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkJoinSeriesImageFilter.h"

/* Btk includes */
#include "btkAffineRegistration.h"
#include "btkDiffusionGradientTable.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkDiffusionSequenceFileHelper.h"
#include "btkFileHelper.h"
#include "btkWeightedSumOfImagesFilter.h"

int main( int argc, char * argv[] )
{
    try{

        //
        // Command line parser
        //

        TCLAP::CmdLine cmd("Writes a dwi sequence as a single image B0 + the diffusion "
                           "images. The new B0 is the mean of all B0 images in the original sequence, "
                           "or a user-provided B0.", ' ', "Unversioned");

        TCLAP::ValueArg<std::string> inputArg("i","input","Original DWI sequence.", true,"","string",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Normalized DWI Sequence. ", true,"","string",cmd);
        TCLAP::ValueArg<std::string> maskArg("m","mask","Image mask for registration. ", false,"","string",cmd);

        cmd.parse( argc, argv );

        std::string inputFileName  = inputArg.getValue();
        std::string outputFileName = outputArg.getValue();
        std::string maskFileName      = maskArg.getValue();

        //
        // Typedefs
        //
        typedef btk::DiffusionSequence TSequence;
        typedef itk::Image< double,3> TImage;
        typedef itk::Image< unsigned char,3> TMask;

        typedef btk::AffineRegistration<TImage> AffineRegistration;
        typedef itk::ExtractImageFilter< TSequence, TImage > ExtractImageFilter;
        typedef itk::JoinSeriesImageFilter< TImage, TSequence > JoinerType;
        typedef btk::ResampleImageFilter< TImage,TImage,double > ResampleImageFilter;
        typedef btk::WeightedSumOfImagesFilter< TImage > WeightedSumOfImagesFilter;

        // TODO : Do a filter or add a normalization method to btkDiffusionSequence class...

        //
        // Processing
        //

        /////////////////////////////////////////////////////////////
        //
        // Read sequence
        //
        TSequence::Pointer    sequence    = btk::DiffusionSequenceHelper::ReadSequence(inputFileName);
        TSequence::RegionType region4D    = sequence -> GetLargestPossibleRegion();
        TSequence::SizeType   input4Dsize = region4D.GetSize();
        TSequence::IndexType  index       = region4D.GetIndex();
        TSequence::SizeType   image3Dsize = input4Dsize;
        image3Dsize[3]=0;

        std::vector< unsigned short>          BValues       = sequence->GetBValues();
        std::vector< btk::GradientDirection > GradientTable = sequence->GetGradientTable();
        std::vector< TImage::Pointer > B0;
        std::vector< TImage::Pointer > B0_resampled;

        /////////////////////////////////////////////////////////////
        //
        // Extraction of B0s
        //
        btkCoutMacro("B0s extraction :");

        unsigned int nb_loop =0;

        for (unsigned int i = 0; i < input4Dsize[3]; i++)
        {
            if (BValues[i] == 0)
            {
                index[3] = i;

                TSequence::RegionType imageRegion;
                imageRegion.SetIndex( index );
                imageRegion.SetSize( image3Dsize );

                ExtractImageFilter::Pointer extractor = ExtractImageFilter::New();
                extractor -> SetInput( sequence );
                extractor -> SetExtractionRegion( imageRegion );
                extractor -> SetDirectionCollapseToSubmatrix( );
                extractor -> Update();
                B0.push_back( extractor ->  GetOutput() );

                nb_loop++;
                std::cout<<"\r -> "<< nb_loop<<" volumes extracted ... ";
            }
        }
        btkCoutMacro("Done.");

        /////////////////////////////////////////////////////////////
        //
        // Registration of B0s
        //
        btkCoutMacro("B0s registration on the first one :");

        // Read mask
        TMask::Pointer  mask;

        if (maskFileName !="")
        {
            mask = btk::ImageHelper< TMask >::ReadImage(maskFileName);
        }

        B0_resampled.push_back(B0[0]);

        nb_loop =0;
        unsigned int i;

        #pragma omp parallel for private(i) schedule(dynamic)
        for(unsigned int i=1; i<B0.size(); i++ )
        {
            AffineRegistration::Pointer registration = AffineRegistration::New();
            ResampleImageFilter::Pointer resampler = ResampleImageFilter::New();

            if (maskFileName != "")
            {
                registration -> SetFixedImageMask( mask );
            }

            registration -> SetFixedImage( B0[0] );
            registration -> SetMovingImage(B0[i] );
            registration -> SetFixedImageRegion( B0[0] -> GetLargestPossibleRegion() );
            registration -> Update();

            resampler -> UseReferenceImageOn();
            resampler -> SetReferenceImage(B0[0]);
            resampler -> SetInput(B0[i]);
            resampler -> SetTransform( registration -> GetTransform() );
            resampler -> Update();

            B0_resampled.push_back( resampler -> GetOutput() );

            nb_loop++;
            std::cout<<"\r -> "<<nb_loop <<" volumes registred (total "<<B0.size()-1<<") ... "<<std::flush;
        }
        btkCoutMacro("Done.");

        /////////////////////////////////////////////////////////////
        //
        // Calculates the mean image
        //
        std::cout<<"Calculating the mean image ... ";

        std::vector< float > weights = std::vector< float >(B0_resampled.size(), 1.0f/static_cast< float >(B0_resampled.size()));

        WeightedSumOfImagesFilter::Pointer sumFilter = WeightedSumOfImagesFilter::New();
        sumFilter->SetInputs(B0_resampled);
        sumFilter->SetWeights(weights);
        sumFilter->Update();

        TImage::Pointer B0_mean = sumFilter->GetOutput();

        // Correcting of B-Values and B-vector (Removes repetitive zero entries)
        while (BValues[1]==0)
        {
            BValues.erase(BValues.begin()+1);
            GradientTable.erase(GradientTable.begin()+1);
        }
        btkCoutMacro("Done.");

        ////////////////////////////////////////////////////////////
        //
        // Join images
        //
        std::cout << "Joining images ... ";

        TSequence::SpacingType  spacing  = sequence -> GetSpacing();

        JoinerType::Pointer joiner = JoinerType::New();
        joiner -> SetOrigin( 0.0 );
        joiner -> SetSpacing( spacing[3] );
        joiner -> SetInput( 0, B0_mean );

        unsigned j = 1;
        for (unsigned int i = B0_resampled.size(); i < input4Dsize[3]; i++)
        {
            index[3] = i;

            TSequence::RegionType imageRegion;
            imageRegion.SetIndex( index );
            imageRegion.SetSize( image3Dsize );

            ExtractImageFilter::Pointer extractor = ExtractImageFilter::New();
            extractor -> SetInput( sequence );
            extractor -> SetExtractionRegion( imageRegion );
            extractor -> SetDirectionCollapseToSubmatrix( );
            extractor -> Update();

            joiner -> SetInput( j, extractor -> GetOutput() );
            j++;
        }
        joiner -> Update();
        btkCoutMacro(" Done.");

        //////////////////////////////////////////////////////////////////////////
        ///
        // Write Outputs : Normalized sequence, B-Values and B-vector
        //

        // Write modified sequence
        btk::ImageHelper<  btk::DiffusionSequence>::WriteImage(joiner->GetOutput(), outputFileName);

        // Write new B-values with removed zero entries
        std::string bvalFileName = btk::FileHelper::GetRadixOf(outputFileName) + ".bval";
        btk::DiffusionSequenceFileHelper::WriteBValues(BValues,bvalFileName);

        // Write new B-vector with removed zero entries
        std::string bvecFileName = btk::FileHelper::GetRadixOf(outputFileName) + ".bvec";
        btk::DiffusionSequenceFileHelper::WriteGradientTable(GradientTable,bvecFileName);


    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

    return EXIT_SUCCESS;
}
