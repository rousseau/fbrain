/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 28/08/2013
  Author(s):Frederic Champ (champ(at)unistra.fr)
  
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
#include "itkExtractImageFilter.h"
#include "itkJoinImageFilter.h"

//Local includes
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkImageHelper.h"
#include "btkFileHelper.h"
#include "btkDiffusionSequenceFileHelper.h"
#include "btkDiffusionDataset.h"
#include "btkWeightedEstimationFilter.h"
#include "btkDwiReconstructionFilter.h"


// Image and sequence definitions


// Filters definitions
typedef btk::DiffusionSequence                          TSequence;
typedef itk::Image<short,3>                             TImage;

typedef itk::ExtractImageFilter< TSequence, TImage >    ExtractImageFilter;

typedef btk::WeightedEstimationFilter                   EstimationType;

typedef btk::DwiReconstructionFilter                    ReconstructionType;

typedef itk::ExtractImageFilter<TSequence,TImage>       ExtractorType;
typedef itk::JoinSeriesImageFilter< TImage,TSequence>   JoinerType;


int main(int argc, char *argv[])
{
    try
    {

        //
        // Parse program's arguments
        //

        // Define command line object for program
        TCLAP::CmdLine cmd("btkReconstruction.", ' ', "1.0");


        // Define arguments
        TCLAP::ValueArg< std::string >  inputFileNameArg ("i", "input", "Input", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >  refFileNameArg ("r", "ref", "ref", false, "", "string", cmd);
        TCLAP::ValueArg< std::string >  outFileNameArg ("o", "output", "output", true, "", "string", cmd);
        TCLAP::ValueArg<bool>           volumeRegistrationArg("","volumeRegistration", "Enable volume to volume registration", false, true, "bool", cmd);
        TCLAP::ValueArg<bool>           sliceRegistrationArg("","sliceRegistration", "Enable slice by slice registration", false, true, "bool", cmd);
        TCLAP::ValueArg<double>         radiusArg("","radius", "Only for WSH: radius of neighbor search", false, false, "double", cmd);


        // Parse command line
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string  inputFileName  = inputFileNameArg.getValue();
        std::string  refFileName    = refFileNameArg.getValue();
        std::string  outputFileName = outFileNameArg.getValue();

        bool   volumeRegistration   = volumeRegistrationArg.getValue();
        bool   sliceRegistration    = sliceRegistrationArg.getValue();
        double radius               = radiusArg.getValue();

        if(!btk::FileHelper::FileExist(inputFileName))
        {
            btkException("Input sequence file is missing or path or name is wrong (name.nii.gz is mandatory) !");

        }
         TSequence::Pointer  inputSequence = btk::DiffusionSequenceHelper::ReadSequence(inputFileName);

        TSequence::RegionType   Region = inputSequence->GetLargestPossibleRegion();
        TSequence::IndexType    Index  = Region.GetIndex();
        TSequence::SizeType     Size   = Region.GetSize();


        std::vector< unsigned short>          BValues       = inputSequence->GetBValues();
        std::vector< btk::GradientDirection > GradientTable = inputSequence->GetGradientTable();


        ReconstructionType::Pointer recons = ReconstructionType::New();
        recons -> SetInputSequence(inputSequence.GetPointer());

        if(refFileNameArg.isSet() )
        {
            if(!btk::FileHelper::FileExist(refFileName))
            {
                btkException("Reference image file is missing or path or name is wrong (name.nii.gz is mandatory) !");
            }
            TImage::Pointer     referenceImage = btk::ImageHelper<TImage>::ReadImage(refFileName);
            recons -> SetReferenceImage(referenceImage.GetPointer());

        }
        recons -> SetVolumeRegistration(volumeRegistration);
        recons -> SetSliceRegistration(sliceRegistration);
        recons -> Update();

        btk::DiffusionDataset::Pointer  Dataset = recons->GetOutput();

        EstimationType::Pointer estimation = EstimationType::New();
        estimation -> SetDiffusionDataset(recons->GetOutput());
        estimation -> SetRadius(radius);
        estimation -> Update();

        TSequence::Pointer EstimatedSequence = estimation->GetOutput();

        JoinerType::Pointer OutJoiner = JoinerType::New();
        OutJoiner -> SetOrigin( 0.0 );
        OutJoiner -> SetSpacing( inputSequence -> GetSpacing()[3]);
        OutJoiner -> SetInput(0,Dataset->GetB0());

        for(unsigned int i=0; i< Size[3]-1; i++)
        {
            TSequence::RegionType  ImageRegion = inputSequence -> GetLargestPossibleRegion();
            ImageRegion.SetSize(3,0);
            ImageRegion.SetIndex(3,i);

            ExtractImageFilter::Pointer extractor = ExtractImageFilter::New();
            extractor -> SetInput(EstimatedSequence);
            extractor -> SetExtractionRegion(ImageRegion);
            extractor -> SetDirectionCollapseToSubmatrix( );
            extractor -> Update();

            OutJoiner -> SetInput(i+1,extractor -> GetOutput());
        }
        OutJoiner -> Update();

        TSequence::Pointer outputSequence = OutJoiner -> GetOutput();

        //////////////////////////////////////////////////////////////////////////
        // Write Outputs : Normalized sequence, B-Values and B-vector
        //////////////////////////////////////////////////////////////////////////
        // Write modified sequence
        btk::ImageHelper<  btk::DiffusionSequence>::WriteImage(outputSequence, outputFileName);

        // Write new B-values with removed zero entries
        std::string bvalFileName = btk::FileHelper::GetRadixOf(outputFileName) + ".bval";
        btk::DiffusionSequenceFileHelper::WriteBValues(BValues,bvalFileName);

        // Write new B-vector with removed zero entries
        std::string bvecFileName = btk::FileHelper::GetRadixOf(outputFileName) + ".bvec";
        btk::DiffusionSequenceFileHelper::WriteGradientTable(GradientTable,bvecFileName);

    }
    catch(itk::ExceptionObject &exception)
    {
        std::cerr << "ITK error:" << std::endl;
        std::cerr << exception << std::endl;
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
