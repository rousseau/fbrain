/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/06/2013
  Author(s):Frederic CHAMP (champ(at)unistra.fr)
  
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

// Local includes
#include "btkFileHelper.h"
#include "btkImageHelper.h"
#include "btkIOTransformHelper.h"
#include "btkMacro.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkDiffusionSequenceFileHelper.h"
#include "btkOutlierCorrectionFilter.h"

#include "btkRBFInterpolateImageFunctionS2S.h"

#include "itkGradientImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

// Image and sequence definitions
typedef short                                                     PixelType;
typedef itk::Image< PixelType,3 >                                 TImage;

typedef btk::DiffusionSequence                                    TSequence;
typedef TSequence::IndexType                                      IndexType;
typedef TSequence::SizeType                                       SizeType;
typedef TSequence::RegionType                                     RegionType;

typedef itk::ImageRegionConstIteratorWithIndex< TSequence >       SequenceIterator;

// Filters definitions
typedef btk::OutlierCorrectionFilter<TImage>                      OutlierCorrectionType;

int main(int argc, char *argv[])
{
    try
    {

        //////////////////////////////////////////////////////////////////////////
        //
        // Parse program's arguments
        //

        // Define command line object for program
        TCLAP::CmdLine cmd("Outliers correction on DWIs.", ' ', "1.0");

        // Define arguments
        TCLAP::ValueArg< std::string >  inputSequenceFileNameArg("i", "sequence", "Input diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >  referenceSequenceFileNameArg("r", "reference", "reference diffusion sequence", false, "", "string", cmd);
        TCLAP::ValueArg< std::string >  outputSequenceFileNameArg("o", "output", "Output diffusion sequence", true, "", "string", cmd);
        TCLAP::ValueArg< std::string >  similaritiesFileNameArg("s", "similarities", "Similarities filename", false, "", "string", cmd);
        TCLAP::ValueArg< std::string >  mpeFileNameArg("f", "", "MPE filename ", false, "", "string", cmd);

        TCLAP::ValueArg< std::string >  methodArg("", "method", "method : SH (Spherical Harmonics decomposition), "
                                                  "WSH ( Weighted Spherical Harmonics decomposition) "
                                                  "or RBF (Radial Basis Function)", false, "SH", "string", cmd);

        TCLAP::ValueArg<double>         rGraArg("","rgra", "Only for RBF : gaussian r_gra parameter", false, 0.15, "double", cmd);

        TCLAP::ValueArg<double>          radiusArg("","radius", "Only for WSH: radius of neighbor search", false, 1.0, "double", cmd);


        TCLAP::SwitchArg verboseModeArg("v", "verbose", "Verbose mode", cmd, false);

        // Parse command line
        cmd.parse(argc, argv);

        // Get back arguments' values
        std::string inputSequenceFileName       = inputSequenceFileNameArg.getValue();
        std::string referenceSequenceFileName   = referenceSequenceFileNameArg.getValue();
        std::string outputSequenceFileName      = outputSequenceFileNameArg.getValue();
        std::string similaritiesFileName        = similaritiesFileNameArg.getValue();
        std::string mpeFileName                 = mpeFileNameArg.getValue();
        std::string method                      = methodArg.getValue();

        double rGra = rGraArg.getValue();
        double radius = radiusArg.getValue();

        bool verboseMode = verboseModeArg.getValue();

        //////////////////////////////////////////////////////////////////////////
        //
        // Read images
        //
        // Test if mandatory files exists
        if(!btk::FileHelper::FileExist(inputSequenceFileName))
        {
            btkException("Input sequence file are missing or path or name is wrong (name.nii.gz is mandatory) !");

        }
        TSequence::Pointer  inputSequence = btk::DiffusionSequenceHelper::ReadSequence(inputSequenceFileName);

        std::vector< unsigned short>           BValues       = inputSequence->GetBValues();
        std::vector< btk::GradientDirection >  GradientTable = inputSequence->GetGradientTable();

        RegionType region = inputSequence->GetLargestPossibleRegion();
        SizeType   size = region.GetSize();

        //////////////////////////////////////////////////////////////////////////
        // Correct outliers
        //////////////////////////////////////////////////////////////////////////
        OutlierCorrectionType::Pointer Outlier = OutlierCorrectionType::New();
        Outlier -> SetInputSequence(inputSequence.GetPointer());
        if(similaritiesFileNameArg.isSet())
        {
            Outlier -> WriteSimilarities(similaritiesFileName);
        }
        if(method.compare("RBF")==0)
        {
            throw("RBF interpolation return -nan value");
        }
        Outlier->SetMethod(method);  // method: SH, WSH, RBF
        Outlier -> SetRgra(rGra); // only used for RBF interoplation
        Outlier -> SetVerboseMod(verboseMode);
        Outlier->SetRadius(radius); // radius of the neighbor search: only used for WSH

        try
        {
            Outlier -> Update();
        }
        catch( itk::ExceptionObject & excp )
        {
            btkCerrMacro(excp);
        }

        TSequence::Pointer correctedSequence = Outlier->GetOutput();

        //////////////////////////////////////////////////////////////////////////
        // MPE calculation between estimation sequence and reference sequence
        //////////////////////////////////////////////////////////////////////////
        if(mpeFileNameArg.isSet() && !referenceSequenceFileName.empty() )
        {
            std::vector<std::vector<double> > VectorMPE;
            VectorMPE.resize(size[2]);

            btkCoutMacro("Calulate MPE:")
                    TSequence::Pointer  referenceSequence = btk::DiffusionSequenceHelper::ReadSequence(referenceSequenceFileName);

            for(unsigned int i=0; i< size[2]; i++)
            {
                double mean_MPE =0.0;
                VectorMPE[i].resize(size[3]+1);
                for(unsigned int j=0; j< size[3]; j++)
                {

                    RegionType  inputRegion       = inputSequence->GetLargestPossibleRegion();
                    SizeType    inputSize         = inputRegion.GetSize();

                    IndexType outlierIndex; //= IndexVector[i];
                    outlierIndex[0]=0;
                    outlierIndex[1]=0;
                    outlierIndex[2]=i;
                    outlierIndex[3]=j;

                    SizeType  outlierSize;
                    outlierSize[0] = inputSize[0];
                    outlierSize[1] = inputSize[1];
                    outlierSize[2] = 1;
                    outlierSize[3] = 1;

                    RegionType outlierRegion;
                    outlierRegion.SetIndex(outlierIndex);
                    outlierRegion.SetSize(outlierSize);

                    SequenceIterator referenceIt (referenceSequence,outlierRegion );
                    SequenceIterator correctedIt (correctedSequence,outlierRegion );

                    double MPE = 0.0; // Mean Percentage Error

                    for(referenceIt.GoToBegin(), correctedIt.GoToBegin(); !referenceIt.IsAtEnd() && !correctedIt.IsAtEnd(); ++referenceIt, ++correctedIt)
                    {
                        double value1 = static_cast< double >(referenceIt.Get());
                        double value2 = static_cast< double >(correctedIt.Get());

                        if(value2 <= std::numeric_limits<double>::epsilon())
                        {
                            MPE+=1.0;
                        }
                        else
                        {
                            double valAbs;
                            if(value1 == value2)
                            {
                                valAbs = 0.0;
                            }
                            else
                            {
                                valAbs = std::fabs(value1 - value2)/ value2;
                            }

                            if(valAbs <= 0.0 )
                            {
                                valAbs =0.0;
                            }
                            MPE += valAbs;
                        }
                    }
                    double nb_points = inputSize[0]*inputSize[1];
                    MPE*= (100.0/nb_points);

                    VectorMPE[i][j]=MPE;
                    if(j>0)
                    {
                        mean_MPE+=MPE;
                    }
                }
                mean_MPE/=size[3]-1;
                VectorMPE[i][size[3]]=mean_MPE;
            }

            //////////////////////////////////////////////////////////////////////////
            //
            // Write outputs
            //
            std::string delimiter = ";";

            std::cout<<"  -> Writing MPE in "<<mpeFileName<<"... ";

            std::ofstream fichierout(mpeFileName.c_str(), std::ios::trunc);
            if(fichierout)
            {
                fichierout<<"slice"<<delimiter;
                for (unsigned short i= 1; i < size[3]; i++)
                {
                    fichierout<<"dwi"<<i<<delimiter;
                }
                fichierout<<"mean"<<delimiter<<std::endl;

                for (unsigned short j=0; j < size[2]; j++)
                {
                    fichierout <<j<<delimiter;
                    for (unsigned short i=1; i < size[3]; i++)
                    {
                        fichierout <<VectorMPE[j][i]<<delimiter;
                    }
                    fichierout <<VectorMPE[j][size[3]]<<std::endl;
                }
                fichierout.close();
                btkCoutMacro(" Done.");
            }
            else
                std::cerr << "Unreadable file !" << std::endl;

            std::cout<<std::endl;
        }

        // Write modified sequence
        btk::ImageHelper<  btk::DiffusionSequence>::WriteImage(correctedSequence, outputSequenceFileName);

        // Write new B-values with removed zero entries
        std::string bvalFileName = btk::FileHelper::GetRadixOf(outputSequenceFileName) + ".bval";
        btk::DiffusionSequenceFileHelper::WriteBValues(BValues ,bvalFileName);

        // Write new B-vector with removed zero entries
        std::string bvecFileName = btk::FileHelper::GetRadixOf(outputSequenceFileName) + ".bvec";
        btk::DiffusionSequenceFileHelper::WriteGradientTable(GradientTable,bvecFileName);

    }
    catch(itk::ExceptionObject &exception)
    {

        std::cerr << "ITK error:" << std::endl;
        std::cerr << exception << std::endl;
        exit(EXIT_FAILURE);
    }
    catch(const char* e)
    {
        std::cerr << e << std::endl;
        exit(EXIT_FAILURE);
    }
}
