/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/07/2013
  Author(s): Aïcha Bentaieb (abentaieb@unistra.fr)

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



/* ITK includes */
#include "itkVector.h"
#include "itkListSample.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"
#include "itkNormalVariateGenerator.h"
#include "itkGaussianDistribution.h"

/* OTHERS */
#include <tclap/CmdLine.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>



//---------------------------------------------------------------------------------------------------
// function computing the standard deviation of a vector array
double GetStDev(std::vector<double> array,int number_of_items, double mean)
{
    float numerator = 0;
    float denominator = number_of_items;

    for (int i = 1; i <= number_of_items; i++)
    {
        numerator = numerator + std::pow((array[i] - mean), 2);
    }

    double standard_deviation = sqrt (numerator/denominator);
    return standard_deviation;
}

//--------------------------------------------------------------------------------------------------
// function computing the mean of a data array
double GetMean(std::vector<double> array, int iSize) {
    double sum = array[0];
    for (int i = 1; i < iSize; ++i) {
        sum += array[i];
    }
    return sum/iSize;
}
//--------------------------------------------------------------------------------------------------
// Likelihood of a sequence of integers
double PXS( std::vector<std::vector<double> >patientMeasurement,std::vector<std::vector<double> > controlMeasurement, std::vector<std::vector<double> > distributionParamsPatients,std::vector<std::vector<double> > distributionParamsControl,std::vector<int> sequence )
{
    int J = patientMeasurement.size();
    int N = sequence.size();

    /* distribution parameters of each event considering whether it appeared (patients) or not (controls) */
    double meanE = 0.0;
    double stdevE = 0.0;
    double meanNE = 0.0;
    double stdevNE = 0.0;
    std::vector<double> likelihood; // likelihood per order
    std::vector<double> patientLikelihood; // likelihood of the sequence for each patient
    double P = 0.0; // total likelihood

    /* likelihood computation for the sequence */
    for (int j=0; j<J; j++)
    {
        for(int order = 0; order<N; order++ )
        {
            double PEk = 1.0; // likelihood of the occurance of the event E knowing event k
            double PNEk = 1.0; // likelihood of the non occurance of event E knowing we are at the order k
            /* Events that occured */
            for( int e=0; e<=order; e++)
            {
                int event = sequence[e];
                meanE = distributionParamsPatients[0][event];
                stdevE = distributionParamsPatients[1][event];
                itk::Statistics::GaussianDistribution::Pointer gaussian = itk::Statistics::GaussianDistribution::New();
                gaussian->SetMean(meanE);
                gaussian->SetVariance(stdevE);
                PEk = PEk * (gaussian->EvaluatePDF(patientMeasurement[j][event]));
            }
            /* Events that didnt occur */
            if (order != N-1)
            {
                for( int Ne = order+1; Ne<=N-1; Ne++)
               {
                    int Nevent = sequence[Ne];
                    meanNE = distributionParamsControl[0][Nevent];
                    stdevNE = distributionParamsControl[1][Nevent];
                    itk::Statistics::GaussianDistribution::Pointer gaussianC = itk::Statistics::GaussianDistribution::New();
                    gaussianC->SetMean(meanNE);
                    gaussianC->SetVariance(stdevNE);
                    PNEk = PNEk * (gaussianC->EvaluatePDF(controlMeasurement[j][Nevent]));
               }
            }

            likelihood.push_back( (PEk * PNEk)/N );
        }

        double sum = 0.0;
        for(int i=0; i<5;i++)
        {
            sum += likelihood[i];
        }
        patientLikelihood.push_back( sum );
        likelihood.clear();
    }

    for (int i=0; i<patientLikelihood.size(); i++)
    {
        patientLikelihood[i] = std::log(patientLikelihood[i]);
        P += patientLikelihood[i];
    }


    return P;


}

//--------------------------------------------------------------------------------------------------
void swap(int x, int y)
{
    int temp = x;
    x = y;
    y = temp;

    return;
}
//--------------------------------------------------------------------------------------------------
/* Factorial function */
int fact(int x)
{
    int f = 1;
    if (x == 0)
        return 1;

    for (int i=x; i>=1; i--)
    {
        f *= i ;
    }
    return f;
}

//--------------------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    try
    {
        /* ================================================================================================================
         * Define input files : txt files corresponding to measurement
         * ================================================================================================================
         */

        TCLAP::CmdLine cmd("Sequence prediction", ' ', "1.0", true);
        TCLAP::MultiArg<std::string> patientEventMeasurementArg("p", "Patient_measurement_File","Patient measurement data for each event (possible multi inputs) ", true, "string");
        cmd.add(patientEventMeasurementArg);
        TCLAP::MultiArg<std::string> controlEventMeasurementArg("c", "Controle_measurement_File","Control measurement data for each event (possible multi inputs) ", true, "string");
        cmd.add(controlEventMeasurementArg);
        TCLAP::ValueArg<std::string> outputFile("o", "output", "Txt file with the probabilities of each sequence of event",true," ","string" );
        cmd.add(outputFile);

        cmd.parse(argc,argv);

        // Store input data files in vector
        std::vector<std::string> patientsMeasurementFiles = patientEventMeasurementArg.getValue();
        std::vector<std::string> controlsMeasurementFiles = controlEventMeasurementArg.getValue();
        std::string output = outputFile.getValue();

        std::cout<< "Number of events considered: "<< patientsMeasurementFiles.size()<<"\n";
        for (unsigned int i=0; i<patientsMeasurementFiles.size(); i++)
        {
            std::cout<<"Patient measurement txt file:  "<<i<<" : "<<patientsMeasurementFiles[i]<<"\n";
        }
        std::cout<< "Number of control files: "<<controlsMeasurementFiles.size()<<"\n";
        for(unsigned int i=0; i<controlsMeasurementFiles.size(); i++)
        {
            std::cout<<"Control measurement txt file:  "<<i<<" : "<<controlsMeasurementFiles[i]<<"\n";
        }
        int NumberOfEvents = patientsMeasurementFiles.size();

        // find the number of measurement per event = Number of lines in each file
        std::string file = patientsMeasurementFiles[0];
        std::ifstream inFile(file.c_str());
        int NumberOfMeasurementsPerEvent = std::count(std::istreambuf_iterator<char>(inFile),std::istreambuf_iterator<char>(), '\n');

        std::cout<<"N° of measurements J = "<<NumberOfMeasurementsPerEvent<<std::endl;
        std::cout<<"N°of events N = "<<NumberOfEvents<<std::endl;

        /* =======================================================================================================================
         * Get the measurement for each patient and each control in a vector
         * =======================================================================================================================
         */
        std::string filename;
        std::vector<std::string> measure;
        std::string line;
        std::string curv;

        std::vector<std::vector<double> > M(NumberOfMeasurementsPerEvent, std::vector<double>(NumberOfEvents));// matrix of patient measurement per event
        std::vector<std::vector<double> > C(NumberOfMeasurementsPerEvent, std::vector<double>(NumberOfEvents));// matrix of control measurement per event

        // patients
        for(unsigned int i=0; i<patientsMeasurementFiles.size(); i++)
        {
            filename = patientsMeasurementFiles[i];
            std::ifstream file(filename.c_str());
            if(file)
            {
                while ( std::getline(file, line) )
                {
                    std::istringstream iss(line);

                    iss >> curv >> curv >> curv;
                    measure.push_back(curv);
                }
                file.close();
            }
            else
                std::cerr<< "Impossible!!!!"<<std::endl;

            for(int j=0; j<NumberOfMeasurementsPerEvent;j++)
            {
                std::istringstream buffer(measure[j]);
                buffer >> M[j][i] ;
            }
            measure.clear();
        }
        // control
        for(unsigned int i=0; i<controlsMeasurementFiles.size(); i++)
        {
            filename = controlsMeasurementFiles[i];
            std::ifstream file(filename.c_str());
            if(file)
            {
                while ( std::getline(file, line) )
                {
                    std::istringstream iss(line);
                    iss >> curv >> curv >> curv;
                    measure.push_back(curv);
                }
                file.close();
            }
            else
                std::cerr<< "Impossible!!!!"<<std::endl;

            for(int j=0; j<measure.size();j++)
            {
                std::istringstream buf(measure[j]);
                buf >> C[j][i];
            }
            measure.clear();
        }

        for (int i=0; i<380; i++)
            for(int j=0; j<5; j++)
            {  std::cout<<C[i][j]<<'\t';
                 std::cout<<std::endl;}
        /* =======================================================================================================================
         * Find the distribution of each event and each dataset (patients and controls)
         * =======================================================================================================================
         */

        std::vector<double> Event_i;
        std::vector<double> NEvent_i;
        std::vector<std::vector<double> > distributionParametersPatients(2, std::vector<double>(NumberOfEvents) ); // for each event, store the mean and the stdev
        std::vector<std::vector<double> > distributionParametersControls(2, std::vector<double>(NumberOfEvents) );
        double mean = 0.0;
        double stdev = 0.0;
        double meanNE = 0.0;
        double stdevNE = 0.0;
        for(int n=0; n<NumberOfEvents;n++)
        {
            for(int j=0; j<NumberOfMeasurementsPerEvent; j++)
            {
                Event_i.push_back(M[j][n]);
                NEvent_i.push_back(C[j][n]);
            }
            /* subjects */
            mean = GetMean(Event_i,Event_i.size());
            stdev = GetStDev(Event_i, Event_i.size(), mean);
            distributionParametersPatients[0][n] = mean;
            distributionParametersPatients[1][n] = stdev;
            /* controls*/
            meanNE = GetMean(Event_i,Event_i.size());
            stdevNE = GetStDev(Event_i, Event_i.size(), meanNE);
            distributionParametersControls[0][n] = meanNE;
            distributionParametersControls[1][n] = stdevNE;

        }

        /* =======================================================================================================================
         * Fonteijn algorithm:
         * =======================================================================================================================
         */

        // 1. find all the possible event combinations (S)

        int combinations = fact(NumberOfEvents);
        std::cout << "The n! possible permutations with"<< NumberOfEvents<<" elements: " << combinations<<std::endl;

        std::vector<std::vector<int> > sequences( combinations, std::vector<int>(NumberOfEvents,0) );
        std::vector<int> row = sequences[0];

        for(int i = 1; i <= NumberOfEvents ; i++)
        {
            row[i-1] = i;

        }


        for (int l=0; l<sequences.size(); l++)
        {
            std::next_permutation(row.begin(), row.end());
            for(int i=0; i<row.size(); i++)
            {
                std::cout<<row[i]<<'\t';
                sequences[l][i]=row[i];
            }std::cout<<std::endl;
        }

        // 2. Compute the likelihood of each sequence (S)

        double likelihoodPerSequence[combinations];
        double p = 0.0;
        for(int i=0; i<sequences.size();i++)
        {
            p = PXS( M, C, distributionParametersPatients, distributionParametersControls,sequences[i]);
            likelihoodPerSequence[i] = p;

        }

        // 3. Output the probabilities per sequence in a txt file
        std::ofstream Ofile(output.c_str());
        if(Ofile)
        {
            for(int i=0; i<sequences.size();i++)
            {
                Ofile << i << "," << likelihoodPerSequence[i] << std::endl;
            }
            Ofile.close();
        }
        else
            std::cerr << "ERROR, impossible to open file !" << std::endl;

    }catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }


  return EXIT_SUCCESS;
}
