/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 31/01/2014
  Author(s): François Rousseau
  
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

#include "btkPandoraBoxRegistrationFilters.h"

namespace btk
{
  void PandoraBoxRegistrationFilters::Register3DImages(itkFloatImagePointer & movingImage, itkFloatImagePointer & movingMaskImage, itkFloatImagePointer & referenceImage, itkFloatImagePointer & referenceMaskImage, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, vnl_vector< double > toleranceVector=vnl_vector<double>(), double tolerance=1e-3)
  {   
    PandoraBoxCostFunctionMSE myCostFunction;
    myCostFunction.SetReferenceImage(referenceImage);
    myCostFunction.SetMovingImage(movingImage);
    myCostFunction.SetMovingMask(movingMaskImage);
    myCostFunction.SetReferenceMask(referenceMaskImage);
    
    //myCostFunction.SetCenter();
    Register3DImages(myCostFunction, inputParameters, outputParameters, parameterRange, toleranceVector, tolerance);
  }

  void PandoraBoxRegistrationFilters::Register3DImages(PandoraBoxCostFunction & costFunction, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, vnl_vector< double > toleranceVector=vnl_vector<double>(), double tolerance=1e-3)
  {
    btk:Simplex toto;
    toto.SetConvergenceTolerance(tolerance);
    toto.SetConvergenceToleranceForEachParameter(toleranceVector);
    //std::cout<<"Size tolerance vector : "<<toleranceVector.size()<<std::endl;

    outputParameters = toto.minimize(inputParameters, parameterRange, costFunction);

    //It might appear that new estimate provide higher cost (nan issue maybe)
    if( costFunction(inputParameters) < costFunction(outputParameters) )
      outputParameters = inputParameters;

  }

  void PandoraBoxRegistrationFilters::GenerateRandomParameters(vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange)
  {
    //srand (time(NULL)); //Need to be done outside this function to be sure to get different random numbers!
    outputParameters = inputParameters;
    for(unsigned int j=0; j < inputParameters.size(); j++)
      outputParameters(j) = inputParameters(j) + 2 * ((rand() * 1.0 / RAND_MAX) - 0.5) * parameterRange(j);
  }

  void PandoraBoxRegistrationFilters::GenerateUniformlyDistributedParameters(vnl_vector< double > & inputParameters, std::vector< vnl_vector< double > > & outputParameters, vnl_vector< double > & parameterRange, int samplingRate)
  {
    outputParameters.clear();
    vnl_vector< double > tmpParameters(inputParameters.size(),0);

    //Loop over all possible initialization of rotation parameters
    for(unsigned int i=0; i < samplingRate; i++)
      for(unsigned int j=0; j < samplingRate; j++)
        for(unsigned int k=0; k < samplingRate; k++)
        {
          tmpParameters(0) = inputParameters(0) + (1.0 * i / (samplingRate-1) - 0.5) * 2 * parameterRange(0);
          tmpParameters(1) = inputParameters(1) + (1.0 * j / (samplingRate-1) - 0.5) * 2 * parameterRange(1);
          tmpParameters(2) = inputParameters(2) + (1.0 * k / (samplingRate-1) - 0.5) * 2 * parameterRange(2);

          for(unsigned int l=3; l < inputParameters.size(); l++)
            tmpParameters(l) = inputParameters(l);

          outputParameters.push_back(tmpParameters);
        }
  }

  void PandoraBoxRegistrationFilters::MultiStart3DRegistration(PandoraBoxCostFunction & costFunction, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, unsigned int numberOfPerturbations, vnl_vector< double > toleranceVector=vnl_vector<double>(), double tolerance=1e-3)
  {
    //Generer un ensemble de possibilite
    //Ensuite calculer le recalage pour chacun d'entre eux
    //Classer et selectionner

    vnl_vector< double > tmpParameters(inputParameters.size(),0);
    vnl_vector< double > currentParameters(inputParameters.size(),0);
    //Init random generator C++11 style
    //std::default_random_engine generator;
    //std::uniform_real_distribution<double> distribution(-1.0,1.0);
    //srand (time(NULL)); //Need to be done outside this function to be sure to get different random numbers!

    for(unsigned int i=0; i < numberOfPerturbations; i++)
    {
      std::cout<<"Multistart registration, iteration "<<i+1<<std::endl;
      //Perturb the input parameters according to min and max range for each parameter
      for(unsigned int j=0; j < inputParameters.size(); j++)
        //tmpParameters(j) = inputParameters(j) + distribution(generator) * parameterRange(j) ;
        tmpParameters(j) = inputParameters(j) + 2 * ((rand() * 1.0 / RAND_MAX) - 0.5) * parameterRange(j);

      //Then run registration process
      Register3DImages(costFunction, tmpParameters, currentParameters, parameterRange, toleranceVector, tolerance);

      //Do we get a better result?
      //If yes, set the output parameters as best estimates found
      if(i==0)
        outputParameters = currentParameters;
      else
      {
        double currentCost = costFunction(currentParameters);
        double bestCost    = costFunction(outputParameters);
        if(currentCost < bestCost)
          outputParameters = currentParameters;
      }
    }
  }

  void PandoraBoxRegistrationFilters::Coarse3DRegistration(PandoraBoxCostFunction & costFunction, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, unsigned int numberOfIncrements, vnl_vector< double > toleranceVector=vnl_vector<double>(), double tolerance=1e-3)
  {
    vnl_vector< double > tmpParameters(inputParameters.size(),0);
    vnl_vector< double > currentParameters(inputParameters.size(),0);

    //Loop over all possible initialization of rotation parameters
    for(unsigned int i=0; i < numberOfIncrements; i++)
      for(unsigned int j=0; j < numberOfIncrements; j++)
        for(unsigned int k=0; k < numberOfIncrements; k++)
        {
          tmpParameters(0) = inputParameters(0) + (1.0 * i / (numberOfIncrements-1) - 0.5) * 2 * parameterRange(0);
          tmpParameters(1) = inputParameters(1) + (1.0 * j / (numberOfIncrements-1) - 0.5) * 2 * parameterRange(1);
          tmpParameters(2) = inputParameters(2) + (1.0 * k / (numberOfIncrements-1) - 0.5) * 2 * parameterRange(2);

          for(unsigned int l=3; l < inputParameters.size(); l++)
            tmpParameters(l) = inputParameters(l);

          std::cout<<tmpParameters<<std::endl;

        }


  }

} // namespace btk
