/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 15/01/2013
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

#include "btkGreedyFeatureSelectionAlgorithm.h"


// VNL includes
#include "vnl/vnl_vector.h"

// STL includes
#include "algorithm"


namespace btk
{

GreedyFeatureSelectionAlgorithm::GreedyFeatureSelectionAlgorithm() : m_CurrentIteration(0), Superclass()
{
    // ----
}

//----------------------------------------------------------------------------------------

void GreedyFeatureSelectionAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void GreedyFeatureSelectionAlgorithm::Update()
{
    this->Initialize();

    // 1. Build an empty weighting vector as initialization
    unsigned int numberOfParameters = m_InputParameters->rows()/3;
    vnl_vector< short > w(numberOfParameters, 0);
    vnl_vector< double > e(numberOfParameters, 0.0);


    // 2. While parameters can be added without increasing the cost function
    bool                        stability = false;
    unsigned int numberOfActiveParameters = 0;

    double minCost = 1e6 * m_CostFunction->Evaluate();


    while(numberOfActiveParameters <= m_MaxNumberOfParameters && !stability)
    {
        std::stringstream message;

        // 2.1. Find a parameter to add (with the lesser cost)
        unsigned int index_add = numberOfParameters;

        message << "Initial cost = " << minCost << std::endl;

        std::vector< double > costsadd(numberOfParameters, minCost);

        double oldCost = minCost;

        #pragma omp parallel for default(shared) schedule(dynamic)
        for(unsigned int i = 0; i < numberOfParameters; i++)
        {
            if(w(i) != 1) // do not check already activated paramters
            {
                costsadd[i] = 1e6 * m_CostFunction->EvaluateActivation(i);
            }
        } // for each parameter

        // Search for minimal solution (if it exists)
        std::vector< double >::iterator itadd = std::min_element(costsadd.begin(), costsadd.end());

        if(*itadd < minCost)
        {
            minCost = *itadd;
            index_add = std::distance(costsadd.begin(), itadd);
        }


        // 2.2. Check the stability
        if(index_add == numberOfParameters)
        {
            stability = true;

            message << "\tCannot find a parameter to add (initial cost function: " << minCost << ")" << std::endl;
        }
        else // index_add != numberOfParameters
        {
            w(index_add)              = 1;
            e(index_add)              = oldCost - minCost;
            numberOfActiveParameters += 1;
            m_CostFunction->ActivateParameters(index_add);

            message << "\tFound one parameter to add (" << index_add+1 << ") with cost equal to " << minCost << std::endl;
            message << "\tUnexplained variance: " << minCost << " (" << 100.0 * (minCost / (minCost + e.sum())) << " %)" << std::endl;

//            if(100.0 * (minCost / (minCost + e.sum())) < 1.0)
//            {
//                stability = true;
//            }


            // 2.3 Remove parameters until the cost function cannot be decreased
            bool removeStability = false;

            while(!removeStability)
            {
                unsigned int index = numberOfParameters;

                std::vector< double > costssub(numberOfParameters, minCost);

                oldCost = minCost;

                #pragma omp parallel for default(shared) schedule(dynamic)
                for(unsigned int i = 0; i < numberOfParameters; i++)
                {
                    if(w(i) == 1 && i != index_add) // check only activated parameters not added recently
                    {
                        costssub[i] = 1e6 * m_CostFunction->EvaluateDesactivation(i);
                    }
                } // for each parameter

                // Search for minimal solution (if it exists)
                std::vector< double >::iterator itsub = std::min_element(costssub.begin(), costssub.end());

                if(*itsub < minCost)
                {
                    minCost = *itsub;
                    index = std::distance(costssub.begin(), itsub);
                }

                // 2.4. Check the stability
                if(index == numberOfParameters)
                {
                    removeStability = true;

                    message << "\tCannot find a parameter to remove (initial cost function: " << minCost << ")" << std::endl;
                }
                else // index != numberOfParameters
                {
                    w(index)                  = 0;
                    e(index)                  = oldCost - minCost;
                    numberOfActiveParameters -= 1;
                    m_CostFunction->DesactivateParameters(index);

                    message << "\tFound one parameter to remove (" << index+1 << ") with cost equal to " << minCost << std::endl;
                    message << "\tUnexplained variance: " << minCost << " (" << 100.0 * (minCost / (minCost + e.sum())) << " %)" << std::endl;

//                    if(100.0 * (minCost / (minCost + e.sum())) < 1.0)
//                    {
//                        removeStability = true;
//                    }
                }
            } // while parameters can be removed
        }

        m_CurrentMessage = message.str();
        this->InvokeEvent(itk::IterationEvent());
        m_CurrentIteration++;
    } // while a parameter can be added


    // Set energy vector as percent of explained variance
    double energySum = e.sum() + minCost;
    for(unsigned int i = 0; i < e.size(); i++)
    {
        e(i) /= energySum;
    }


    // 3. Set outputs
    m_WeightsVector = new vnl_vector< short >(w);
    m_EnergyVector = new vnl_vector< double >(e);
    m_UnexplainedVariance = minCost;
}

} // namespace btk
