/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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

#include "btkSimulatedAnnealingOptimizer.h"

namespace btk
{


//-------------------------------------------------------------------------------------------------
SimulatedAnnealingOptimizer::SimulatedAnnealingOptimizer()
    :m_CurrentValue(0.0),m_Temperature(5000),m_Iteration(2000)
{
    Superclass::m_Stop = false;
}
//-------------------------------------------------------------------------------------------------
void SimulatedAnnealingOptimizer::StartOptimization()
{
    srand(time(NULL));
    if(this->m_CostFunction.IsNull() || this->m_Stop)
    {
        return;
    }

    //std::cout<<"Optimization..."<<std::endl;

    m_StopConditionDescription.str("");
    m_StopConditionDescription<<this->GetNameOfClass()<<": ";

    this->InvokeEvent(itk::StartEvent());
    m_Stop = false;

    unsigned int n = this->m_CostFunction->GetNumberOfParameters();
    m_RangeNeighboor.SetSize(n);

    m_RangeNeighboor[0] = 1.0;
    m_RangeNeighboor[1] = 1.0;
    m_RangeNeighboor[2] = 1.0;

//    m_RangeNeighboor[3] = 0.0;
//    m_RangeNeighboor[4] = 0.0;
//    m_RangeNeighboor[5] = 0.0;

//    m_RangeNeighboor[6] = 0.2;
//    m_RangeNeighboor[7] = 0.2;
//    m_RangeNeighboor[8] = 0.2;

    m_RangeNeighboor[3] = 1.0;
    m_RangeNeighboor[4] = 1.0;
    m_RangeNeighboor[5] = 1.0;


    ParametersType x, xn;
    x.SetSize(n);
    xn.SetSize(n);
    x = this->GetInitialPosition();
    xn =  x;
    this->SetCurrentPosition(x);
    double cost = this->m_CostFunction->GetValue(x);
    this->m_CurrentValue = cost;
    double costn;
    double temp = m_Temperature;
    unsigned int it = 0;



    while(it< m_Iteration && cost > 0.0001)
    {
        xn = GenerateRandomNeighboor(x);
        costn = this->m_CostFunction->GetValue(xn);
        temp = 0.7 *temp; // for decrease the temperature at each step
        double Dcost = costn-cost;
        double p = exp(-(Dcost/temp));

        if(costn < cost)
        {
            x = xn;
            this->SetCurrentPosition(x);
            cost = costn;
            this->m_CurrentValue = cost;
        }
        else
        {
            if(Random() < p)
            {
                x = xn;
                this->SetCurrentPosition(x);
                cost = costn;
                this->m_CurrentValue = cost;
            }
        }

        it++;
    }

    this->InvokeEvent(itk::EndEvent());

    //std::cout<<"Best cost : "<<this->m_CurrentValue<<" best parameters : "<<this->GetCurrentPosition()<<std::endl;





}
//-------------------------------------------------------------------------------------------------
SimulatedAnnealingOptimizer::ParametersType
SimulatedAnnealingOptimizer::GenerateRandomNeighboor(ParametersType _x)
{
    ParametersType xn;
    unsigned int n = this->m_CostFunction->GetNumberOfParameters();
    xn.SetSize(n);
    xn.Fill(0.0);

    for(unsigned int i =0; i<n; i++)
    {
        double p = Random();

        p = p * m_RangeNeighboor[i];

        xn[i] = (_x[i] - m_RangeNeighboor[i]/2.0) + p;


    }

    return xn;

}


//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
}
