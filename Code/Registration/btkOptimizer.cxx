/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 13/12/2012
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
#ifndef _btkOptimizer_cxx
#define _btkOptimizer_cxx

#include "btkOptimizer.h"

namespace btk
{
//---------------------------------------------------------------------------
Optimizer::Optimizer()
{
    m_Stop = false;
    m_StopConditionDescription<< this->GetNameOfClass()<< ":";

    m_CurrentValue = 0;
    m_CurrentCost = 0;

    m_MaximumNumberOfIterations = 100;



}

//---------------------------------------------------------------------------
Optimizer::~Optimizer()
{

}
//---------------------------------------------------------------------------
void Optimizer::StartOptimization(void)
{
    if(m_CostFunction.IsNull())
    {
        return;
    }

    m_StopConditionDescription.str("");
    m_StopConditionDescription<<this->GetNameOfClass()<<": ";

    this->InvokeEvent(itk::StartEvent());
    m_Stop = false;

    // TODO: here do the optimization


    //double val = this->m_CostFunction->GetValue(param);
    //this->m_InitialPosition;
    //this->m_CurrentPosition;
    //this->m_CurrentCost

}
//---------------------------------------------------------------------------
void Optimizer::SetCostFunction(SingleValuedCostFunction *costFunction)
{
    this->m_CostFunction = costFunction;
}
//---------------------------------------------------------------------------
const std::string Optimizer::GetStopConditionDescription() const
{
    return m_StopConditionDescription.str();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
}
#endif
