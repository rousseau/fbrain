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

#ifndef BTK_SUPERRESOLUTIONCOSTFUNCTIONITKWRAPPER_TXX
#define BTK_SUPERRESOLUTIONCOSTFUNCTIONITKWRAPPER_TXX

#include "btkSuperResolutionCostFunctionITKWrapper.h"


namespace btk
{
//-------------------------------------------------------------------------------------------------
template < class TImage >
SuperResolutionCostFunctionITKWrapper< TImage >::SuperResolutionCostFunctionITKWrapper()
{
    m_CostFunction = NULL;
    m_CostFunction = CostFunctionType::New();
}
//-------------------------------------------------------------------------------------------------
template < class TImage >
SuperResolutionCostFunctionITKWrapper< TImage >::~SuperResolutionCostFunctionITKWrapper()
{
    if(m_CostFunction !=NULL)
    {
        delete m_CostFunction;
        m_CostFunction = NULL;
    }
}

//-------------------------------------------------------------------------------------------------
template< typename TImage >
typename SuperResolutionCostFunctionITKWrapper< TImage >::MeasureType
SuperResolutionCostFunctionITKWrapper<TImage>::GetValue(const ParametersType &parameters) const
{
    MeasureType cost = this->m_CostFunction->f(parameters);

    return ( cost );
}

//-------------------------------------------------------------------------------------------------
template< typename TImage >
void SuperResolutionCostFunctionITKWrapper< TImage >::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
   this->m_CostFunction->GetGradient(parameters,derivative);
}

//-------------------------------------------------------------------------------------------------
template< typename TImage >
unsigned int SuperResolutionCostFunctionITKWrapper< TImage >::GetNumberOfParameters() const
{
    return m_NumberOfParameters;
}

//-------------------------------------------------------------------------------------------------
}

#endif // BTK_SUPERRESOLUTIONCOSTFUNCTIONITKWRAPPER_TXX
