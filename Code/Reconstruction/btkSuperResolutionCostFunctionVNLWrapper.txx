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

#ifndef BTKSUPERRESOLUTIONCOSTFUNCTIONVNLWRAPPER_TXX
#define BTKSUPERRESOLUTIONCOSTFUNCTIONVNLWRAPPER_TXX

#include "btkSuperResolutionCostFunctionVNLWrapper.hxx"

namespace btk
{

//-------------------------------------------------------------------------------------------------
template < class TImage >
SuperResolutionCostFunctionVNLWrapper< TImage >::SuperResolutionCostFunctionVNLWrapper(unsigned int dim)
    :vnl_cost_function(dim)
{
    m_CostFunction  = NULL;
    m_CostFunction = btk::SuperResolutionCostFunction< TImage >::New();
    m_CostFunction->SetNumberOfParameters(dim);
}
//-------------------------------------------------------------------------------------------------
template < class TImage >
SuperResolutionCostFunctionVNLWrapper< TImage >::~SuperResolutionCostFunctionVNLWrapper()
{

}
//-------------------------------------------------------------------------------------------------
template < class TImage >
double
SuperResolutionCostFunctionVNLWrapper< TImage >::f(const vnl_vector< double >& _x)
{
    double value = this->m_CostFunction->GetValue(_x);

    return value;
}

//-------------------------------------------------------------------------------------------------
template < class TImage >
void
SuperResolutionCostFunctionVNLWrapper< TImage >::gradf(const vnl_vector< double >& _x, vnl_vector< double >& _g)
{
    this->m_CostFunction->GetGradient(_x,_g);
}

//-------------------------------------------------------------------------------------------------


}

#endif // BTKSUPERRESOLUTIONCOSTFUNCTIONVNLWRAPPER_TXX
