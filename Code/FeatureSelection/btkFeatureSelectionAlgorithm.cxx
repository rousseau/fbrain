/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 07/01/2012
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

#include "btkFeatureSelectionAlgorithm.h"


namespace btk
{

FeatureSelectionAlgorithm::FeatureSelectionAlgorithm() : m_InputParameters(NULL), m_WeightsVector(NULL), m_EnergyVector(NULL), m_MaxNumberOfParameters(10), m_CurrentIteration(0), m_CurrentMessage("")
{
    // ----
}

//----------------------------------------------------------------------------------------

FeatureSelectionAlgorithm::~FeatureSelectionAlgorithm()
{
    if(m_WeightsVector != NULL)
    {
        delete m_WeightsVector;
        m_WeightsVector = NULL;
    }

    if(m_EnergyVector != NULL)
    {
        delete m_EnergyVector;
        m_EnergyVector = NULL;
    }
}

//----------------------------------------------------------------------------------------

void FeatureSelectionAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void FeatureSelectionAlgorithm::Initialize()
{
    if(m_InputParameters != NULL)
    {
        m_CostFunction->SetInputParameters(m_InputParameters);
        m_CostFunction->Initialize();
    }
    else
    {
        btkException("No input parameters are given !");
    }
}

} // namespace btk
