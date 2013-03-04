/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 01/03/2013
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

#include "btkNormalProbabilityDensity.h"

// STD includes
#include "cmath"

namespace btk
{

NormalProbabilityDensity::NormalProbabilityDensity() : Superclass(), m_Mu(0.0), m_Sigma(1.0)
{
    this->Initialize();
}

//----------------------------------------------------------------------------------------

NormalProbabilityDensity::NormalProbabilityDensity(double mu, double sigma) : Superclass(), m_Mu(mu), m_Sigma(sigma)
{
    this->Initialize();
}

//----------------------------------------------------------------------------------------

void NormalProbabilityDensity::Initialize()
{
    // Precompute some constants
    m_normalizationConstant = 1.0 / (m_Sigma * m_sqrt2pi);
    m_Sigma2                = m_Sigma * m_Sigma;
}

//----------------------------------------------------------------------------------------

double NormalProbabilityDensity::Evaluate(double x)
{
    double deviationToMean = x - m_Mu;
    return m_normalizationConstant * std::exp( -0.5 * (deviationToMean*deviationToMean) / m_Sigma2 );
}

//----------------------------------------------------------------------------------------

double NormalProbabilityDensity::Simulate()
{
    // Box-Muller algorithm
    double u = static_cast< double >(rand()) / static_cast< double >(RAND_MAX);
    double v = static_cast< double >(rand()) / static_cast< double >(RAND_MAX);

    return m_Sigma * std::sqrt(-2.0 * std::log(u)) * std::cos(m_2pi * v) + m_Mu;
}

} // namespace btk
