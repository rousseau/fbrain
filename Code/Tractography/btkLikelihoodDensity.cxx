/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 07/03/2013
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

#include "btkLikelihoodDensity.h"

// STL includes
#include "cfloat"


namespace btk
{

LikelihoodDensity::LikelihoodDensity() : m_Model(NULL), m_Signal(NULL)
{
    // ----
}

//----------------------------------------------------------------------------------------

LikelihoodDensity::LikelihoodDensity(DiffusionSignal::Pointer signal, DiffusionModel::Pointer model) : m_Model(model), m_Signal(signal)
{
    // ----
}

//----------------------------------------------------------------------------------------

double LikelihoodDensity::Evaluate(GradientDirection vk, PhysicalPoint xk, GradientDirection mu)
{
    // 1. Mean directon
    // we have it in parameters


    // 2. Compute alpha = arcos(vk.mu)
    double alpha = -std::acos(vk*mu);


    // 3. Compute rotation axis = (vk/\mu)
    GradientDirection axis = itk::CrossProduct(vk,mu);

    // 4. Build rotation matrix
    double sinAlpha = std::sin(alpha);
    double cosAlpha = std::cos(alpha);
    double axisx2   = axis[0]*axis[0];
    double axisy2   = axis[1]*axis[1];
    double axisz2   = axis[2]*axis[2];
    double axisxy   = axis[0]*axis[1];
    double axisxz   = axis[0]*axis[2];
    double axisyz   = axis[1]*axis[2];
    double R11 = axisx2 + (1.0-axisx2) * cosAlpha;
    double R12 = axisxy*(1.0-cosAlpha) - axis[2]*sinAlpha;
    double R13 = axisxz*(1.0-cosAlpha) + axis[1]*sinAlpha;
    double R22 = axisy2 + (1.0-axisy2)*cosAlpha;
    double R23 = axisyz*(1.0-cosAlpha) - axis[0]*sinAlpha;
    double R33 = axisz2 + (1.0-axisz2)*cosAlpha;


    // 5. Rotate all gradients coordinates
    std::vector< GradientDirection > g;
    std::vector< GradientDirection > directions = m_Signal->GetGradientTable();

    for(unsigned int i = 0; i < directions.size(); i++)
    {
        GradientDirection tmp1 = directions[i];
        GradientDirection tmp2(tmp1[0]*R11 + tmp1[1]*R12 + tmp1[2]*R13,
                               tmp1[0]*R12 + tmp1[1]*R22 + tmp1[2]*R23,
                               tmp1[0]*R13 + tmp1[1]*R23 + tmp1[2]*R33);
        g.push_back(tmp2);
    } // for dIt


    // 6. Compute densities
    double density = 0.0;

    DiffusionSignal::PixelType S = m_Signal->SignalAt(xk);
    std::vector< float >       M = m_Model->SignalAt(xk, g);

    std::vector< double > pseudoResidualsStdDeviation = m_Signal->GetPseudoResidualsStdDeviation();

    for(unsigned int i = 0; i < g.size(); i++)
    {
        double   mesuredSignal = S[i];
        double estimatedSignal = M[i];

        density += this->EvaluateNormalCenteredLogDensity(pseudoResidualsStdDeviation[i], mesuredSignal - estimatedSignal);
    } // for i

    return density - std::log(g.size());
}

//----------------------------------------------------------------------------------------

inline double LikelihoodDensity::EvaluateNormalCenteredLogDensity(double sigma, double x)
{
    double value = DBL_MIN;

    if(sigma > 0.0)
    {
        double fraction = x / sigma;
        value = -std::log(sigma) - m_logSqrt2Pi - 0.5 * fraction * fraction;
    }

    return value;
}

}
