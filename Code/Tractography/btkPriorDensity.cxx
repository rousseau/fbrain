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

#include "btkPriorDensity.h"


namespace btk
{

PriorDensity::PriorDensity() : m_Concentration(30.0)
{
    this->Initialize();
}

//----------------------------------------------------------------------------------------

PriorDensity::PriorDensity(double concentration) : m_Concentration(concentration)
{
    this->Initialize();
}

//----------------------------------------------------------------------------------------

void PriorDensity::Initialize()
{
    m_NormalizationCoefficient = std::log(m_Concentration / (m_2pi * (std::exp(m_Concentration) - std::exp(-m_Concentration))));
    m_InverseConcentration     = 1.0 / m_Concentration;
    m_SimulationTerm           = std::exp(-m_Concentration);
    m_SimulationCoefficient    = m_Concentration * ((2.0/m_Concentration) * (std::exp(m_Concentration) - std::exp(-m_Concentration))/2.0);
}

//----------------------------------------------------------------------------------------

double PriorDensity::Evaluate(GradientDirection vk, GradientDirection vkm1)
{
    return m_NormalizationCoefficient + (vkm1 * vk * m_Concentration);
}

//----------------------------------------------------------------------------------------

GradientDirection PriorDensity::Simulate(GradientDirection vkm1)
{
    // Simulate the vMF distribution with mean (0,0,1)
    double        y = static_cast< double >(std::rand()) / static_cast< double >(RAND_MAX);
    double        w = m_InverseConcentration * std::log( m_SimulationTerm + m_SimulationCoefficient * y);
    double    theta = ( static_cast< double >(std::rand()) / static_cast< double >(RAND_MAX)) * m_2pi;
    double constant = std::sqrt(1-w*w);

    GradientDirection x(constant*std::cos(theta), constant*std::sin(theta), w);

    // Rotate if necessary
    if(vkm1[0] != 0 || vkm1[1] != 0 || vkm1[2] != 1) // modal direction is not (0,0,1)
    {
        SphericalDirection mu = vkm1.GetSphericalDirection();

        // Rotations from (0,0,1) to current modal direction
        double cosTheta = std::cos(mu[0]);
        double sinTheta = std::sin(mu[0]);
        double cosPhi   = std::cos(mu[1]);
        double sinPhi   = std::sin(mu[1]);

        // Ry(theta)
        // | cos(theta)  0  sin(theta)  |   |x|
        // |     0       1       0      | . |y|
        // |-sin(theta)  0   cos(theta) |   |z|
        double Ryx = cosTheta * x[0] + sinTheta * x[2];
        double Ryy = x[1];
        double Ryz = -sinTheta * x[0] + cosTheta * x[2];

        // Rz(phi)
        // | cos(phi)  -sin(phi)  0 |   |x|
        // | sin(phi)   cos(phi)  0 | . |y|
        // |    0          0      1 |   |z|
        x[0] = cosPhi * Ryx - sinPhi * Ryy;
        x[1] = sinPhi * Ryx + cosPhi * Ryy;
        x[2] = Ryz;

        x.UpdateSphericalCoordinates();
    }

    return x;
}

} // namespace btk
