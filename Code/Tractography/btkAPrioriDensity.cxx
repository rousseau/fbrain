/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

15 june 2010
< pontabry at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#include "btkAPrioriDensity.h"


// STL includes
#include "cmath"


namespace btk
{

APrioriDensity::APrioriDensity(Real K)
{
    m_K = K;
    m_C = std::log(K) - std::log(4.0 * M_PI) - std::log(std::sinh(K));

    m_1OnK = 1.0 / m_K;
    m_c    = 2.0 / m_K * std::sinh(m_K);
    m_Kc   = m_K*m_c;
    m_emK  = std::exp(-m_K);
    m_2PI  = 2.0 * M_PI;
}

Real APrioriDensity::compute(Direction uk, Direction ukm1)
{
    Vector vk   = uk.toVector();
    Vector vkm1 = ukm1.toVector();

    return m_C + m_K * (vk.x()*vkm1.x() + vk.y()*vkm1.y() + vk.z()*vkm1.z());
}

Direction APrioriDensity::simulate(Direction mean)
{
    // Sample random scalar
    Real y = (Real)std::rand() / (Real)RAND_MAX;
    Real w = m_1OnK * std::log(m_emK + m_Kc * y);

    // Sample random angle (to get random unit vector)
    Real angle = ((Real)std::rand() / (Real)RAND_MAX) * m_2PI;

    // Concatenate to obtain unit vector with vmf distribution with mean = (0,0,1)
    Real cst = std::sqrt(1-w*w);
    Vector X(cst*std::cos(angle), cst*std::sin(angle), w);

    Real theta, phi;
    Vector vmean = mean.toVector();

    // Rotate if necessary
    if(!(vmean.x() == 0 && vmean.y() == 0 && vmean.z() == 1)) // mean not (0,0,1)
    {
        // Rotations from (0,0,1) to mean
        Real cosTheta = std::cos(mean.theta());
        Real sinTheta = std::sin(mean.theta());
        Real cosPhi   = std::cos(mean.phi());
        Real sinPhi   = std::sin(mean.phi());

        // Ry(theta)
        // | cos(theta)  0  sin(theta)  |   |x|
        // |     0       1       0      | . |y|
        // |-sin(theta)  0   cos(theta) |   |z|
        Vector tmp(
                cosTheta * X.x() + sinTheta * X.z(),
                X.y(),
                -sinTheta * X.x() + cosTheta * X.z()
        );

        // Rz(phi)
        // | cos(phi)  -sin(phi)  0 |   |x|
        // | sin(phi)   cos(phi)  0 | . |y|
        // |    0          0      1 |   |z|
        Vector tmp2(
                cosPhi * tmp.x() - sinPhi * tmp.y(),
                sinPhi * tmp.x() + cosPhi * tmp.y(),
                tmp.z()
        );

        theta = tmp2.toSphericalCoordinates().theta();
        phi   = tmp2.toSphericalCoordinates().phi();
    }
    else // vkm1 is (0,0,1)
    {
        theta = X.toSphericalCoordinates().theta();
        phi   = X.toSphericalCoordinates().phi();
    }

    return Direction(theta,phi);
}

Real APrioriDensity::getConcentration() const
{
    return m_K;
}

} // namespace btk

