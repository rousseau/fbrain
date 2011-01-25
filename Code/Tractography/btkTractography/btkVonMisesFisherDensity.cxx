/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

23 april 2010
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


#include "btkVonMisesFisherDensity.h"


#ifndef NDEBUG
    #include "assert.h"
#endif // NDEBUG

// STL includes
#include "cmath"
#include "cstdlib"
#include "ctime"
#include "iostream"

// Local includes
#include "btkVector.h"


namespace btk
{

VonMisesFisherDensity::VonMisesFisherDensity(Real K) : SimDensity()
{
    #ifndef NDEBUG
        assert(K >= 0.);
    #endif // NDEBUG

    m_K  = K;
    m_c  = 2./m_K * std::sinh(m_K);
    m_C  = m_K / (4.*M_PI * std::sinh(m_K));
}

void VonMisesFisherDensity::setMeanDirection(Direction u)
{
    m_mean = u;
}

Real VonMisesFisherDensity::compute(Direction u)
{
    Vector v     =  u.toVector();
    Vector mean  = m_mean.toVector();
    Real vectors = mean.x()*v.x() + mean.y()*v.y() + mean.z()*v.z();

    return m_C * std::exp(m_K * vectors);
}

Real VonMisesFisherDensity::compute(Direction u, Real kappa)
{
    Vector v     =  u.toVector();
    Vector mean  = m_mean.toVector();
    Real vectors = mean.x()*v.x() + mean.y()*v.y() + mean.z()*v.z();

    Real C = kappa / (4.*M_PI * std::sinh(kappa));
//    return m_C * std::exp(m_K * vectors);
    return C * std::exp(kappa * vectors);
}

Real VonMisesFisherDensity::compute(Direction u, Direction mean, Real kappa)
{
    Vector v = u.toVector();
    Vector m = mean.toVector();
    Real vectors = m.x()*v.x() + m.y()*v.y() + m.z()*v.z();

    Real C = kappa / (4.*M_PI * std::sinh(kappa));
//    return m_C * std::exp(m_K * vectors);
    return C * std::exp(kappa * vectors);
}

Direction VonMisesFisherDensity::simulate()
{
    // Sample random scalar
    Real y = (Real)std::rand() / (Real)RAND_MAX;
    Real w = 1.0/m_K * std::log(std::exp(-m_K) + m_K * m_c * y);

    // Sample random angle (to get random unit vector)
    Real angle = ((Real)std::rand() / (Real)RAND_MAX) * 2.0 * M_PI;

    // Concatenate to obtain unit vector with vmf distribution
    Real cst = std::sqrt(1-w*w);
    Vector X(cst*std::cos(angle), cst*std::sin(angle), w);

    Real theta, phi;
    Vector mean = m_mean.toVector();

    if(!(mean.x() == 0 && mean.y() == 0 && mean.z() == 1)) // mean not (0,0,1)
    {
        // Rotations from (0,0,1) to mean
        Real cosTheta = std::cos(m_mean.theta());
        Real sinTheta = std::sin(m_mean.theta());
        Real cosPhi   = std::cos(m_mean.phi());
        Real sinPhi   = std::sin(m_mean.phi());

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

Direction VonMisesFisherDensity::simulate(Direction mean, Real kappa)
{
    // Sample random scalar
    Real y = (Real)std::rand() / (Real)RAND_MAX;
//    Real w = 1.0/m_K * std::log(std::exp(-m_K) + m_K * m_c * y);
    Real c = 2./kappa * std::sinh(kappa);
    Real w = 1.0/kappa * std::log(std::exp(-kappa) + kappa * c * y);

    // Sample random angle (to get random unit vector)
    Real angle = ((Real)std::rand() / (Real)RAND_MAX) * 2.0 * M_PI;

    // Concatenate to obtain unit vector with vmf distribution
    Real cst = std::sqrt(1-w*w);
    Vector X(cst*std::cos(angle), cst*std::sin(angle), w);

    Real theta, phi;
    Vector m = mean.toVector();

    if(!(m.x() == 0 && m.y() == 0 && m.z() == 1)) // mean not (0,0,1)
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

} // namespace btk

