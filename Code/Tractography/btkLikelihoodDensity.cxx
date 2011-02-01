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


#include "btkLikelihoodDensity.h"


// Debug
#ifndef NDEBUG
    #include "cassert"
#endif // NDEBUG

// STL includes
#include "cmath"


namespace btk
{

LikelihoodDensity::LikelihoodDensity(/*NormalDensity d, */Signal *signal, SHModel *model) /*: m_d(d)*/
{
    #ifndef NDEBUG
        assert(model);
        assert(signal);
    #endif // NDEBUG

    std::cout << "\tLikelihood density..." << std::flush;

    m_model  = model;
    m_signal = signal;

    m_directions = m_signal->getDirections();
    m_sigmas     = m_signal->getSigmas();

	m_logSqrt2PI = std::log(std::sqrt(2. * M_PI));

    std::cout << "done." << std::endl;
}

Real LikelihoodDensity::compute(Direction uk, Point xk, Direction mean)
{
    // 1. Mean directon
    // we have it in parameters


    // 2. Compute alpha = arcos(vk.mean)
    Vector vk   = uk.toVector();
    Vector vmean = mean.toVector();
    Real alpha   = std::acos(vk.x()*vmean.x() + vk.y()*vmean.y() + vk.z()*vmean.z());


    // 3. Compute rotation axis = (vk/\mean)
    Vector axis(vk.y()*vmean.z() - vk.z()*vmean.y(),
                vk.z()*vmean.x() - vk.x()*vmean.z(),
                vk.x()*vmean.y() - vk.y()*vmean.x());


    // 4. Build rotation matrix
    Real sinAlpha = std::sin(alpha);
    Real cosAlpha = std::cos(alpha);
    Real axisx2   = axis.x()*axis.x();
    Real axisy2   = axis.y()*axis.y();
    Real axisz2   = axis.z()*axis.z();
    Real axisxy   = axis.x()*axis.y();
    Real axisxz   = axis.x()*axis.z();
    Real axisyz   = axis.y()*axis.z();
    Real R11 = axisx2 + (1.-axisx2) * cosAlpha;
    Real R12 = axisxy*(1.-cosAlpha) - axis.z()*sinAlpha;
    Real R13 = axisxz*(1.-cosAlpha) + axis.y()*sinAlpha;
    Real R22 = axisy2 + (1.-axisy2)*cosAlpha;
    Real R23 = axisyz*(1.-cosAlpha) - axis.x()*sinAlpha;
    Real R33 = axisz2 + (1.-axisz2)*cosAlpha;


    // 5. Rotate all gradients coordinates
    std::vector<Direction> g;
    std::vector<Direction>::iterator dIt;

    for(dIt = m_directions->begin(); dIt != m_directions->end(); dIt++)
    {
        Vector tmp1 = dIt->toVector();
        Vector tmp2(tmp1.x()*R11 + tmp1.y()*R12 + tmp1.z()*R13,
                    tmp1.x()*R12 + tmp1.y()*R22 + tmp1.z()*R23,
                    tmp1.x()*R13 + tmp1.y()*R23 + tmp1.z()*R33);

        SphericalCoordinates s = tmp2.toSphericalCoordinates();
        Direction d(s.theta(), s.phi());
        g.push_back(d);
    } // for dIt


    // 6. Compute densities
    Real density = /*0*/1;
    Matrix S     = m_signal->signalAt(xk);
/*
    for(unsigned int i=0; i<g.size(); i++)
    {
        Real mesuredSignal   = S(i,0);
        Real estimatedSignal = m_model->signalAt(g[i], xk);

//		density += std::log(this->computeNormalDensity(m_sigmas->at(i), mesuredSignal - estimatedSignal));
        density += this->computeNormalDensity(m_sigmas->at(i), mesuredSignal - estimatedSignal);
    } // for i
*/

    return density;
}

Real LikelihoodDensity::computeNormalDensity(Real sigma, Real x)
{
    if(sigma != 0)
    {
//        Real coefficient = 1./(sigma * m_sqrt2PI);
//        Real fraction    = x / sigma;
//        Real exponent    = -0.5 * fraction * fraction;
//
//        return coefficient * std::exp(exponent);

        Real fraction = x / sigma;

        return ( -std::log(sigma) - m_logSqrt2PI - 0.5 * fraction * fraction );
    }
    else
        return 1;
}

} // namespace btk

