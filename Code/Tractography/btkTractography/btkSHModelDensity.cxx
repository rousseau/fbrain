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


#include "btkSHModelDensity.h"


// Debug
#ifndef NDEBUG
    #include "cassert"
#endif // NDEBUG


// STL includes
#include "cmath"
#include "cstdlib"
#include "ctime"
#include "iostream"


namespace btk
{

SHModelDensity::SHModelDensity(SHModel *model) : SimDensity()
{
    #ifndef NDEBUG
        assert(model);
    #endif // NDEBUG

    m_model = model;

    m_pi = 1.0 / (2.0 * M_PI * M_PI);
}

Real SHModelDensity::compute(Direction u)
{
/*
    return ( (m_model->odfAt(u, m_p) - m_min) / (m_max - m_min) )  / m_norm;
*/
    return this->compute(m_p, u);
}

Real SHModelDensity::compute(Point p, Direction u)
{
    return ( (m_model->odfAt(u, p) - m_min) / (m_max - m_min) )  / m_norm;
}

Direction SHModelDensity::simulate()
{

    Real u, theta, phi;
    Direction direction;

    do
    {
        // u ~ U(0,1)
        u = (Real)std::rand() / (Real)RAND_MAX;

        // (theta,phi) ~ pi
        theta = M_PI * ((Real)std::rand() / (Real)RAND_MAX);
        phi   = 2.0*M_PI * ((Real)rand() / (Real)RAND_MAX);

        direction = Direction(theta, phi);

    } while(u >= this->compute(direction) / (m_c * m_pi));

    return direction;

//    return this->simulate(m_p);
}

Direction SHModelDensity::simulate(Point p)
{
//////////////////////////////////////////////////
    // Get ODF at this position
    Matrix Psi = m_model->odfAt(p);

    // Search max, min
    Real max = Psi(0,0), min = Psi(0,0);
    for(unsigned int i=1; i<Psi.Rows(); i++)
    {
        if(max < Psi(i,0))
            max = Psi(i,0);

        if(min > Psi(i,0))
            min = Psi(i,0);
    } // for i

    // Search for normalisation factor
    Real norm = 0;
    for(unsigned int i=0; i<Psi.Rows(); i++)
        norm += (Psi(i,0) - min) / (max - min);

    // c is the max of ODF (for rejection sampling)
    Real c = (1.0 / norm) / m_pi;
////////////////////////////////////////////////////

    Real u, theta, phi;
    Direction direction;

    do
    {
        // u ~ U(0,1)
        u = (Real)std::rand() / (Real)RAND_MAX;

        // (theta,phi) ~ pi
        theta = M_PI * ((Real)std::rand() / (Real)RAND_MAX);
        phi   = 2.0*M_PI * ((Real)rand() / (Real)RAND_MAX);

        direction = Direction(theta, phi);

    } while(u >= this->compute(p, direction) / (c * m_pi));

    return direction;
}

void SHModelDensity::setPoint(Point p)
{
    m_p = p;

    // Get ODF at this position
    Matrix Psi = m_model->odfAt(p);

    // Search max, min
    m_max = Psi(0,0), m_min = Psi(0,0);
    for(unsigned int i=1; i<Psi.Rows(); i++)
    {
        if(m_max < Psi(i,0))
            m_max = Psi(i,0);

        if(m_min > Psi(i,0))
            m_min = Psi(i,0);
    } // for i

    // Search for normalisation factor
    m_norm = 0;
    for(unsigned int i=0; i<Psi.Rows(); i++)
        m_norm += (Psi(i,0) - m_min) / (m_max - m_min);

    // c is the max of ODF (for rejection sampling)
    m_c = (1.0 / m_norm) / m_pi;
}

} // namespace btk

