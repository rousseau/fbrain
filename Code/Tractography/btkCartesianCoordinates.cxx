/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

22 mars 2010
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

#include "btkCartesianCoordinates.h"


namespace btk
{

CartesianCoordinates::CartesianCoordinates()
{
    m_x = 0;
    m_y = 0;
    m_z = 0;
}

CartesianCoordinates::CartesianCoordinates(Real x, Real y, Real z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

Real CartesianCoordinates::x() const
{
    return m_x;
}

Real CartesianCoordinates::y() const
{
    return m_y;
}

Real CartesianCoordinates::z() const
{
    return m_z;
}

SphericalCoordinates CartesianCoordinates::toSphericalCoordinates()
{
    Real rau   = std::sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
    Real theta, phi;

    if(m_x == 0 && m_y == 0)    // on z-axis
    {
        theta = 0;
        phi = 0;
    }
    else
    {
        theta = std::acos(m_z / rau);
        phi   = std::acos(m_x / std::sqrt(m_x*m_x + m_y*m_y));

        if(m_y < 0)
            phi = 2.0 * M_PI - phi;
    }

    return SphericalCoordinates(theta, phi, rau);
}

CartesianCoordinates CartesianCoordinates::operator*(Real factor)
{
    return CartesianCoordinates(m_x*factor, m_y*factor, m_z*factor);
}

std::ostream &operator<<(std::ostream &os, const CartesianCoordinates& c)
{
    return os << "(" << c.m_x << ", " << c.m_y << ", " << c.m_z << ")";
}

} // namespace btk

