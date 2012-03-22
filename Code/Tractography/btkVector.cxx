/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

14 april 2010
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

#include "btkVector.h"


// Local includes
#include "btkSphericalCoordinates.h"


namespace btk
{

Vector::Vector() : CartesianCoordinates()
{
    // ----
}

Vector::Vector(Real x, Real y, Real z) : CartesianCoordinates(x,y,z)
{
    // ----
}

bool Vector::isNull()
{
    return (this->x() == 0) && (this->y() == 0) && (this->z() == 0);
}

Direction Vector::toDirection()
{
    SphericalCoordinates s = this->toSphericalCoordinates();

    return Direction(s.theta(), s.phi());
}

Vector Vector::operator*(Real factor)
{
    return Vector(this->x()*factor, this->y()*factor, this->z()*factor);
}

Vector Vector::operator*(std::vector<Real> factors)
{
    return Vector(m_x*factors[0], m_y*factors[1], m_z*factors[2]);
}

Real Vector::angleWith(Vector v)
{
    Real normThis = std::sqrt( this->x()*this->x() + this->y()*this->y() + this->z()*this->z() );
    Real normV    = std::sqrt( v.x()*v.x() + v.y()*v.y() + v.z()*v.z() );
    Real dotProd  = this->x()*v.x() + this->y()*v.y() + this->z()*v.z();

    return std::acos( dotProd / (normThis * normV) );
}

void Vector::Normalize()
{
    Real norm = std::sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
    m_x /= norm;
    m_y /= norm;
    m_z /= norm;
}

} // namespace btk

