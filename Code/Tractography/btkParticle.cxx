/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

31 march 2010
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

#include "btkParticle.h"


// STL includes
#include "cassert"
#include "cmath"

// ITK includes
#include "itkPoint.h"


namespace btk
{

Particle::Particle(Point p)
        : m_points(1,p)
{
    m_outside = false;
    m_active  = true;
}

void Particle::setWeight(Real w)
{
    assert(w >= 0);
    assert(w <= 1);

    m_weight.push_back(w);
}

char Particle::addToPath(Vector v, Mask::Pointer mask)
{
    Point p = this->lastPoint()+v;

    char isIn = this->IsInside(p, mask);

    if(isIn == 0 || isIn == 2)
    {
        m_points.push_back(p);
        m_vectors.push_back(v);
        m_outside = true;
    }
    else if(isIn == 1)
    {
        m_points.push_back(p);
        m_vectors.push_back(v);
    }

    return isIn;
}

Real Particle::weight() const
{
    return m_weight.back();
}

void Particle::SetLastWeight(Real w)
{
    assert(w >= 0);
    assert(w <= 1);

    m_weight[m_weight.size()-1] = w;
}

Point Particle::lastPoint() const
{
    return m_points.back();
}

void Particle::SetLastPoint(Point p)
{
    m_points[m_points.size()-1] = p;
}

Vector Particle::lastVector() const
{
    return m_vectors.back();
}

void Particle::SetLastVector(Vector v)
{
    m_vectors[m_vectors.size()-1] = v;
}

Point Particle::getPoint(unsigned int i) const
{
    assert(i>=0);
    assert(i<m_points.size());

    return m_points[i];
}

Vector Particle::getVector(unsigned int i) const
{
    assert(i>=0);
    assert(i<m_vectors.size());

    return m_vectors[i];
}

Real Particle::getWeight(unsigned int i) const
{
    assert(i>=0);
    assert(i<m_weight.size());

    return m_weight[i];
}

unsigned int Particle::length() const
{
    assert(m_points.size() == m_vectors.size()+1);
    assert(m_vectors.size() == m_weight.size());

    return m_vectors.size();
}

bool Particle::isOutside()
{
    return m_outside;
}

bool Particle::isActive()
{
    return m_active;
}

void Particle::setActive(bool active)
{
    m_active = active;
}

std::ostream &operator<<(std::ostream &os, Particle p)
{
    os << p.m_weight.back() << " " << p.m_points.size() << " ";

    std::vector<Point>::iterator it;
    for(it = p.m_points.begin(); it != p.m_points.end(); it++)
        os << (*it) << " ";

    os << std::endl;

    return os;
}

char Particle::IsInside(Point p, Mask::Pointer mask)
{
    Mask::IndexType index;
    index[0] = std::floor(p.x() + 0.5);
    index[1] = std::floor(p.y() + 0.5);
    index[2] = std::floor(p.z() + 0.5);

    MaskRegion region = mask->GetLargestPossibleRegion();

    char isIn = 0;

    if(index[0] >= 0 && index[0] < (unsigned int)region.GetSize(0) &&
       index[1] >= 0 && index[1] < (unsigned int)region.GetSize(1) &&
       index[2] >= 0 && index[2] < (unsigned int)region.GetSize(2))
        isIn = mask->GetPixel(index);

    return isIn;
}

} // namespace btk

