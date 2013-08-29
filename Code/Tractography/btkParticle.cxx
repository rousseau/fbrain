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

#include "btkParticle.h"


namespace btk
{

Particle::Particle() : m_IsActive(true)
{
    // ----
}

//----------------------------------------------------------------------------------------

Particle::Particle(const Self &p)
{
    m_IsActive   = p.m_IsActive;
    m_Points     = p.m_Points;
    m_Directions = p.m_Directions;
    m_Weights    = p.m_Weights;
}

//----------------------------------------------------------------------------------------

Particle::Particle(PhysicalPoint p) : m_IsActive(true)
{
    m_Points.push_back(p);
}

//----------------------------------------------------------------------------------------

Particle::PhysicalPoint Particle::GetLastPoint() const
{
    return m_Points.back();
}

//----------------------------------------------------------------------------------------

GradientDirection Particle::GetLastDirection() const
{
    return m_Directions.back();
}

//----------------------------------------------------------------------------------------

double Particle::GetLastWeight() const
{
    assert(m_Weights.size() > 0);
    return m_Weights.back();
}

//----------------------------------------------------------------------------------------

bool Particle::IsActive() const
{
    return m_IsActive;
}

//----------------------------------------------------------------------------------------

void Particle::AddToPath(GradientDirection v, PhysicalPoint p, double w)
{
    m_Directions.push_back(v);
    m_Points.push_back(p);
    m_Weights.push_back(w);
}

//----------------------------------------------------------------------------------------

void Particle::Desactivate()
{
    m_IsActive = false;
}

//----------------------------------------------------------------------------------------

Particle::PhysicalPoint Particle::GetPointAtStep(unsigned int i) const
{
    assert(i < m_Points.size());
    return m_Points[i];
}

//----------------------------------------------------------------------------------------

double Particle::GetWeightAtStep(unsigned int i) const
{
    assert(i < m_Weights.size());
    return m_Weights[i];
}

//----------------------------------------------------------------------------------------

unsigned int Particle::GetPathLength() const
{
    return m_Points.size();
}

//----------------------------------------------------------------------------------------

GradientDirection Particle::GetVectorAtStep(unsigned int i) const
{
    assert(i < m_Directions.size());
    return m_Directions[i];
}

//----------------------------------------------------------------------------------------

void Particle::NormalizeLastWeightWith(double normalizationCoefficient)
{
    assert(m_Weights.size() > 0);
    m_Weights[m_Weights.size()-1] /= normalizationCoefficient;
}

//----------------------------------------------------------------------------------------

void Particle::Resample(PhysicalPoint p, double w)
{
    assert(m_Points.size() > 0);
    assert(m_Weights.size() > 0);
    m_Points[m_Points.size()-1]   = p;
    m_Weights[m_Weights.size()-1] = w;
}

//----------------------------------------------------------------------------------------

void Particle::AddLikelihood(double likelihood)
{
//    m_LikelihoodLog.push_back(std::log(likelihood));
    m_LikelihoodLog.push_back(likelihood);
}

//----------------------------------------------------------------------------------------

double Particle::GetLikelihoodAtStep(unsigned int i) const
{
    assert(i < m_LikelihoodLog.size());
    return m_LikelihoodLog[i];
}

} // namespace btk
