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


#include "btkImportanceDensity.h"


// STL includes
#include "cassert"
#include "cmath"


namespace btk
{

ImportanceDensity::ImportanceDensity(VonMisesFisherDensity d, SHModel *model, Real angleThreshold) : m_d(d)
{
    assert(model);


    m_model = model;

    m_angleThreshold = angleThreshold / 2.0;
}

Direction ImportanceDensity::computeMeanDirection(Point xk, Direction ukm1)
{
    std::vector<Direction> maxima = m_model->getMaxDirectionsAt(xk);

    if(maxima.size() > 1)
    {
        Vector v = ukm1.toVector();
        std::vector<Direction> mus;
        std::vector<Real> values;
        Real sum = 0;

        // we keep only maxima in the solid angle
        for(std::vector<Direction>::iterator it=maxima.begin(); it!=maxima.end(); it++)
        {
            Vector u = it->toVector();

            Real scalProd = u.x()*v.x() + u.y()*v.y() + u.z()*v.z();
            Real normU    = u.x()*u.x() + u.y()*u.y() + u.z()*u.z();
            Real normV    = v.x()*v.x() + v.y()*v.y() + v.z()*v.z();
            Real alpha    = std::acos( scalProd / std::sqrt(normU * normV) );

            if(alpha <= m_angleThreshold)
            {
                mus.push_back(*it);
                values.push_back(m_model->odfAt(*it, xk));
                sum += values.back();
            }
        } // for each maxima

        if(mus.size() > 1)
        {
            // simulate x~U(0,1)
            Real x = std::rand() / (Real)RAND_MAX;

            // compare to intervals and choose the mean direction
            bool found = false;
            unsigned int i = 0;
            Real interval = 0;

            while(!found && i<values.size())
            {
                interval += values[i]/sum;

                if(x < interval)
                    found = true;
                else
                    i++;
            }

            return mus[i];
        }
        else if(mus.size() == 1)
            return mus.front();
        else
            return ukm1;
    }
    else if(maxima.size() == 1)
        return maxima.front();
    else // maxima.size() == 0
        return ukm1;
}

Direction ImportanceDensity::simulate(Direction mean, Real kappa)
{
    return m_d.simulate(mean, kappa);
}

Real ImportanceDensity::compute(Direction vk, Direction mean, Real kappa)
{
    return m_d.compute(vk, mean, kappa);
}

Real ImportanceDensity::computeConcentration(Direction mu, Point xk)
{
    Matrix Psi = m_model->odfAt(xk);

    Real min = Psi(0,0), max = Psi(0,0), sum = 0;
    for(unsigned int i=0; i<Psi.Rows(); i++)
    {
        if(min > Psi(i,0))
            min = Psi(i,0);

        if(max < Psi(i,0))
            max = Psi(i,0);

        sum += Psi(i,0);
    }

    Real psi  = ((m_model->odfAt(mu, xk)-min)/(max-min));
    Real psi2 = psi*psi;

    Real sintheta  = std::sin(mu.theta());
    Real sin2theta = sintheta*sintheta;

    Real h = 0.1;
    Real psi_theta_theta =
        ( ((m_model->odfAt(Direction(mu.theta()+h,mu.phi()), xk)-min)/(max-min))/sum - 2.0*psi + ((m_model->odfAt(Direction(mu.theta()-h, mu.phi()), xk)-min)/(max-min))/sum )
        /
        (h*h);
    Real psi_phi_phi =
        ( ((m_model->odfAt(Direction(mu.theta(),mu.phi()+h), xk)-min)/(max-min))/sum - 2.0*psi + ((m_model->odfAt(Direction(mu.theta(), mu.phi()-h), xk)-min)/(max-min))/sum )
        /
        (h*h);

    Real e = psi - psi_theta_theta;
    Real G = psi2 * sin2theta;
    Real g = psi * sin2theta - psi_phi_phi;
    Real E = psi2;

    Real H = 0.5 * ( (e*G + g*E) / (E*G) );

    return std::exp(H*0.017);
}

} // namespace btk

