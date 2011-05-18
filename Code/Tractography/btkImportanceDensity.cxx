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

ImportanceDensity::ImportanceDensity(Signal *signal, SHModel *model, Real angleThreshold)
{
    assert(model);

    m_signal         = signal;
    m_model          = model;
    m_angleThreshold = angleThreshold / 2.0;
    m_2PI            = 2.0 * M_PI;
    m_log4PI         = std::log(4.0 * M_PI);
}

Direction ImportanceDensity::computeMeanDirection(Point xk, Direction ukm1)
{
    std::vector<Direction> maxima = m_model->getMaxDirectionsAt(xk);

    if(maxima.size() > 1)
    {
        Vector vkm1 = ukm1.toVector();
        std::vector<Direction> mus;
        std::vector<Real> values;
        Real sum = 0;

        // we keep only maxima in the solid angle
        for(std::vector<Direction>::iterator it=maxima.begin(); it!=maxima.end(); it++)
        {
            Vector v = it->toVector();

            Real scalProd = v.x()*vkm1.x() + v.y()*vkm1.y() + v.z()*vkm1.z();
            Real normU    = v.x()*v.x() + v.y()*v.y() + v.z()*v.z();
            Real normV    = vkm1.x()*vkm1.x() + vkm1.y()*vkm1.y() + vkm1.z()*vkm1.z();
            Real alpha    = std::acos( scalProd / std::sqrt(normU * normV) );

            if(alpha <= m_angleThreshold)
            {
                // TODO : proportions en fonction de l'état précédent (angles) et plus amplitude
                mus.push_back(*it);
                values.push_back(m_model->odfAt(*it, xk));
                sum += values.back();
            }
        } // for each maxima

        if(mus.size() > 1)
        {
            // simulate x~U(0,1)
            Real x = ((Real)std::rand() / (Real)RAND_MAX) * sum;

            // compare to intervals and choose the mean direction
            bool found = false;
            unsigned int i = 0;
            Real interval = 0;

            while(!found && i<values.size())
            {
                interval += values[i];

                if(x < interval)
                    found = true;
                else
                    i++;
            }

            return mus[i];
        }
        else if(mus.size() == 1)
            return mus.front();
        else // mus.size() == 0
        {
//            return ukm1;
            Direction nullDir;
            nullDir.setNull();

            return nullDir;
        }
    }
    else if(maxima.size() == 1)
        return maxima.front();
    else // maxima.size() == 0
    {
//        return ukm1;
        Direction nullDir;
        nullDir.setNull();

        return nullDir;
    }
}

Direction ImportanceDensity::simulate(Direction mean, Real kappa)
{
    Vector X;

    do
    {
        // Sample random scalar
        Real y = (Real)std::rand() / (Real)RAND_MAX;
        Real w = 1.0/kappa * std::log(std::exp(-kappa) + kappa * 2.0/kappa * std::sinh(kappa) * y);

        // Sample random angle (to get random unit vector)
        Real angle = ((Real)std::rand() / (Real)RAND_MAX) * m_2PI;

        // Concatenate to obtain unit vector with vmf distribution with mean = (0,0,1)
        Real cst = std::sqrt(1-w*w);
        X = Vector(cst*std::cos(angle), cst*std::sin(angle), w);
    } while(X.toDirection().theta() > m_angleThreshold);

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

Real ImportanceDensity::compute(Direction uk, Direction mean, Real kappa)
{
    Vector vk   = uk.toVector();
    Vector vmean = mean.toVector();

    return std::log(kappa) - m_log4PI - std::log(std::sinh(kappa)) + kappa * (vk.x()*vmean.x() + vk.y()*vmean.y() + vk.z()*vmean.z());
}

Real ImportanceDensity::computeConcentration(Direction mu, Point xk)
{
    Matrix S = m_signal->signalAt(xk);
    Matrix M = m_model->signalAt(xk);

    Vector vmu = mu.toVector();
    std::vector<Direction> *directions = m_signal->getDirections();

    Real sumSquare = 0;
    Real min = MAX_REAL, max = MIN_REAL;
    unsigned int n = 0;

    for(unsigned int g=0; g<S.Rows(); g++)
    {
        Vector vg = directions->at(g).toVector();

        Real scalarProduct = vg.x()*vmu.x() + vg.y()*vmu.y() + vg.z()*vmu.z();
        Real normVG        = vg.x()*vg.x() + vg.y()*vg.y() + vg.z()*vg.z();
        Real normVMU       = vmu.x()*vmu.x() + vmu.y()*vmu.y() + vmu.z()*vmu.z();
        Real angle         = std::acos( scalarProduct / std::sqrt(normVG * normVMU) );

        if(angle <= m_angleThreshold)
        {
            Real diff = M(g,0) - S(g,0);
            sumSquare += diff*diff;

            if(min > M(g,0))
            {
                if(M(g,0) > S(g,0))
                    min = S(g,0);
                else
                    min = M(g,0);
            }

            if(max < M(g,0))
            {
                if(M(g,0) < S(g,0))
                    max = S(g,0);
                else
                    max = M(g,0);
            }

            n++;
        }
    }

    Real nrmse = std::sqrt(sumSquare / n) / (max - min);

    // TODO : calcul de concentration
    return 435.1455649013809 * std::exp(-25.74351897609373 * nrmse + 44.9437993453943 * nrmse * nrmse);
}

} // namespace btk

