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


#include "btkDTFPLikelihoodDensity.h"


// STL includes
#include "cassert"
#include "cmath"

// ITK includes
#include "itkDiffusionTensor3D.h"
#include "itkFixedArray.h"
#include "itkMatrix.h"

// default 0.25
#define Tau 0.25


namespace btk
{

DTFPLikelihoodDensity::DTFPLikelihoodDensity(DTFPSignal *signal)
{
    assert(signal);


    m_signal = signal;

    m_directions = m_signal->getDirections();
    m_sigmas     = m_signal->getSigmas();

	m_logSqrt2PI = std::log(std::sqrt(2. * M_PI));
}

Real DTFPLikelihoodDensity::compute(Direction uk, Point xk, Direction mean)
{
    Real density = 0;
    itk::DiffusionTensor3D<Real> tensor = m_signal->DiffusionTensorAt(xk);

    itk::FixedArray<Real,3> eigenValues;
    itk::Matrix<Real,3,3> eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues,eigenVectors);
//    Pr(eigenValues[0]); Pr(eigenValues[1]); Pr(eigenValues[2]);
//    Pr(eigenVectors(0,0)); Pr(eigenVectors(0,1)); Pr(eigenVectors(0,2));
//    Pr(eigenVectors(1,0)); Pr(eigenVectors(1,1)); Pr(eigenVectors(1,2));
//    Pr(eigenVectors(2,0)); Pr(eigenVectors(2,1)); Pr(eigenVectors(2,2));

    Real cl = (eigenValues[2] - eigenValues[1]) / std::sqrt(eigenValues[2]*eigenValues[2] + eigenValues[1]*eigenValues[1] + eigenValues[0]*eigenValues[0]);

    if(cl > Tau)
    {
        Real lortho = (eigenValues[1] + eigenValues[0]) / 2.0;
        Real lbar   = tensor.GetTrace() / 3.0;
        Matrix S    = m_signal->signalAt(xk);

        unsigned int i=0;
        for(std::vector<Vector>::iterator g = m_directions->begin(); g != m_directions->end(); g++, i++)
        {
            Real mesuredSignal   = S(i,0);

            Vector v = uk.toVector();
            Real scalarProd = v.x()*g->x() + v.y()*g->y() + v.z()*g->z();
            Real estimatedSignal = m_signal->BaselineAt(xk) * std::exp( -(Real)m_signal->GetBValue() * ( lortho + 3.0*(scalarProd*scalarProd) * (lbar - lortho) ) );

            density += this->computeNormalDensity(m_sigmas->at(i), std::log(mesuredSignal) - std::log(estimatedSignal));
        } // for i
    }
    else // cl <= Tau
    {
        Real moySigma = 0;
        for(std::vector<Real>::iterator sigma = m_sigmas->begin(); sigma != m_sigmas->end(); sigma++)
            moySigma += (*sigma)/m_sigmas->size();

        Vector v = uk.toVector();
        Real scalarProd = v.x()*eigenVectors(0,0) + v.y()*eigenVectors(0,1) + v.z()*eigenVectors(0,2);
        Real coeff      = 1.0 / ( moySigma * std::sqrt(2.0*M_PI) );
        Real exponent   = ( std::acos(scalarProd) - M_PI/2.0 ) / moySigma;
        density = std::log( coeff * std::exp(-0.5 * (exponent*exponent)) * 1.0/(2.0*M_PI) );
    }

    return density;
}

Real DTFPLikelihoodDensity::computeNormalDensity(Real sigma, Real x)
{
    if(sigma != 0)
    {
        Real fraction = x / sigma;

        return ( -std::log(sigma) - m_logSqrt2PI - 0.5 * fraction * fraction );
    }
    else
        return 1;
}

} // namespace btk

