/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 12/02/2010
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


#include "btkSphericalHarmonics.h"


// STL includes
#include "cmath"

// Local includes
#include "btkLegendrePolynomial.h"
#include "btkMathFunctions.h"


namespace btk
{

double SphericalHarmonics::ComputeBasis(btk::SphericalDirection u, unsigned int l, int m)
{
    // Order coefficient
    double orderCoefficient = 0.0;

    if(l == 0)
    {
        orderCoefficient = SphericalHarmonics::m_CoefficientOrder0;
    }
    else if(l == 2)
    {
        orderCoefficient = SphericalHarmonics::m_CoefficientOrder2;
    }
    else if(l == 4)
    {
        orderCoefficient = SphericalHarmonics::m_CoefficientOrder4;
    }
    else if(l == 6)
    {
        orderCoefficient = SphericalHarmonics::m_CoefficientOrder6;
    }
    else // if(l == 8)
    {
        orderCoefficient = SphericalHarmonics::m_CoefficientOrder8;
    }

    // Degree coefficient and spherical harmonics
    double degreeCoefficient = 0.0;
    double value = 0.0;

    if(m < 0)
    {
        degreeCoefficient = M_SQRT2 * std::sqrt( (double)btk::MathFunctions::factorial(l+m) / (double)btk::MathFunctions::factorial(l-m) );
                    value = (double)btk::LegendrePolynomial::Compute(l,(unsigned int)-m,u[0]) * (double)std::sin(-m*u[1]);
    }
    else if(m > 0)
    {
        degreeCoefficient = M_SQRT2 * std::sqrt( (double)btk::MathFunctions::factorial(l-m) / (double)btk::MathFunctions::factorial(l+m) );
                    value = (double)btk::LegendrePolynomial::Compute(l,(unsigned int)m,u[0]) * (double)std::cos(m*u[1]);
    }
    else // m = 0
    {
        degreeCoefficient = 1.0;
                    value = (double)btk::LegendrePolynomial::Compute(l,(unsigned int)0,u[0]);
    }

    return orderCoefficient * degreeCoefficient * value;
}



} // namespace btk

