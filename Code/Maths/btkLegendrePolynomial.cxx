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

#include "btkLegendrePolynomial.h"


// STL includes
#include "cmath"


namespace btk
{

float LegendrePolynomial::ComputeInZero(unsigned int l)
{
    float odd = 1;

    for(unsigned int k=1; k<l; k+=2)
        odd *= k;


    float even = 2;

    for(unsigned int k=2; k<=l; k+=2)
        even *= k;


    float frac = odd / even;

    if((l/2)%2 != 0) // l/2 is odd
        return -frac;
    else // l/2 is even
        return frac;
}

//----------------------------------------------------------------------------------------

float LegendrePolynomial::Compute(unsigned int l, unsigned int m, float theta)
{
    float pmm   = 1;
    float x     = std::cos(theta);
    float somx2 = std::sin(theta);
    float fact  = 1;

    // pmm(n) = pmm(n-1) * -fact * sin(theta)
    for(unsigned int i=1; i<=m; i++)
    {
        pmm   = pmm * (-fact * somx2);
        fact += 2;
    } // for i

    float Plm;

    if(l == m)
        Plm = pmm;
    else // l <> m
    {
        float pmmp1 = x * (2. * m + 1.) * pmm;

        if(l == m+1)
            Plm = pmmp1;

        else // l < m OR l > m+1
        {
            float pk = 0;

            for(unsigned int k=m+2; k<=l; k++)
            {
                pk   = (x * (2. * k - 1.) * pmmp1 - (k + m - 1.) * pmm) / (k - m);
                pmm   = pmmp1;
                pmmp1 = pk;
            } // for k

            Plm = pk;
        } // else l < m OR l > m+1
    } // else l <> m

    return Plm;
}

} // namespace btk
