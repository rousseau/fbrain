/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 24/07/2012
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

#include "btkLegendrePolynomialTest.h"


// STL includes
#include "cmath"

// Local includes
#include "btkLegendrePolynomial.h"


#define EPSILON 1e-9


namespace btk
{

void LegendrePolynomialTest::setUp()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void LegendrePolynomialTest::tearDown()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void LegendrePolynomialTest::testComputeInZero()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, btk::LegendrePolynomial::ComputeInZero(0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5, btk::LegendrePolynomial::ComputeInZero(2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/8.0, btk::LegendrePolynomial::ComputeInZero(4), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-5.0/16.0, btk::LegendrePolynomial::ComputeInZero(6), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(35.0/128.0, btk::LegendrePolynomial::ComputeInZero(8), EPSILON);
}

//-----------------------------------------------------------------------------------------------------------

void LegendrePolynomialTest::testCompute()
{
    double theta = 0.276959;
    double     x = std::cos(theta);
    double    x2 = x*x;
    double     s = std::sqrt(1.0-x2);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, btk::LegendrePolynomial::Compute(0, 0, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * (3.0*x2 - 1.0), btk::LegendrePolynomial::Compute(2, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0*x*s, btk::LegendrePolynomial::Compute(2, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0 * (1.0-x2), btk::LegendrePolynomial::Compute(2, 2, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/8.0 * (35.0*x2*x2 - 30.0*x2 + 3.0), btk::LegendrePolynomial::Compute(4, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0/2.0 * (7.0*x2*x - 3.0*x)*s, btk::LegendrePolynomial::Compute(4, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(15.0/2.0 * (7.0*x2 - 1.0) * (1.0 - x2), btk::LegendrePolynomial::Compute(4, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0*x*s*s*s, btk::LegendrePolynomial::Compute(4, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0 * s*s*s*s, btk::LegendrePolynomial::Compute(4, 4, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/16.0 * (231.0*x2*x2*x2 - 315.0*x2*x2 + 105.0*x2 - 5), btk::LegendrePolynomial::Compute(6, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.0/8.0 * x * (33.0*x2*x2 - 30.0*x2 + 5.0)*s, btk::LegendrePolynomial::Compute(6, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0/8.0 * s*s * (33.0*x2*x2 - 18.0*x2 + 1.0), btk::LegendrePolynomial::Compute(6, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/2.0 * (11.0*x2 - 3.0)*x*s*s*s, btk::LegendrePolynomial::Compute(6, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(945.0/2.0 * s*s*s*s * (11.0 * x2 - 1.0), btk::LegendrePolynomial::Compute(6, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0*x*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0 * s*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 6, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * (6435.0*x2*x2*x2*x2 - 12012*x2*x2*x2 + 6930*x2*x2 -1260*x2 + 35), btk::LegendrePolynomial::Compute(8, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0/16.0 *x*s * (715*x2*x2*x2 - 1001*x2*x2 + 385*x2 - 35), btk::LegendrePolynomial::Compute(8, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/16.0 * s*s * (143.0*x2*x2*x2 - 143.0*x2*x2 + 33.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3465.0/8.0 * x*s*s*s * (39*x2*x2 - 26.0*x2 + 3.0), btk::LegendrePolynomial::Compute(8, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0/8.0 * s*s*s*s * (65.0*x2*x2 - 26.0*x2 + 1.0), btk::LegendrePolynomial::Compute(8, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * x*s*s*s*s*s * (5.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * s*s*s*s*s*s * (15.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 6, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *x*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 7, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *s*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 8, theta), EPSILON);


    theta = 1.0;
        x = std::cos(theta);
       x2 = x*x;
        s = std::sqrt(1.0-x2);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, btk::LegendrePolynomial::Compute(0, 0, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * (3.0*x2 - 1.0), btk::LegendrePolynomial::Compute(2, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0*x*s, btk::LegendrePolynomial::Compute(2, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0 * (1.0-x2), btk::LegendrePolynomial::Compute(2, 2, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/8.0 * (35.0*x2*x2 - 30.0*x2 + 3.0), btk::LegendrePolynomial::Compute(4, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0/2.0 * (7.0*x2*x - 3.0*x)*s, btk::LegendrePolynomial::Compute(4, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(15.0/2.0 * (7.0*x2 - 1.0) * (1.0 - x2), btk::LegendrePolynomial::Compute(4, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0*x*s*s*s, btk::LegendrePolynomial::Compute(4, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0 * s*s*s*s, btk::LegendrePolynomial::Compute(4, 4, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/16.0 * (231.0*x2*x2*x2 - 315.0*x2*x2 + 105.0*x2 - 5), btk::LegendrePolynomial::Compute(6, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.0/8.0 * x * (33.0*x2*x2 - 30.0*x2 + 5.0)*s, btk::LegendrePolynomial::Compute(6, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0/8.0 * s*s * (33.0*x2*x2 - 18.0*x2 + 1.0), btk::LegendrePolynomial::Compute(6, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/2.0 * (11.0*x2 - 3.0)*x*s*s*s, btk::LegendrePolynomial::Compute(6, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(945.0/2.0 * s*s*s*s * (11.0 * x2 - 1.0), btk::LegendrePolynomial::Compute(6, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0*x*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0 * s*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 6, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * (6435.0*x2*x2*x2*x2 - 12012*x2*x2*x2 + 6930*x2*x2 -1260*x2 + 35), btk::LegendrePolynomial::Compute(8, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0/16.0 *x*s * (715*x2*x2*x2 - 1001*x2*x2 + 385*x2 - 35), btk::LegendrePolynomial::Compute(8, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/16.0 * s*s * (143.0*x2*x2*x2 - 143.0*x2*x2 + 33.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3465.0/8.0 * x*s*s*s * (39*x2*x2 - 26.0*x2 + 3.0), btk::LegendrePolynomial::Compute(8, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0/8.0 * s*s*s*s * (65.0*x2*x2 - 26.0*x2 + 1.0), btk::LegendrePolynomial::Compute(8, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * x*s*s*s*s*s * (5.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * s*s*s*s*s*s * (15.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 6, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *x*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 7, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *s*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 8, theta), EPSILON);


    theta = 2.326514;
        x = std::cos(theta);
       x2 = x*x;
        s = std::sqrt(1.0-x2);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, btk::LegendrePolynomial::Compute(0, 0, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * (3.0*x2 - 1.0), btk::LegendrePolynomial::Compute(2, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0*x*s, btk::LegendrePolynomial::Compute(2, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0 * (1.0-x2), btk::LegendrePolynomial::Compute(2, 2, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/8.0 * (35.0*x2*x2 - 30.0*x2 + 3.0), btk::LegendrePolynomial::Compute(4, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0/2.0 * (7.0*x2*x - 3.0*x)*s, btk::LegendrePolynomial::Compute(4, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(15.0/2.0 * (7.0*x2 - 1.0) * (1.0 - x2), btk::LegendrePolynomial::Compute(4, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0*x*s*s*s, btk::LegendrePolynomial::Compute(4, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0 * s*s*s*s, btk::LegendrePolynomial::Compute(4, 4, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/16.0 * (231.0*x2*x2*x2 - 315.0*x2*x2 + 105.0*x2 - 5), btk::LegendrePolynomial::Compute(6, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.0/8.0 * x * (33.0*x2*x2 - 30.0*x2 + 5.0)*s, btk::LegendrePolynomial::Compute(6, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0/8.0 * s*s * (33.0*x2*x2 - 18.0*x2 + 1.0), btk::LegendrePolynomial::Compute(6, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/2.0 * (11.0*x2 - 3.0)*x*s*s*s, btk::LegendrePolynomial::Compute(6, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(945.0/2.0 * s*s*s*s * (11.0 * x2 - 1.0), btk::LegendrePolynomial::Compute(6, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0*x*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0 * s*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 6, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * (6435.0*x2*x2*x2*x2 - 12012*x2*x2*x2 + 6930*x2*x2 -1260*x2 + 35), btk::LegendrePolynomial::Compute(8, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0/16.0 *x*s * (715*x2*x2*x2 - 1001*x2*x2 + 385*x2 - 35), btk::LegendrePolynomial::Compute(8, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/16.0 * s*s * (143.0*x2*x2*x2 - 143.0*x2*x2 + 33.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3465.0/8.0 * x*s*s*s * (39*x2*x2 - 26.0*x2 + 3.0), btk::LegendrePolynomial::Compute(8, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0/8.0 * s*s*s*s * (65.0*x2*x2 - 26.0*x2 + 1.0), btk::LegendrePolynomial::Compute(8, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * x*s*s*s*s*s * (5.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * s*s*s*s*s*s * (15.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 6, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *x*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 7, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *s*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 8, theta), EPSILON);


    theta = 3.02155;
        x = std::cos(theta);
       x2 = x*x;
        s = std::sqrt(1.0-x2);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, btk::LegendrePolynomial::Compute(0, 0, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5 * (3.0*x2 - 1.0), btk::LegendrePolynomial::Compute(2, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0*x*s, btk::LegendrePolynomial::Compute(2, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0 * (1.0-x2), btk::LegendrePolynomial::Compute(2, 2, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/8.0 * (35.0*x2*x2 - 30.0*x2 + 3.0), btk::LegendrePolynomial::Compute(4, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0/2.0 * (7.0*x2*x - 3.0*x)*s, btk::LegendrePolynomial::Compute(4, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(15.0/2.0 * (7.0*x2 - 1.0) * (1.0 - x2), btk::LegendrePolynomial::Compute(4, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0*x*s*s*s, btk::LegendrePolynomial::Compute(4, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0 * s*s*s*s, btk::LegendrePolynomial::Compute(4, 4, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/16.0 * (231.0*x2*x2*x2 - 315.0*x2*x2 + 105.0*x2 - 5), btk::LegendrePolynomial::Compute(6, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(21.0/8.0 * x * (33.0*x2*x2 - 30.0*x2 + 5.0)*s, btk::LegendrePolynomial::Compute(6, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(105.0/8.0 * s*s * (33.0*x2*x2 - 18.0*x2 + 1.0), btk::LegendrePolynomial::Compute(6, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/2.0 * (11.0*x2 - 3.0)*x*s*s*s, btk::LegendrePolynomial::Compute(6, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(945.0/2.0 * s*s*s*s * (11.0 * x2 - 1.0), btk::LegendrePolynomial::Compute(6, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0*x*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0 * s*s*s*s*s*s, btk::LegendrePolynomial::Compute(6, 6, theta), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * (6435.0*x2*x2*x2*x2 - 12012*x2*x2*x2 + 6930*x2*x2 -1260*x2 + 35), btk::LegendrePolynomial::Compute(8, 0, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0/16.0 *x*s * (715*x2*x2*x2 - 1001*x2*x2 + 385*x2 - 35), btk::LegendrePolynomial::Compute(8, 1, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(315.0/16.0 * s*s * (143.0*x2*x2*x2 - 143.0*x2*x2 + 33.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 2, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3465.0/8.0 * x*s*s*s * (39*x2*x2 - 26.0*x2 + 3.0), btk::LegendrePolynomial::Compute(8, 3, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(10395.0/8.0 * s*s*s*s * (65.0*x2*x2 - 26.0*x2 + 1.0), btk::LegendrePolynomial::Compute(8, 4, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * x*s*s*s*s*s * (5.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 5, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(135135.0/2.0 * s*s*s*s*s*s * (15.0*x2 - 1.0), btk::LegendrePolynomial::Compute(8, 6, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *x*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 7, theta), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(2027025.0 *s*s*s*s*s*s*s*s, btk::LegendrePolynomial::Compute(8, 8, theta), EPSILON);
}

} // namespace btk
