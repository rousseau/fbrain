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

#include "btkSphericalHarmonicsTest.h"


// STL includes
#include "cmath"

// Local includes
#include "btkSphericalHarmonics.h"


#define EPSILON 1e-7


namespace btk
{

void SphericalHarmonicsTest::setUp()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void SphericalHarmonicsTest::tearDown()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void SphericalHarmonicsTest::testComputeBasis()
{
    btk::SphericalDirection u(0,0);

    double cosTheta = std::cos(u[0]);
    double sinTheta = std::sin(u[0]);

    double cosTheta2 = cosTheta*cosTheta;
    double cosTheta3 = cosTheta2*cosTheta;
    double cosTheta4 = cosTheta3*cosTheta;
    double cosTheta5 = cosTheta4*cosTheta;
    double cosTheta6 = cosTheta5*cosTheta;
    double cosTheta7 = cosTheta6*cosTheta;
    double cosTheta8 = cosTheta7*cosTheta;

    double sinTheta2 = sinTheta*sinTheta;
    double sinTheta3 = sinTheta2*sinTheta;
    double sinTheta4 = sinTheta3*sinTheta;
    double sinTheta5 = sinTheta4*sinTheta;
    double sinTheta6 = sinTheta5*sinTheta;
    double sinTheta7 = sinTheta6*sinTheta;
    double sinTheta8 = sinTheta7*sinTheta;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/2.0 * std::sqrt(1.0/M_PI), btk::SphericalHarmonics::ComputeBasis(u, 0, 0), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/4.0 * std::sqrt(15.0/(2.0*M_PI)) * std::sin(2.0*u[1]) * sinTheta2, btk::SphericalHarmonics::ComputeBasis(u, 2, -2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/2.0 * std::sqrt(15.0/(2.0*M_PI)) * std::sin(u[1]) * sinTheta * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 2, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/4.0 * std::sqrt(5.0/M_PI) * (3.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 2, 0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/2.0 * std::sqrt(15.0/(2.0*M_PI)) * std::cos(u[1]) * sinTheta * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 2, 1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/4.0 * std::sqrt(15.0/(2.0*M_PI)) * std::cos(2.0*u[1]) * sinTheta2, btk::SphericalHarmonics::ComputeBasis(u, 2, 2), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/16.0 * std::sqrt(35.0/(2.0*M_PI)) * std::sin(4.0*u[1]) * sinTheta4, btk::SphericalHarmonics::ComputeBasis(u, 4, 4), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(35.0/M_PI) * std::sin(3.0*u[1]) * sinTheta3 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 4, 3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/8.0 * std::sqrt(5.0/(2.0*M_PI)) * std::sin(2.0*u[1]) * sinTheta2 * (7.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 4, 2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(1.0/M_PI) * std::sin(u[1]) * sinTheta * (7.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 4, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/16.0 * std::sqrt(1.0/M_PI) * (35.0*cosTheta4 - 30.0*cosTheta2 + 3.0), btk::SphericalHarmonics::ComputeBasis(u, 4, 0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(1.0/M_PI) * std::cos(u[1]) * sinTheta * (7.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 4, 1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/8.0 * std::sqrt(5.0/(2.0*M_PI)) * std::cos(2.0*u[1]) * sinTheta2 * (7.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 4, 2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(35.0/M_PI) * std::cos(3.0*u[1]) * sinTheta3 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 4, 3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/16.0 * std::sqrt(35.0/(2.0*M_PI)) * std::cos(4.0*u[1]) * sinTheta4, btk::SphericalHarmonics::ComputeBasis(u, 4, 4), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(3003.0/M_PI) * std::sin(6.0*u[1]) *sinTheta6, btk::SphericalHarmonics::ComputeBasis(u, 6, -6), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/32.0 * std::sqrt(1001.0/M_PI) * std::sin(5.0*u[1]) * sinTheta5 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 6, -5), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/32.0 * std::sqrt(91.0/(2.0*M_PI)) * std::sin(4.0*u[1]) * sinTheta4 * (11.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, -4), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/32.0 * std::sqrt(1365.0/M_PI) * std::sin(3.0*u[1]) * sinTheta3 * (11.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, -3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(1365.0/M_PI) * std::sin(2.0*u[1]) * sinTheta2 * (33.0*cosTheta4 - 18.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, -2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/16.0 * std::sqrt(273.0/(2.0*M_PI)) * std::sin(u[1]) * sinTheta * (33.0*cosTheta5 - 30.0*cosTheta3 + 5.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/32.0 * std::sqrt(13.0/M_PI) * (231.0*cosTheta6 - 315.0*cosTheta4 + 105.0*cosTheta2 - 5.0), btk::SphericalHarmonics::ComputeBasis(u, 6, 0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/16.0 * std::sqrt(273.0/(2.0*M_PI)) * std::cos(u[1]) * sinTheta * (33.0*cosTheta5 - 30.0*cosTheta3 + 5.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, 1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(1365.0/M_PI) * std::cos(2.0*u[1]) * sinTheta2 * (33.0*cosTheta4 - 18.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, 2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/32.0 * std::sqrt(1365.0/M_PI) * std::cos(3.0*u[1]) * sinTheta3 * (11.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, 3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/32.0 * std::sqrt(91.0/(2.0*M_PI)) * std::cos(4.0*u[1]) * sinTheta4 * (11.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, 4), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/32.0 * std::sqrt(1001.0/M_PI) * std::cos(5.0*u[1]) * sinTheta5 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 6, 5), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(3003.0/M_PI) * std::cos(6.0*u[1]) *sinTheta6, btk::SphericalHarmonics::ComputeBasis(u, 6, 6), EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/256.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::sin(8.0*u[1]) * sinTheta8, btk::SphericalHarmonics::ComputeBasis(u, 8, -8), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::sin(7.0*u[1]) * sinTheta7 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 8, -7), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * std::sqrt(7293.0/M_PI) * std::sin(6.0*u[1]) * sinTheta6 * (15.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, -6), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17017.0/(2.0*M_PI)) * std::sin(5.0*u[1]) * sinTheta5 * (5.0*cosTheta3 - 1.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, -5), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(1309.0/(2.0*M_PI)) * std::sin(4.0*u[1]) * sinTheta4 * (65.0*cosTheta4 - 26.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, -4), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/64.0 * std::sqrt(19635.0/(2.0*M_PI)) * std::sin(3.0*u[1]) * sinTheta3 * (39.0*cosTheta5 - 26.0*cosTheta3 + 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, -3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(595.0/M_PI) * std::sin(2.0*u[1]) * sinTheta2 * (143.0*cosTheta6 - 143.0*cosTheta4 +33.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, -2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17.0/(2.0*M_PI)) * std::sin(u[1]) * sinTheta * (715.0*cosTheta7 - 1001.0*cosTheta5 + 385.0*cosTheta3 -35.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/256.0 * std::sqrt(17.0/M_PI) * (6435.0*cosTheta8 - 12012.0*cosTheta6 + 6930.0*cosTheta4 - 1260.0*cosTheta2 + 35.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 0), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17.0/(2.0*M_PI)) * std::cos(u[1]) * sinTheta * (715.0*cosTheta7 - 1001.0*cosTheta5 + 385.0*cosTheta3 -35.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, 1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(595.0/M_PI) * std::cos(2.0*u[1]) * sinTheta2 * (143.0*cosTheta6 - 143.0*cosTheta4 +33.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 2), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/64.0 * std::sqrt(19635.0/(2.0*M_PI)) * std::cos(3.0*u[1]) * sinTheta3 * (39.0*cosTheta5 - 26.0*cosTheta3 + 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, 3), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(1309.0/(2.0*M_PI)) * std::cos(4.0*u[1]) * sinTheta4 * (65.0*cosTheta4 - 26.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 4), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17017.0/(2.0*M_PI)) * std::cos(5.0*u[1]) * sinTheta5 * (5.0*cosTheta3 - 1.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, 5), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * std::sqrt(7293.0/M_PI) * std::cos(6.0*u[1]) * sinTheta6 * (15.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 6), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::cos(7.0*u[1]) * sinTheta7 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 8, 7), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/256.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::cos(8.0*u[1]) * sinTheta8, btk::SphericalHarmonics::ComputeBasis(u, 8, 8), EPSILON);


    u = btk::SphericalDirection(1.02554,2.5695);

    cosTheta = std::cos(u[0]);
    sinTheta = std::sin(u[0]);

    cosTheta2 = cosTheta*cosTheta;
    cosTheta3 = cosTheta2*cosTheta;
    cosTheta4 = cosTheta3*cosTheta;
    cosTheta5 = cosTheta4*cosTheta;
    cosTheta6 = cosTheta5*cosTheta;
    cosTheta7 = cosTheta6*cosTheta;
    cosTheta8 = cosTheta7*cosTheta;

    sinTheta2 = sinTheta*sinTheta;
    sinTheta3 = sinTheta2*sinTheta;
    sinTheta4 = sinTheta3*sinTheta;
    sinTheta5 = sinTheta4*sinTheta;
    sinTheta6 = sinTheta5*sinTheta;
    sinTheta7 = sinTheta6*sinTheta;
    sinTheta8 = sinTheta7*sinTheta;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/2.0 * std::sqrt(1.0/M_PI), btk::SphericalHarmonics::ComputeBasis(u, 0, 0), EPSILON);

//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/4.0 * std::sqrt(15.0/(2.0*M_PI)) * std::sin(2.0*u[1]) * sinTheta2, btk::SphericalHarmonics::ComputeBasis(u, 2, -2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/2.0 * std::sqrt(15.0/(2.0*M_PI)) * std::sin(u[1]) * sinTheta * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 2, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/4.0 * std::sqrt(5.0/M_PI) * (3.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 2, 0), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/2.0 * std::sqrt(15.0/(2.0*M_PI)) * std::cos(u[1]) * sinTheta * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 2, 1), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/4.0 * std::sqrt(15.0/(2.0*M_PI)) * std::cos(2.0*u[1]) * sinTheta2, btk::SphericalHarmonics::ComputeBasis(u, 2, 2), EPSILON);

//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/16.0 * std::sqrt(35.0/(2.0*M_PI)) * std::sin(4.0*u[1]) * sinTheta4, btk::SphericalHarmonics::ComputeBasis(u, 4, 4), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(35.0/M_PI) * std::sin(3.0*u[1]) * sinTheta3 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 4, 3), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/8.0 * std::sqrt(5.0/(2.0*M_PI)) * std::sin(2.0*u[1]) * sinTheta2 * (7.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 4, 2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(1.0/M_PI) * std::sin(u[1]) * sinTheta * (7.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 4, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/16.0 * std::sqrt(1.0/M_PI) * (35.0*cosTheta4 - 30.0*cosTheta2 + 3.0), btk::SphericalHarmonics::ComputeBasis(u, 4, 0), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(1.0/M_PI) * std::cos(u[1]) * sinTheta * (7.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 4, 1), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/8.0 * std::sqrt(5.0/(2.0*M_PI)) * std::cos(2.0*u[1]) * sinTheta2 * (7.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 4, 2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/8.0 * std::sqrt(35.0/M_PI) * std::cos(3.0*u[1]) * sinTheta3 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 4, 3), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/16.0 * std::sqrt(35.0/(2.0*M_PI)) * std::cos(4.0*u[1]) * sinTheta4, btk::SphericalHarmonics::ComputeBasis(u, 4, 4), EPSILON);

//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(3003.0/M_PI) * std::sin(6.0*u[1]) *sinTheta6, btk::SphericalHarmonics::ComputeBasis(u, 6, -6), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/32.0 * std::sqrt(1001.0/M_PI) * std::sin(5.0*u[1]) * sinTheta5 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 6, -5), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/32.0 * std::sqrt(91.0/(2.0*M_PI)) * std::sin(4.0*u[1]) * sinTheta4 * (11.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, -4), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/32.0 * std::sqrt(1365.0/M_PI) * std::sin(3.0*u[1]) * sinTheta3 * (11.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, -3), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(1365.0/M_PI) * std::sin(2.0*u[1]) * sinTheta2 * (33.0*cosTheta4 - 18.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, -2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/16.0 * std::sqrt(273.0/(2.0*M_PI)) * std::sin(u[1]) * sinTheta * (33.0*cosTheta5 - 30.0*cosTheta3 + 5.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/32.0 * std::sqrt(13.0/M_PI) * (231.0*cosTheta6 - 315.0*cosTheta4 + 105.0*cosTheta2 - 5.0), btk::SphericalHarmonics::ComputeBasis(u, 6, 0), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/16.0 * std::sqrt(273.0/(2.0*M_PI)) * std::cos(u[1]) * sinTheta * (33.0*cosTheta5 - 30.0*cosTheta3 + 5.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, 1), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(1365.0/M_PI) * std::cos(2.0*u[1]) * sinTheta2 * (33.0*cosTheta4 - 18.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, 2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/32.0 * std::sqrt(1365.0/M_PI) * std::cos(3.0*u[1]) * sinTheta3 * (11.0*cosTheta3 - 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 6, 3), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/32.0 * std::sqrt(91.0/(2.0*M_PI)) * std::cos(4.0*u[1]) * sinTheta4 * (11.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 6, 4), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/32.0 * std::sqrt(1001.0/M_PI) * std::cos(5.0*u[1]) * sinTheta5 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 6, 5), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/64.0 * std::sqrt(3003.0/M_PI) * std::cos(6.0*u[1]) *sinTheta6, btk::SphericalHarmonics::ComputeBasis(u, 6, 6), EPSILON);

//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/256.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::sin(8.0*u[1]) * sinTheta8, btk::SphericalHarmonics::ComputeBasis(u, 8, -8), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::sin(7.0*u[1]) * sinTheta7 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 8, -7), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * std::sqrt(7293.0/M_PI) * std::sin(6.0*u[1]) * sinTheta6 * (15.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, -6), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17017.0/(2.0*M_PI)) * std::sin(5.0*u[1]) * sinTheta5 * (5.0*cosTheta3 - 1.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, -5), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(1309.0/(2.0*M_PI)) * std::sin(4.0*u[1]) * sinTheta4 * (65.0*cosTheta4 - 26.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, -4), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/64.0 * std::sqrt(19635.0/(2.0*M_PI)) * std::sin(3.0*u[1]) * sinTheta3 * (39.0*cosTheta5 - 26.0*cosTheta3 + 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, -3), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(595.0/M_PI) * std::sin(2.0*u[1]) * sinTheta2 * (143.0*cosTheta6 - 143.0*cosTheta4 +33.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, -2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17.0/(2.0*M_PI)) * std::sin(u[1]) * sinTheta * (715.0*cosTheta7 - 1001.0*cosTheta5 + 385.0*cosTheta3 -35.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, -1), EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/256.0 * std::sqrt(17.0/M_PI) * (6435.0*cosTheta8 - 12012.0*cosTheta6 + 6930.0*cosTheta4 - 1260.0*cosTheta2 + 35.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 0), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17.0/(2.0*M_PI)) * std::cos(u[1]) * sinTheta * (715.0*cosTheta7 - 1001.0*cosTheta5 + 385.0*cosTheta3 -35.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, 1), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(595.0/M_PI) * std::cos(2.0*u[1]) * sinTheta2 * (143.0*cosTheta6 - 143.0*cosTheta4 +33.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 2), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.0/64.0 * std::sqrt(19635.0/(2.0*M_PI)) * std::cos(3.0*u[1]) * sinTheta3 * (39.0*cosTheta5 - 26.0*cosTheta3 + 3.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, 3), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/128.0 * std::sqrt(1309.0/(2.0*M_PI)) * std::cos(4.0*u[1]) * sinTheta4 * (65.0*cosTheta4 - 26.0*cosTheta2 + 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 4), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(17017.0/(2.0*M_PI)) * std::cos(5.0*u[1]) * sinTheta5 * (5.0*cosTheta3 - 1.0*cosTheta), btk::SphericalHarmonics::ComputeBasis(u, 8, 5), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0/128.0 * std::sqrt(7293.0/M_PI) * std::cos(6.0*u[1]) * sinTheta6 * (15.0*cosTheta2 - 1.0), btk::SphericalHarmonics::ComputeBasis(u, 8, 6), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0/64.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::cos(7.0*u[1]) * sinTheta7 * cosTheta, btk::SphericalHarmonics::ComputeBasis(u, 8, 7), EPSILON);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0/256.0 * std::sqrt(12155.0/(2.0*M_PI)) * std::cos(8.0*u[1]) * sinTheta8, btk::SphericalHarmonics::ComputeBasis(u, 8, 8), EPSILON);
}

} // namespace btk
