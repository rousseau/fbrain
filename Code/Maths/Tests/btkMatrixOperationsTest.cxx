/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/07/2012
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

#include "btkMatrixOperationsTest.h"


#define EPSILON     1e-5
#define EPSILON_BIS 1e-1


namespace btk
{

void MatrixOperationsTest::setUp()
{
    // Null matrix
    O = MatrixOperations::Matrix(4,4);
    O.Fill(0);
    O(3,3) = 1;

    // Idendity matrix
    I = MatrixOperations::Matrix(4,4);
    I.SetIdentity();

    // General matrix
    M1 = MatrixOperations::Matrix(4,4);
    M1(0,0) = 2; M1(0,1) = 5; M1(0,2) = 1; M1(0,3) = 0;
    M1(1,0) = 0; M1(1,1) = 3; M1(1,2) = 9; M1(1,3) = 0;
    M1(2,0) = 4; M1(2,1) = 2; M1(2,2) = 6; M1(2,3) = 0;
    M1(3,0) = 0; M1(3,1) = 0; M1(3,2) = 0; M1(3,3) = 1;

    M2 = MatrixOperations::Matrix(4,4);
    M2(0,0) = 0.89; M2(0,1) = 0.3; M2(0,2) = 0;    M2(0,3) = 0;
    M2(1,0) = 0.4;  M2(1,1) = 1.2; M2(1,2) = 0.25; M2(1,3) = 0;
    M2(2,0) = 0.23; M2(2,1) = 0.1; M2(2,2) = 1;    M2(2,3) = 0;
    M2(3,0) = 0;    M2(3,1) = 0;   M2(3,2) = 0;    M2(3,3) = 1;
}

//-----------------------------------------------------------------------------------------------------------

void MatrixOperationsTest::tearDown()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void MatrixOperationsTest::testNorm()
{
    CPPUNIT_ASSERT_DOUBLES_EQUAL(MatrixOperations::Norm(O), 1.0, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(MatrixOperations::Norm(I), 1.0, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(MatrixOperations::Norm(M1), 12.0493728236864, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(MatrixOperations::Norm(M2), 1.51768453824806, EPSILON);
}

//-----------------------------------------------------------------------------------------------------------

void MatrixOperationsTest::testSqrt()
{
    MatrixOperations::Matrix Osqrt_target(4,4); Osqrt_target.Fill(0); Osqrt_target(3,3) = 1;
    MatrixOperations::Matrix Osqrt = MatrixOperations::Sqrt(O, EPSILON);
    CPPUNIT_ASSERT_EQUAL(Osqrt, Osqrt_target);

    MatrixOperations::Matrix Isqrt_target(4,4); Isqrt_target.SetIdentity();
    MatrixOperations::Matrix Isqrt = MatrixOperations::Sqrt(I, EPSILON);
    CPPUNIT_ASSERT_EQUAL(Isqrt, Isqrt_target);

    MatrixOperations::Matrix M1sqrt = MatrixOperations::Sqrt(M1, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(0,0), 1.777995186215568, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(0,1), 1.394613259510546, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(0,2),-0.430201439142169, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(0,3), 0.0,               EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(1,0),-0.536312337764582, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(1,1), 1.854310954962167, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(1,2), 2.027622763130861, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(1,3), 0.0,               EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(2,0), 0.960755932254226, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(2,1), 0.152631537493201, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(2,2), 2.470594949587503, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(2,3), 0.0,               EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(3,0), 0.0,               EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(3,1), 0.0,               EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(3,2), 0.0,               EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1sqrt(3,3), 1.0,               EPSILON);

    MatrixOperations::Matrix M2sqrt = MatrixOperations::Sqrt(M2, EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(0,0), 0.928626288333147,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(0,1), 0.149540944053739,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(0,2),-0.009404288153312,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(0,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(1,0), 0.192177971154113,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(1,1), 1.080017167804240,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(1,2), 0.121169214388569,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(1,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(2,0), 0.115400400160132,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(2,1), 0.039815740654381,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(2,2), 0.998128657338684,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(2,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(3,0), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(3,1), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(3,2), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2sqrt(3,3), 1.0,              EPSILON);
}

//-----------------------------------------------------------------------------------------------------------

void MatrixOperationsTest::testExponential()
{
    MatrixOperations::Matrix Oexp = MatrixOperations::Exponential(O);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(0,0), 1.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(0,1), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(0,2), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(0,3), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(1,0), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(1,1), 1.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(1,2), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(1,3), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(2,0), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(2,1), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(2,2), 1.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(2,3), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(3,0), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(3,1), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(3,2), 0.0,             EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Oexp(3,3), 2.71828182845904,EPSILON);

    MatrixOperations::Matrix Iexp = MatrixOperations::Exponential(I);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(0,0), 2.718281828459045,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(0,1), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(0,2), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(0,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(1,0), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(1,1), 2.718281828459045,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(1,2), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(1,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(2,0), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(2,1), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(2,2), 2.718281828459045,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(2,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(3,0), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(3,1), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(3,2), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(Iexp(3,3), 2.718281828459045,EPSILON);

    MatrixOperations::Matrix M1exp = MatrixOperations::Exponential(M1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(0,0), 1.10895400324639e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(0,1), 1.30732374478052e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(0,2), 2.52364437362355e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(0,3), 0.0,                 EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(1,0), 1.71089215310983e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(1,1), 2.01675578782176e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(1,2), 3.89298567840378e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(1,3), 0.0,                 EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(2,0), 1.54011672894503e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(2,1), 1.81560356915075e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(2,2), 3.50451680874633e+04,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(2,3), 0.0,                 EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(3,0), 0.0,                 EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(3,1), 0.0,                 EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(3,2), 0.0,                 EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1exp(3,3), 2.71828182845904,    EPSILON_BIS);

    MatrixOperations::Matrix M2exp = MatrixOperations::Exponential(M2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(0,0), 2.607497055239036,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(0,1), 0.877966296064564,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(0,2), 0.106578075003882,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(0,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(1,0), 1.252331585589062,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(1,1), 3.550254919507046,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(1,2), 0.770717207555226,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(1,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(2,0), 0.664581247649188,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(2,1), 0.406338712025662,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(2,2), 2.763156233456653,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(2,3), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(3,0), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(3,1), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(3,2), 0.0,              EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2exp(3,3), 2.71828182845904, EPSILON);
}

//-----------------------------------------------------------------------------------------------------------

void MatrixOperationsTest::testLogarithm()
{
    MatrixOperations::Matrix     Ilog = MatrixOperations::Logarithm(I,EPSILON);
    MatrixOperations::Matrix IlogTrue(4,4); IlogTrue.Fill(0);
    CPPUNIT_ASSERT_EQUAL(Ilog, IlogTrue);

    MatrixOperations::Matrix M1log = MatrixOperations::Logarithm(M1,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(0,0), 1.609371909146170,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(0,1), 1.347314321936097,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(0,2),-0.922198681580306,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(0,3), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(1,0),-0.901256631404012,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(1,1), 1.538360046114096,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(1,2), 1.614034811221365,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(1,3), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(2,0), 0.817488430698831,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(2,1),-0.142023726064148,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(2,2), 1.976232024142993,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(2,3), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(3,0), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(3,1), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(3,2), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M1log(3,3), 0.0,              EPSILON_BIS);

    MatrixOperations::Matrix M2log = MatrixOperations::Logarithm(M2,EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(0,0),-0.176000693769881,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(0,1), 0.301832757864699,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(0,2),-0.038109149959553,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(0,3), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(1,0), 0.373226662183941,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(1,1), 0.123190106037123,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(1,2), 0.237553943235413,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(1,3), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(2,0), 0.234453846359700,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(2,1), 0.059961159331376,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(2,2),-0.005878408615922,EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(2,3), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(3,0), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(3,1), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(3,2), 0.0,              EPSILON_BIS);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(M2log(3,3), 0.0,              EPSILON_BIS);
}

} // namespace btk
