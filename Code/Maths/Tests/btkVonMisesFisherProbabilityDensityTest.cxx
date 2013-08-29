/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 04/03/2013
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

#include "btkVonMisesFisherProbabilityDensityTest.h"


// Local includes
#include "btkMacro.h"
#include "btkVonMisesFisherProbabilityDensity.h"


namespace btk
{

void VonMisesFisherProbabilityDensityTest::setUp()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void VonMisesFisherProbabilityDensityTest::tearDown()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void VonMisesFisherProbabilityDensityTest::testEvaluate()
{
    VonMisesFisherProbabilityDensity density;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, density.GetMu()[0], DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, density.GetMu()[1], DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.GetMu()[2], DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.GetKappa(), DBL_EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.184065499616596, density.Evaluate(GradientDirection(0.0,0.0,1.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.067713913137896, density.Evaluate(GradientDirection(0.707106781186547,0.707106781186547,0.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.067713913137896, density.Evaluate(GradientDirection(1.0,0.0,0.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.067713913137896, density.Evaluate(GradientDirection(0.0,1.0,0.0)), 1e-15);

    density = VonMisesFisherProbabilityDensity(GradientDirection(0.0,0.0,1.0), 0.1);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, density.GetMu()[0], DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, density.GetMu()[1], DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.GetMu()[2], DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.1, density.GetKappa(), DBL_EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.087800300268453, density.Evaluate(GradientDirection(0.0,0.0,1.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.079444996997689, density.Evaluate(GradientDirection(0.707106781186547,0.707106781186547,0.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.079444996997689, density.Evaluate(GradientDirection(1.0,0.0,0.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.079444996997689, density.Evaluate(GradientDirection(0.0,1.0,0.0)), 1e-15);


    density = VonMisesFisherProbabilityDensity(GradientDirection(0.707106781186547, 0.707106781186547, 0.0), 0.6589);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.707106781186547, density.GetMu()[0], 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.707106781186547, density.GetMu()[1], 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, density.GetMu()[2], 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6589, density.GetKappa(), DBL_EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.074098257445687, density.Evaluate(GradientDirection(0.0,0.0,1.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.143207125987788, density.Evaluate(GradientDirection(0.707106781186547,0.707106781186547,0.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.118073188520780, density.Evaluate(GradientDirection(1.0,0.0,0.0)), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.118073188520780, density.Evaluate(GradientDirection(0.0,1.0,0.0)), 1e-15);
}

//-----------------------------------------------------------------------------------------------------------

void VonMisesFisherProbabilityDensityTest::testSimulate()
{
    srand(time(NULL));
    btkTicTocInit();


    VonMisesFisherProbabilityDensity density;

    std::vector< GradientDirection > simulations(1000000);
    GradientDirection mean(0.0, 0.0, 0.0);

    btkTic();
    for(unsigned int n = 0; n < simulations.size(); n++)
    {
        simulations[n] = density.Simulate();
        mean[0]       += simulations[n][0];
        mean[1]       += simulations[n][1];
        mean[2]       += simulations[n][2];
    }
    btkToc();

    double Rbar = std::sqrt(mean[0]*mean[0] + mean[1]*mean[1] + mean[0]*mean[0]) / simulations.size();
    double kappa = ( Rbar * (3.0 - Rbar*Rbar) ) / ( 1 - Rbar*Rbar );

    mean.Normalize();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[0], mean[0], 1e-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[1], mean[1], 1e-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[2], mean[2], 1e-2);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetKappa(), kappa, 1e-2);


    density = VonMisesFisherProbabilityDensity(GradientDirection(0.0,0.0,1.0), 0.1);

    std::vector< GradientDirection > simulations2(1000000);
    mean = GradientDirection(0.0,0.0,0.0);

    btkTic();
    for(unsigned int n = 0; n < simulations2.size(); n++)
    {
        simulations2[n] = density.Simulate();
        mean[0]        += simulations2[n][0];
        mean[1]        += simulations2[n][1];
        mean[2]        += simulations2[n][2];
    }
    btkToc();

    Rbar = std::sqrt(mean[0]*mean[0] + mean[1]*mean[1] + mean[0]*mean[0]) / simulations2.size();
    kappa = ( Rbar * (3.0 - Rbar*Rbar) ) / ( 1 - Rbar*Rbar );

    mean.Normalize();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[0], mean[0], 1e-1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[1], mean[1], 1e-1);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[2], mean[2], 1e-2);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetKappa(), kappa, 1e-2);


    density = VonMisesFisherProbabilityDensity(GradientDirection(0.707106781186547, 0.707106781186547, 0.0), 12);

    std::vector< GradientDirection > simulations3(1000000);
    mean = GradientDirection(0.0,0.0,0.0);

    btkTic();
    for(unsigned int n = 0; n < simulations3.size(); n++)
    {
        simulations3[n] = density.Simulate();
        mean[0]        += simulations3[n][0];
        mean[1]        += simulations3[n][1];
        mean[2]        += simulations3[n][2];
    }
    btkToc();

    Rbar = std::sqrt(mean[0]*mean[0] + mean[1]*mean[1] + mean[0]*mean[0]) / simulations3.size();
    kappa = ( Rbar * (3.0 - Rbar*Rbar) ) / ( 1 - Rbar*Rbar );

    mean.Normalize();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[0], mean[0], 1e-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[1], mean[1], 1e-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu()[2], mean[2], 1e-2);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetKappa(), kappa, 1e-2);
}

} // namespace btk
