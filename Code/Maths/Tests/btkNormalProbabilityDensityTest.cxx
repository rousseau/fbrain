/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 01/03/2013
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

#include "btkNormalProbabilityDensityTest.h"


// Local includes
#include "btkMacro.h"
#include "btkNormalProbabilityDensity.h"


namespace btk
{

void NormalProbabilityDensityTest::setUp()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void NormalProbabilityDensityTest::tearDown()
{
    // ----
}

//-----------------------------------------------------------------------------------------------------------

void NormalProbabilityDensityTest::testEvaluate()
{
    NormalProbabilityDensity density;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, density.GetMu(), DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, density.GetSigma(), DBL_EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.398942280401433, density.Evaluate(0.0), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.241970724519143, density.Evaluate(-1.0), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.241970724519143, density.Evaluate(1.0), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.388077183263685, density.Evaluate(0.235), 1e-15);


    density = NormalProbabilityDensity(-0.86, 0.6589);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.86, density.GetMu(), DBL_EPSILON);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6589, density.GetSigma(), DBL_EPSILON);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.258325226256689, density.Evaluate(0.0), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.591953072477219, density.Evaluate(-1.0), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.011264495476859, density.Evaluate(1.0), 1e-15);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.152186964129487, density.Evaluate(0.235), 1e-15);
}

//-----------------------------------------------------------------------------------------------------------

void NormalProbabilityDensityTest::testSimulate()
{
    srand(time(NULL));
    btkTicTocInit();


    NormalProbabilityDensity density;

    std::vector< double > simulations(1000000, 0.0);
    double mean = 0.0;

    btkTic();
    for(unsigned int n = 0; n < simulations.size(); n++)
    {
        simulations[n] = density.Simulate();
        mean          += simulations[n];
    }
    btkToc();

    mean /= simulations.size();

    double variance = 0.0;

    for(unsigned int n = 0; n < simulations.size(); n++)
    {
        variance += (simulations[n]-mean) * (simulations[n]-mean);
    }

    variance /= simulations.size()-1;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu(), mean, 1e-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetSigma(), std::sqrt(variance), 1e-2);


    density = NormalProbabilityDensity(-3.24, 2.61);

    simulations = std::vector< double >(1000000, 0.0);
    mean = 0.0;

    btkTic();
    for(unsigned int n = 0; n < simulations.size(); n++)
    {
        simulations[n] = density.Simulate();
        mean          += simulations[n];
    }
    btkToc();

    mean /= simulations.size();

    variance = 0.0;

    for(unsigned int n = 0; n < simulations.size(); n++)
    {
        variance += (simulations[n]-mean) * (simulations[n]-mean);
    }

    variance /= simulations.size()-1;

    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetMu(), mean, 1e-2);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(density.GetSigma(), std::sqrt(variance), 1e-2);
}

} // namespace btk
