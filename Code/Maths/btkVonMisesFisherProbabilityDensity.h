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

#ifndef BTK_VON_MISES_FISHER_PROBABILITY_DENSITY_H
#define BTK_VON_MISES_FISHER_PROBABILITY_DENSITY_H

// Local includes
#include "btkMacro.h"
#include "btkProbabilityDensity.h"
#include "btkGradientDirection.h"


namespace btk
{

class VonMisesFisherProbabilityDensity : public ProbabilityDensity< GradientDirection >
{
    public:
        typedef VonMisesFisherProbabilityDensity        Self;
        typedef ProbabilityDensity< GradientDirection > Superclass;

        /**
         * @brief Constructor (mean to (0,0) and concentration to 1).
         */
        VonMisesFisherProbabilityDensity();

        /**
         * @brief Constructor.
         * @param mu Mean direction.
         * @param kappa Concentration.
         */
        VonMisesFisherProbabilityDensity(GradientDirection mu, double kappa);

        btkGetMacro(Mu, GradientDirection);

        btkGetMacro(Kappa, double);

        /**
         * @brief Evaluate the distribution in x.
         * @param x Point of evaluation.
         * @return The value of the distribution in x.
         */
        double Evaluate(GradientDirection x);

        /**
         * @brief Simulate the distribution.
         * @return Simulated point.
         */
        GradientDirection Simulate();

    private:
        /**
         * @brief Initialize the density.
         */
        void Initialize();

    private:
        /**
         * @brief Mean direction of the distribution.
         */
        GradientDirection m_Mu;

        /**
         * @brief Concentration of the distribution.
         */
        double m_Kappa;

        /**
         * @brief Precomputed constant.
         */
        static const double m_2pi = 6.283185307179586;

        /**
         * @brief Precomputed normalization constant.
         */
        double m_NormalizationCoefficient;

        /**
         * @brief Precomputed constant.
         */
        double m_SimulationTerm;

        /**
         * @brief Precomputed constant.
         */
        double m_SimulationCoefficient;

        /**
         * @brief Precomputed constant.
         */
        double m_InverseKappa;
};

} // namespace btk

#endif // BTK_VON_MISES_FISHER_PROBABILITY_DENSITY_H
