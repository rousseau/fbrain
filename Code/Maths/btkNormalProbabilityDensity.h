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

#ifndef BTK_NORMAL_PROBABILITY_DENSITY_H
#define BTK_NORMAL_PROBABILITY_DENSITY_H


// Local includes
#include "btkMacro.h"
#include "btkProbabilityDensity.h"


namespace btk
{

class NormalProbabilityDensity : public ProbabilityDensity< double >
{
    public:
        typedef NormalProbabilityDensity     Self;
        typedef ProbabilityDensity< double > Superclass;

        /**
         * @brief Constructor (standard normal distribution).
         */
        NormalProbabilityDensity();

        /**
         * @brief Constructor.
         * @param mu Mean of the normal distribution.
         * @param sigma Standard deviation of the normal distribution.
         */
        NormalProbabilityDensity(double mu, double sigma);

        btkGetMacro(Mu, double);

        btkGetMacro(Sigma, double);

        /**
         * @brief Evaluate the distribution in x.
         * @param x Point of evaluation.
         * @return The value of the distribution in x.
         */
        double Evaluate(double x);

        /**
         * @brief Evaluate the distribution in x.
         * @param sigma Standard deviation of the normal law.
         * @param x Point of evaluation.
         * @return The value of the distribution in x.
         */
        double Evaluate(double sigma, double x);

        /**
         * @brief Simulate the distribution.
         * @return Simulated point.
         */
        double Simulate();

    private:
        /**
         * @brief Initialize the density by precomputing some constants.
         */
        void Initialize();

    private:
        /**
         * @brief Mean of the normal distribution.
         */
        double m_Mu;

        /**
         * @brief Standard deviation of the normal distribution.
         */
        double m_Sigma;

        /**
         * @brief Precomputed constant.
         */
        static const double m_sqrt2pi = 2.506628274631000;

        /**
         * @brief Precomputed constant.
         */
        static const double m_2pi = 6.283185307179586;

        /**
         * @brief Precomputed normalization constant.
         */
        double m_normalizationConstant;

        /**
         * @brief Precomputed variance.
         */
        double m_Sigma2;
};

} // namespace btk

#endif // BTK_NORMAL_PROBABILITY_DENSITY_H
