/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

23 april 2010
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


#ifndef BTK_VONMISESFISHERDENSITY_H
#define BTK_VONMISESFISHERDENSITY_H

    // Local includes
    #include "btkSimDensity.h"


    namespace btk
    {

        /**
         * @class VonMisesFisherDensity
         * @brief Von Mises-Fisher density
         * @author Julien Pontabry
         */
        class VonMisesFisherDensity : public SimDensity
        {
            public:
                /**
                 * @brief Constructor
                 * Build von Mises-Fisher probability density with the given concentration
                 * @param K Concentration
                 */
                VonMisesFisherDensity(Real K);

                /**
                 * @brief Compute probability density
                 * Compute probability density for a given direction
                 * @param u Direction
                 * @return Value of the probability law for this direction
                 */
                Real compute(Direction u);

                /**
                 * @brief Compute probability density
                 * Compute probability density for a given direction
                 * @param u Direction
                 * @param kappa Concentration
                 * @return Value of the probability law for this direction
                 */
                Real compute(Direction u, Real kappa);

                /**
                 * @brief Compute probability density
                 * Compute probability density for a given direction
                 * @param u Direction
                 * @param mean Mean direction
                 * @param kappa Concentration
                 * @return Value of the probability law for this direction
                 */
                Real compute(Direction u, Direction mean, Real kappa);

                /**
                 * @brief Simulate a direction
                 * @return Simulated direction
                 */
                Direction simulate();

                /**
                 * @brief Simulate a direction
                 * @param mean Mean direction
                 * @param kappa Concentration
                 * @return Simulated direction
                 */
                Direction simulate(Direction mean, Real kappa);

                /**
                 * @brief Set the mean direction
                 * Set the mean direction for probability law
                 * @param u Mean direction
                 */
                void setMeanDirection(Direction u);

            private:
                Real m_K;   /**< Concentration */
                Real m_c;   /**< Coefficient (for rejection sampling) */
                Real m_C;   /**< Coefficient of vmf */

                Direction m_mean;   /**< Mean direction */
        };

    } // namespace btk

#endif // BTK_VONMISESFISHERDENSITY_H

