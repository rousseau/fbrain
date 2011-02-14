/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

15 june 2010
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


#ifndef BTK_APRIORIDENSITY_H
#define BTK_APRIORIDENSITY_H

    // Local includes
    #include "btkTypes.h"
    #include "btkDirection.h"


    namespace btk
    {

        /**
         * @class APrioriDensity
         * @brief A priori density
         * @author Julien Pontabry
         */
        class APrioriDensity
        {
            public:
                /**
                 * @brief Constructor
                 * Build an a priori density for particle filter using von Mises-Fisher density
                 * @param K Concentration of the vMF density
                 */
                APrioriDensity(Real K);

                /**
                 * @brief Compute density
                 * Compute density with given last vector and current vector
                 * @param uk Current vector (wich had been simulated)
                 * @param ukm1 Last vector
                 */
                Real compute(Direction uk, Direction ukm1);

                /**
                 * @brief Simulate the density
                 * @param mean Mean direction of the vMF density
                 * @return A simulated direction
                 */
                Direction simulate(Direction mean);

            private:
                Real m_K;   /**< Concentration */
                Real m_C;   /**< Coefficient (for vMF computing) */

                Real m_1OnK;    /**< Precomputed constant */
                Real m_c;       /**< Precomputed constant */
                Real m_Kc;      /**< Precomputed constant */
                Real m_emK;     /**< Precomputed constant */
                Real m_2PI;     /**< Precomputed constant */
        };

    } // namespace btk

#endif /* BTK_APRIORIDENSITY_H */

