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


#ifndef BTK_DTFPLIKELIHOODDENSITY_H
#define BTK_DTFPLIKELIHOODDENSITY_H

    // STL includes
    #include "vector"
    #include "string"

    // Local includes
    #include "btkTypes.h"
    #include "btkPoint.h"
    #include "btkDirection.h"
    #include "btkDTFPSignal.h"


    namespace btk
    {

    /**
     * @class LikelihoodDensity
     * @brief Likelihood density
     * @author Julien Pontabry
     * Likelihood density for particle filter (using MRI data and associated Q-ball HS model)
     */
    class DTFPLikelihoodDensity
    {
        public:
            /**
             * @brief Constructor
             * Build likelihood density with normal density, MRI data and associated Q-ball HS model
             * @param d Normal density
             * @param signal Signal data
             * @param model Q-ball HS model
             */
            DTFPLikelihoodDensity(DTFPSignal *signal);

            /**
             * @brief Compute density
             * Compute density with given point, current vector and mean direction
             * @param uk Direction we want to compute probability
             * @param xk Point in euclidian space
             * @param mean Mean direction
             * @return Likelihood probability
             */
            Real compute(Direction uk, Point xk, Direction mean);

        private:
            /**
             * @brief Compute a normal density
             * @param sigma Standard deviation of the normal density
             * @param x Variable of the normal density
             * @return Value of the normal density N(x | 0,sigma)
             */
			Real computeNormalDensity(Real sigma, Real x);


        private:
            DTFPSignal  *m_signal;     /**< Signal data */

            std::vector<Vector> *m_directions;    /**< Gradients directions */
            std::vector<Real>      *m_sigmas;        /**< noise's standard deviations */

			Real m_logSqrt2PI;  /**< Precomputed constant */
    };

    } // namespace btk

#endif /* BTK_DTFPLIKELIHOODDENSITY_H */

