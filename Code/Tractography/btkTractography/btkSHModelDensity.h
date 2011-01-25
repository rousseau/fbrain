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


#ifndef BTK_SH_MODEL_DENSITY_H
#define BTK_SH_MODEL_DENSITY_H

    // Local includes
    #include "btkSimDensity.h"
    #include "btkTypes.h"
    #include "btkPoint.h"
    #include "btkSHModel.h"


    namespace btk
    {

        /**
         * @class ODFDensity
         * @brief ODF Density
         * @author Julien Pontabry
         * ODF probability density using Q-ball HS model
         */
        class SHModelDensity : public SimDensity
        {
            public:
                /**
                 * @brief Constructor
                 * @param model Q-ball HS model of MRI data
                 */
                SHModelDensity(SHModel *model);

                /**
                 * @brief Compute density
                 * Compute density for a given direction
                 * @param u Given direction
                 * @return Value of density for the given direction
                 */
                Real compute(Direction u);

                /**
                 * @brief Compute density
                 * Compute density for a given direction
                 * @param p Location in MRI image
                 * @param u Given direction
                 * @return Value of density for the given direction
                 */
                Real compute(Point p, Direction u);

                /**
                 * @brief Simulate a direction
                 * Simulate a direction following ODF probability
                 * @return Simulated direction
                 */
                Direction simulate();

                /**
                 * @brief Simulate a direction
                 * Simulate a direction following ODF probability
                 * @param p Location in MRI image
                 * @return Simulated direction
                 */
                Direction simulate(Point p);

                /**
                 * @brief Set the point in model image
                 * @param p Point in euclidian space
                 */
                void setPoint(Point p);

            private:
                SHModel *m_model;  /**< Spherical harmonics model */

                Real m_pi;      /**< Pi number */
                Real m_c;       /**< C constant (for acceptation-rejet simulation) */
                Real m_norm;    /**< Normalisation constant */
                Real m_min;     /**< Min value (for min-max normalisation) */
                Real m_max;     /**< Max value (for min-max normalisation) */

                Point m_p;  /**< Point in model image */
        };

    } // namespace btk

#endif // BTK_SH_MODEL_DENSITY_H

