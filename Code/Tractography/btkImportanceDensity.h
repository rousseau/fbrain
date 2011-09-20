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


#ifndef BTK_IMPORTANCEDENSITY_H
#define BTK_IMPORTANCEDENSITY_H

    // Local includes
    #include "btkTypes.h"
    #include "btkPoint.h"
    #include "btkDirection.h"
    #include "btkSignal.h"
    #include "btkSHModel.h"


    namespace btk
    {

    /**
     * @class ImportanceDensity
     * @brief Importance density
     * @author Julien Pontabry
     * @ingroup Tractography
     */
    class ImportanceDensity
    {
        public:
            /**
             * @brief Constructor
             * Importante density need von Mises-Fisher density and Q-ball HS model
             * @param signal Diffusion signal
             * @param model Q-ball HS model of MRI data
             * @param angleThreshold Angle threshold
             */
            ImportanceDensity(Signal *signal, SHModel *model, Real angleThreshold);

            /**
             * @brief Simulate a direction
             * Simulate a direction to follow
             * @param mean Mean direction
             * @param kappa Concentration of vmf distribution
             * @return Simulated direction
             */
            Direction simulate(Direction mean, Real kappa);

            /**
             * @brief Compute density
             * Compute density for a given direction
             * @param vk Direction
             * @param mean Mean direction
             * @param kappa Concentration of vmf distribution
             * @return Value of probability law corresponding to direction
             */
            Real compute(Direction uk, Direction mean, Real kappa);

            /**
             * @brief Compute mean direction
             * Compute mean direction needed for von Mises-Fisher density (with given data)
             * @param xk Point in euclidian space
             * @param ukm1 Last direction
             * @return Mean direction
             */
            Direction computeMeanDirection(Point xk, Direction ukm1);

            /**
             * @brief Compute concentration corresponding to the direction mu
             * The concentration is proportionnal to the mean curvature of the ODF at mu.
             * @param mu Direction
             * @param xk Point in MRI volume
             * @return The importance parameter for vmf
            */
            Real computeConcentration(Direction mu, Point xk);

        private:
            SHModel *m_model;       /**< Q-ball HS model of MRI data */
            Signal  *m_signal;      /**< Diffusion signal */
            Real m_angleThreshold;  /**< Threshold angle of importance function */

            Real m_2PI;     /**< Precomputed constant */
            Real m_log4PI;  /**< Precomputed constant */
    };

    } // namespace btk

#endif /* BTK_IMPORTANCEDENSITY_H */

