/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 februar 2010
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


#ifndef BTK_SH_MODEL_ESTIMATOR_H
#define BTK_SH_MODEL_ESTIMATOR_H

    // STL includes
    #include "string"
    #include "vector"

    // Local includes
    #include "btkTypes.h"
    #include "btkDirection.h"


    namespace btk
    {

    /**
     * @class SHModelEstimator
     * @brief Estimate diffusion model
     * @author Julien Pontabry
     * @ingroup Tractography
     */
    class SHModelEstimator
    {
        public:
            /**
             * @brief Constructor
             * @param signalFileName Signal sequence's filename
             * @param directionsFileName Directions' filename
             * @param maskFileName Image mask filename
             * @param order Model's order
             * @param lambda Regularization term
             */
            SHModelEstimator(const std::string &signalFileName, const std::string &directionsFileName, const std::string &maskFileName, unsigned int order, Real lambda, char displayMode);

            /**
             * @brief Constructor
             * @param signal Diffusion signal image
             * @param directions Gradient directions
             * @param mask Image Mask image
             * @param order Model's order
             * @param lambda Regularization term
             */
            SHModelEstimator(Sequence::Pointer signal, std::vector<Direction> *directions, Mask::Pointer mask, unsigned int order, Real lambda, char displayMode);

            /**
             * @brief Destructor
             */
            ~SHModelEstimator();

            /**
             * @brief Estimate model
             */
            void estimate();

            /**
             * @brief Save model
             * Save model sequence in nifti file
             */
            void save();

            /**
             * @brief Get the model image representation
             * @return The model coefficients
             */
            Sequence::Pointer GetModel();

        private:
            /**
             * @brief Read needed files
             * @param signalFileName Signal sequence's filename
             * @param directionsFileName Directions' filename
             * @param maskFileName Image mask filename
             */
            void readFiles(const std::string &signalFileName, const std::string &directionsFileName, const std::string &maskFileName);

            /**
             * @brief Commpute SH basis matrix
             */
            void computeSHBasisMatrix();

            /**
             * @brief Compute Laplace-Beltrami smoothing matrix
             */
            void computeLaplaceBeltramiMatrix();

            /**
             * @brief Compute Omega matrix
             */
            void computeOmegaMatrix();

        private:
            std::vector<Direction> *m_directions;   /**< Gradient's directions */
            Sequence::Pointer       m_signal;       /**< Signal sequence */
            Sequence::Pointer       m_model;        /**< Model sequence */
            Mask::Pointer           m_mask;         /**< Image mask */

            unsigned int m_order;   /**< Max order of SH decomposition */
            unsigned int m_R;       /**< Number of even SH basis */
            Real         m_lambda;  /**< Regularization parameter */

            Real m_2PI; /**< Precomputed constant */

            Matrix *m_Y;       /**< SH basis matrix */
            Matrix *m_L;       /**< Laplace-Beltrami smoothing matrix */
            Matrix *m_P;       /**< Legendre matrix */
            Matrix *m_Omega;   /**< Omega matrix */

            char m_displayMode;
    };

    } // namespace btk

#endif // BTK_SH_MODEL_ESTIMATOR_H

