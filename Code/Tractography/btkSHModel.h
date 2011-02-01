/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

12 februar 2010
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


#ifndef BTK_SH_MODEL_H
#define BTK_SH_MODEL_H

    // STL includes
    #include "string"
    #include "vector"

    // Local includes
    #include "btkTypes.h"
    #include "btkDirection.h"
    #include "btkPoint.h"

    namespace btk
    {

    class SHModel
    {
        public:
            /**
             * @brief Constructor
             * @param filename Filename of model ITK image
             */
            SHModel(const std::string &filename);

            SHModel(Sequence::Pointer model, std::vector<Direction> *originalDirections);

            /**
             * @brief Destructor
             */
            ~SHModel();

            /**
             * @brief Get signal value at given direction and position
             * The value is interpolated
             * @param u Direction in wich we want the signal value
             * @param p Point in the image where the signal value is requested
             * @return Signal value at this point and direction
             */
            Real signalAt(Direction u, Point p);

            /**
             * @brief Get signal value at given direction and position
             * The value is interpolated
             * @param p Point in the image where the signal value is requested
             * @return Signal value at this point and direction
             */
            Matrix signalAt(Point p);

            /**
             * @brief Get ODF value at given direction and position
             * The value is interpolated
             * @param u Direction in wich we want the ODF value
             * @param p Point in the image where the ODF value is requested
             * @return ODF value at this point and direction
             */
            Real odfAt(Direction u, Point p);

            /**
             * @brief Get ODF values at given position
             * Directions are gradient directions
             * @param p Point in the image where ODF values are requested
             * @return Matrix containing ODF values at that position
             */
            Matrix odfAt(Point p);

            /**
             * @brief Get the direction where ODF is maximal
             * @param p Position in the image
             * @return Direction where at this position, ODF is maximal
             */
            Direction getMaxDirectionAt(Point p);

            /**
             * @brief Get the direction where ODF is maximal
             * @param p Position in the image
             * @return Directions where at these locations, ODF is maximal
             */
            std::vector<Direction> getMaxDirectionsAt(Point p);

            Image::PointType getImageOrigin()
            {
                return m_model[0]->GetOrigin();
            }

            Image::SpacingType getImageSpacing()
            {
                return m_model[0]->GetSpacing();
            }

            Image::SizeType getImageSize()
            {
                return m_model[0]->GetLargestPossibleRegion().GetSize();
            }

			std::vector<Direction> &getDirections()
			{
				return *m_directions;
			}

            itk::Matrix<Real,3,3> GetDirection()
            {
                return m_model[0]->GetDirection();
            }

        private:
            /**
             * @brief Read image file and directions file
             * @param filename Filename of model ITK image
             */
            Sequence::Pointer readFiles(const std::string &filename);

            /**
             * @brief Compute Legendre matrix
             * This matrix is needed for ODF evaluation
             */
            void computeLegendreMatrix();

            /**
             * @brief Commpute SH basis matrix
             */
            void computeSHBasisMatrix();

            void computeSHBasisOriMatrix();

            /**
             * @brief Build coefficient to get sharp ODF
             */
            void buildSharperODFMatrix();

        private:
            Image::Pointer             *m_model;        /**< Data array containing model image */
            ImageInterpolator::Pointer *m_interp;       /**< Interpolators' images */
            std::vector<Direction>     *m_directions;   /**< Directions */
            std::vector<Direction>     *m_originalDirections;

            unsigned int m_order;   /**< Order of model */
            unsigned int m_R;       /**< Number of even spherical harmonics basis */
            unsigned int m_reso;    /**< Resolution of ODF */
            Real m_maxTheta;        /**< Theta's maximal value */
            Real m_maxPhi;          /**< Phi's maximal value */
            Real m_pasTheta;        /**< Theta's pas */
            Real m_pasPhi;          /**< Phi's pas */

            Matrix *m_P;     /**< Legendre matrix */
            Matrix *m_Y;     /**< SH basis matrix */
            Matrix *m_Yori;
            Matrix *m_Sharp; /**< Sharper ODF matrix */

            Real m_4PI;
            Real m_2PI;
    };

    } // namespace btk

#endif // BTK_SH_MODEL_H

