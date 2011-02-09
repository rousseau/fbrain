/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

6 september 2010
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


#ifndef BTK_SIGNAL
#define BTK_SIGNAL

    // STL includes
    #include "string"
    #include "vector"

    // Local includes
    #include "btkTypes.h"
    #include "btkPoint.h"
    #include "btkDirection.h"


    namespace btk
    {

    class Signal
    {
        public:
            /**
             * @brief Constructor
             * @param filename Signal filename
             */
            Signal(const std::string &filename, const std::string &sigmasFilename, const std::string &dirFileName);

            Signal(Sequence::Pointer signal, std::vector<Real> *sigmas, std::vector<Direction> *directions, char displayMode);

            /**
             * @brief Destructor
             */
            ~Signal();

            /**
             * @brief Get signal at specified position
             * @param p Position
             * @return Signal at specified position p
             */
            Matrix signalAt(Point p);

            /**
             * @brief Get signal standard deviations
             * @return Vector of standard deviations
             */
            std::vector<Real> *getSigmas();

            /**
             * @brief Get image's sizes
             * @return Size of the 3D image
             */
            Image::SizeType getSize();

            /**
             * @brief Get image's origin
             * @return Origin of the 3D image
             */
            Image::PointType getOrigin();

            /**
             * @brief Get image's spacing
             * @return Spacing of the 3D image
             */
            Image::SpacingType getSpacing();

            /**
             * @brief Get gradient directions
             * @return Pointer to vector of gradient directions
             */
            std::vector<Direction> *getDirections();

        private:
            /**
              * @brief Read some files
              * Files will be red
              * @param filename Signal image's filename
              * @param sigmasFilename sigmas' filename
              */
            Sequence::Pointer readFiles(const std::string &filename, const std::string &sigmasFilename, const std::string &dirFileName);

        private:
            Image::Pointer             *m_signal;       /**< Data array containing signal image */
            ImageInterpolator::Pointer *m_interp;       /**< Image's interpolators */
            std::vector<Direction>     *m_directions;   /**< Gradient directions */

            std::vector<Real> *m_sigmas;

            unsigned int m_N;   /**< Number of images */

            char m_displayMode;
    };

    } // namespace btk

#endif // BTK_SIGNAL

