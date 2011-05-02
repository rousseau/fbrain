/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

14 april 2010
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

#ifndef BTK_DIRECTION_H
#define BTK_DIRECTION_H

    // Local includes
    #include "btkTypes.h"
    #include "btkSphericalCoordinates.h"
    #include "btkVector.h"


    namespace btk
    {

    class Vector;

    /**
     * @class Direction
     * @brief Direction in 3D space
     * @author Julien Pontabry
     */
    class Direction : public SphericalCoordinates
    {
        public:
            /**
             * @brief Constructor
             * Build a new null direction
             */
            Direction();

            /**
             * @brief Constructor
             * Build a new direction with the given parameters
             * @param azimuth Azimuth coordinate
             * @param elevation Elevation coordinate
             */
            Direction(Real azimuth, Real elevation);

            /**
             * @brief Conversion to vector
             * Convert the direction un a cartesian vector in 3D space
             * @return A cartesian vector in 3D space
             */
            Vector toVector();

            /**
             * @brief Set the direction to be null
             */
            void setNull();

            /**
             * @brief Get if the direction is null
             * @return True if the direction is null, false otherwise
             */
            bool isNull();

            /**
             * @brief Display the vector on stream outputs
             * @param os Current output stream
             * @param u Direction to display
             * @return An output stream with the direction displayed on
             */
            friend std::ostream &operator<<(std::ostream &os, const Direction& u);
    };

    } // namespace btk

#endif // BTK_DIRECTION_H
