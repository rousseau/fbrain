/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

22 march 2010
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

#ifndef BTK_CARTESIANCOORDINATES_H
#define BTK_CARTESIANCOORDINATES_H

    // STL includes
    #include "cmath"
    #include "ostream"

    // Local includes
    #include "btkTypes.h"
    #include "btkSphericalCoordinates.h"


    namespace btk
    {

    class SphericalCoordinates;

    /**
     * @class CartesianCoordinates
     * @brief Cartesian coordinates
     * @author Julien Pontabry
     * Cartesian coordinates in a 3D space
     */
    class CartesianCoordinates
    {
        public:
            /**
              * @brief Constructor
              * Build a new null coordinates
              */
            CartesianCoordinates();

            /**
             * @brief Constructor
             * Build new coordinates given 3 parameters
             * @param x Coordinate along x axis
             * @param y Coordinate along y axis
             * @param z Coordinate along z axis
             */
            CartesianCoordinates(Real x, Real y, Real z);

            /**
             * @brief Get x coordinate
             * @return Coordinate along x axis
             */
            Real x() const;

            /**
             * @brief Get y coordinate
             * @return Coordinate along y axis
             */
            Real y() const;

            /**
             * @brief Get z coordinate
             * @return Coordinate along z axis
             */
            Real z() const;

            /**
             * @brief Conversion to spherical coordinates
             * @return Corresponding spherical coordinates
             */
            SphericalCoordinates toSphericalCoordinates();

            /**
              * @brief Multiply coordinates by a factor
              * @return New coordinates
              */
            CartesianCoordinates operator*(Real factor);

            /**
              * @brief Display cartesian coordinates on output stream
              * @param os Output stream
              * @param c Coordinates to display
              */
            friend std::ostream &operator<<(std::ostream &os, const CartesianCoordinates& c);

        private:
            Real m_x; /**< X coordinate */
            Real m_y; /**< Y coordinate */
            Real m_z; /**< Z coordinate */
    };

    } // namespace btk

#endif // BTK_CARTESIANCOORDINATES_H

