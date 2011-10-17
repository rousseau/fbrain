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

#ifndef BTK_POINT_H
#define BTK_POINT_H

    // Local includes
    #include "btkTypes.h"
    #include "btkCartesianCoordinates.h"
    #include "btkVector.h"


    namespace btk
    {

    /**
     * @class Point
     * @brief Point in 3D space
     * @author Julien Pontabry
     * @ingroup Tractography
     */
    class Point : public CartesianCoordinates
    {
        public:
            /**
             * @brief Constructor
             * Build a new point at the origin
             */
            Point();

            /**
             * @brief Constructor
             * Build a new point at the coordinates given in parameters
             * @param x X-axis coordinate
             * @param y Y-axis coordinate
             * @param z Z-axis coordinate
             */
            Point(Real x, Real y, Real z);

            /**
             * @brief Move a point
             * Translate the point with the translation vector v
             * @param v Translation vector
             * @return The translated point in the 3D space
             */
            Point operator+(Vector v);

            /**
             * @brief Verify if this point is equal to the point in paramaters
             * @param p Point to compare
             * @return True if the two point are equal, false otherwise
             */
             bool operator==(Point p);
    };

    } // namespace btk

#endif // BTK_POINT_H

