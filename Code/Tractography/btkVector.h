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

#ifndef BTK_VECTOR_H
#define BTK_VECTOR_H

    // Local includes
    #include "btkTypes.h"
    #include "btkDirection.h"


    namespace btk
    {

        class Direction;

        /**
         * @class Vector
         * @brief Vector in 3D space
         * @author Julien Pontabry
         * @ingroup Tractography
         */
        class Vector : public CartesianCoordinates
        {
            public:
                /**
                 * @brief Constructor
                 * Build a new null vector
                 */
                Vector();

                /**
                 * @brief Constructor
                 * Build a new vector with the given coordinates
                 * @param x X-axis coordinate
                 * @param y Y-axis coordinate
                 * @param z Z-axis coordinate
                 */
                Vector(Real x, Real y, Real z);

                /**
                 * @brief Verify if vector is the null vector
                 * @return True if the vector is the null vector, false otherwise
                 */
                bool isNull();

                /**
                 * @brief Assuming this is an unit vector, convert to a direction
                 * @return Direction corresponding to this unit vector
                 */
                Direction toDirection();

                /**
                  * @brief Multiply coordinates by a factor
                  * @param factor Multiplication factor
                  * @return New coordinates
                  */
                Vector operator*(Real factor);

                /**
                 * @brief Multiply coordinates by a factor, components by components
                 * @param factorX Multiplication factor for X-axis
                 * @param factorY Multiplication factor for Y-axis
                 * @param factorZ Multiplication factor for Z-axis
                 * @return New coordinates
                 */
                 Vector operator*(std::vector<Real> factors);

                /**
                 * @brief Gives the angle between two vectors in radian
                 * @param v Vector we want to know the angle with this
                 * @return An angle measure in radian
                 */
                Real angleWith(Vector v);

                /**
                 * @brief Normalize the vector
                 */
                 void Normalize();
        };

    } // namespace btk

#endif // BTK_VECTOR_H

