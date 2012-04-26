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

#ifndef BTK_SPHERICALCOORDINATES_H
#define BTK_SPHERICALCOORDINATES_H

    // STL includes
    #include "cmath"

    // Local includes
    #include "btkTypes.h"
    #include "btkCartesianCoordinates.h"


    namespace btk
    {

    class CartesianCoordinates;

    /**
     * @class SphericalCoordinates
     * @brief Spherical coordinates
     * @author Julien Pontabry
     * @ingroup Tractography
     * Spherical coordinates in a Q-ball space (unit sphere)
     */
    class SphericalCoordinates
    {
        public:
            /**
              * @brief Constructor
              * Build new null spherical coordinates
              */
            SphericalCoordinates();

            /**
             * @brief Constructor
             * Build spherical coordinates given two parameters theta and phi
             * @param theta
             * @param phi
             * @param rau
             */
            SphericalCoordinates(Real theta, Real phi, Real rho);

            /**
             * @brief Get theta coordinate
             * @return theta
             */
            Real theta() const;

            /**
             * @brief Get phi coordinate
             * @return phi
             */
            Real phi() const;

            /**
              * @brief Get rau coordinate
              * @return rau
              */
            Real rho() const;

            /**
             * @brief Conversion to cartesian coordinates
             * @return Corresponding cartesian coordinates
             */
            CartesianCoordinates toCartesianCoordinates();

            /**
              * @brief Display spherical coordinates on output stream
              * @param os Output stream
              * @param u Coordinates to display
              */
            friend std::ostream &operator<<(std::ostream &os, const SphericalCoordinates& u);

        protected:
            Real m_theta;   /**< Theta */
            Real m_phi;     /**< Phi */
            Real m_rho;     /**< Rau */
    };

    } // namespace btk

#endif // BTK_SPHERICALCOORDINATES_H

