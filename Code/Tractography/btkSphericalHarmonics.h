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


#ifndef BTK_SPHERICAL_HARMONICS_H
#define BTK_SPHERICAL_HARMONICS_H

    // Local includes
    #include "btkTypes.h"
    #include "btkDirection.h"


    namespace btk
    {

    /**
     * @class SphericalHarmonics
     * @brief Spherical harmonics mathematics implementation
     * @author Julien Pontabry
     * @ingroup Tractography
     */
    class SphericalHarmonics
    {
        public:
            /**
             * @brief Compute SH basis
             * Compute spherical harmonic basis in direction, order and degree given
             * @param u Direction of basis
             * @param l Order of spherical harmonic
             * @param m Degree of spherical harmonic
             * @return Spherical harmonic basis' coefficient
             */
            static Real computeBasis(Direction u, unsigned int l, int m);

            /**
             * @brief Legendre polynom in zero
             * Compute Legendre polynom in zero at given order
             * @param l Order
             * @return Value of Legendre polynom at order l for parameter x
             */
            static Real legendrePolynomialInZero(unsigned int l);

    #ifndef NDEBUG
        public:
    #else
        private:
    #endif // NDEBUG
            /**
             * @brief Factorial function
             * @param n Number to compute the factorial
             * @return Factorial of n
             */
            static unsigned int factorial(unsigned int n);

            /**
             * @brief Associated Legendre polynom
             * Compute associated Legendre polynom with given order and degree
             * @param l Order
             * @param m Degree
             * @param theta Parameters of Legendre polynomial (x = cos(theta))
             * @return Value of associated Legendre polynomial for cos(theta) at order l and degree m
             */
            static Real legendrePolynomial(unsigned int l, unsigned int m, Real theta);
    };

    } // namespace btk

#endif // BTK_SPHERICAL_HARMONICS_H

