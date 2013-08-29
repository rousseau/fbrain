/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 07/03/2013
  Author(s): Julien Pontabry (pontabry@unistra.fr)

  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.

==========================================================================*/

#ifndef BTK_PRIOR_DENSITY_H
#define BTK_PRIOR_DENSITY_H

// Local includes
#include "btkMacro.h"
#include "btkGradientDirection.h"

namespace btk
{

/**
 * @brief Prior probability density for particle filtering.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class PriorDensity
{
    public:
        /**
         * @brief Constructor.
         */
        PriorDensity();

        /**
         * @brief Constructor.
         * @param concentration Concentration parameter.
         */
        PriorDensity(double concentration);

        btkGetMacro(Concentration, double);

        /**
         * @brief Evaluate the prior density of a direction depending of a previous vector.
         * @param vk Current vector.
         * @param vkm1 Previous vector.
         * @return The value of the prior density.
         */
        double Evaluate(GradientDirection vk, GradientDirection vkm1);

        /**
         * @brief Simulate the prior density with the mean direction
         * @param mean Mean direction.
         * @return The simulated direction.
         */
        GradientDirection Simulate(GradientDirection vkm1);

    private:
        /**
         * @brief Initialize the density.
         */
        void Initialize();

    private:
        /**
         * @brief Concentration parameter.
         */
        double m_Concentration;

        /**
         * @brief Precomputed constant.
         */
        static const double m_2pi = 6.283185307179586;

        /**
         * @brief Precomputed normalization constant.
         */
        double m_NormalizationCoefficient;

        /**
         * @brief Precomputed constant.
         */
        double m_SimulationTerm;

        /**
         * @brief Precomputed constant.
         */
        double m_SimulationCoefficient;

        /**
         * @brief Precomputed constant.
         */
        double m_InverseConcentration;
};

} // namespace btk

#endif // BTK_PRIOR_DENSITY_H
