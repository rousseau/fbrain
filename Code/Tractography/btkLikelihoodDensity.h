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

#ifndef BTK_LIKELIHOOD_DENSITY_H
#define BTK_LIKELIHOOD_DENSITY_H

// ITK incluudes
#include "itkPoint.h"

// Local includes
#include "btkGradientDirection.h"
#include "btkDiffusionModel.h"
#include "btkDiffusionSignal.h"

namespace btk
{

/**
 * @brief Likelihood density for particle filtering.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class LikelihoodDensity
{
    public:
        typedef itk::Point< double,3 > PhysicalPoint;

        /**
         * @brief Constructor.
         */
        LikelihoodDensity();

        /**
         * @brief Constructor.
         * @param signal Diffusion signal.
         * @param model Diffusion model.
         */
        LikelihoodDensity(DiffusionSignal::Pointer signal, DiffusionModel::Pointer model);

        /**
         * @brief Evaluate the likelihood probability density.
         * @param vk Direction vector to evaluate.
         * @param pk Location in the volume.
         * @param mu Mean direction for simulation.
         * @return The value of the likelihood function at pk, in direction vk with a given mean.
         */
        double Evaluate(GradientDirection vk, PhysicalPoint xk, GradientDirection mu);

    private:
        /**
         * @brief Evaluate the logarithm of the normal centered probability density.
         * @param sigma Standard deviation of the probability.
         * @param x Point of evaluation
         * @return The value of the distribution in x.
         */
        inline double EvaluateNormalCenteredLogDensity(double sigma, double x);

    private:
        /**
         * @brief Diffusion model.
         */
        DiffusionModel::Pointer m_Model;

        /**
         * @brief Diffusion signal.
         */
        DiffusionSignal::Pointer m_Signal;

        /**
         * @brief Precomputed constant.
         */
        static const double m_logSqrt2Pi = 0.918938533204673;
};

} // namespace btk

#endif // BTK_LIKELIHOOD_DENSITY_H
