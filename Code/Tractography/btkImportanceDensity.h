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

#ifndef BTK_IMPORTANCE_DENSITY_H
#define BTK_IMPORTANCE_DENSITY_H

// ITK includes
#include "itkPoint.h"

// Local includes
#include "btkMacro.h"
#include "btkGradientDirection.h"
#include "btkDiffusionModel.h"
#include "btkDiffusionSignal.h"

namespace btk
{

/**
 * @brief Importance probability density of particle filtering.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class ImportanceDensity
{
    public:
        typedef itk::Point< double,3 > PhysicalPoint;

        btkGetMacro(Model, DiffusionModel::Pointer);
        btkSetMacro(Model, DiffusionModel::Pointer);

        btkGetMacro(Signal, DiffusionSignal::Pointer);
        btkSetMacro(Signal, DiffusionSignal::Pointer);

        /**
         * @brief Constructor.
         */
        ImportanceDensity();

        /**
         * @brief Constructor.
         * @param model Diffusion model used.
         * @param signal Diffusion signal.
         * @param angleThreshold Threshold angle of the density.
         */
        ImportanceDensity(DiffusionModel::Pointer model, DiffusionSignal::Pointer signal, double angleThreshold);

        /**
         * @brief Compute the mean direction
         * @param pk Current position.
         * @param vkm1 Last direction.
         * @return The mean direction.
         */
        GradientDirection GetMeanDirection(PhysicalPoint pk, GradientDirection vkm1);

        /**
         * @brief Simulate the importance density.
         * @param meanDirection Mean direction of the simulation.
         * @param concentration Concentration parameter.
         * @return Simulated direction.
         */
        GradientDirection Simulate(GradientDirection meanDirection, double concentration);

        /**
         * @brief Evaluate the importance density.
         * @param vk Direction to evaluate.
         * @param vkm1 Previous direction
         * @param concentration Concentration parameter of the von Mises-Fisher density.
         * @return The value of the density probability for vk with given concentration an mean direction.
         */
        double Evaluate(GradientDirection vk, GradientDirection vkm1, double concentration);

        /**
         * @brief Estimate the concentration parameter.
         * @param mu Mean direction.
         * @param pk Position in the image.
         * @return Estimated concentration parameter.
         */
        double EstimateConcentrationParameter(GradientDirection mu, PhysicalPoint pk);

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
         * @brief Threshold angle used for propagation.
         */
        double m_AngleThreshold;

        /**
         * @brief Precomputed constant.
         */
        static const double m_2pi = 6.283185307179586;
};

} // namespace btk

#endif // BTK_IMPORTANCE_DENSITY_H
