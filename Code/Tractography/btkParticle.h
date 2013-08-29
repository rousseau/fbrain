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

#ifndef BTK_PARTICLE_H
#define BTK_PARTICLE_H

// STL includes
#include "vector"

// ITK includes
#include "itkPoint.h"

// Local includes
#include "btkGradientDirection.h"

namespace btk
{

/**
 * @brief Particle of particle filtering.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class Particle
{
    public:
        typedef Particle         Self;
        typedef itk::LightObject Superclass;

        typedef itk::Point< double,3 > PhysicalPoint;

        /**
         * @brief Constructor.
         */
        Particle();

        /**
         * @brief Copy constructor.
         * @param p Particle to copy.
         */
        Particle(const Self &p);

        /**
         * @brief Contructor with initialization to first point.
         * @param p First point of the particle's path.
         */
        Particle(PhysicalPoint p);

        /**
         * @brief Get the last point of the particle's path.
         * @return Last point of the particle's path.
         */
        PhysicalPoint GetLastPoint() const;

        /**
         * @brief Get the last direction of the particle's path.
         * @return Last direction of the particle's path.
         */
        GradientDirection GetLastDirection() const;

        /**
         * @brief Get the last weight of the particle.
         * @return Last weight of the particle.
         */
        double GetLastWeight() const;

        /**
         * @brief Test if particle is active (if propagation has not been stopped).
         * @return True if the particle is active, false otherwise.
         */
        bool IsActive() const;

        /**
         * @brief Add displacement to path.
         * @param v Displacement vector at this step.
         * @param p Position of the particle after displacement.
         * @param w Weight of the the particle at this step.
         */
        void AddToPath(GradientDirection v, PhysicalPoint p, double w);

        /**
         * @brief Desactivate the particle (stop its propagation).
         */
        void Desactivate();

        /**
         * @brief Get the position of the particle at step i.
         * @param i Step of the particle
         * @return Point representing position of the particle at given step.
         */
        PhysicalPoint GetPointAtStep(unsigned int i) const;

        /**
         * @brief Get the weight of the particle at step i.
         * @param i Step of the particle.
         * @return Weight of the particle at the given step.
         */
        double GetWeightAtStep(unsigned int i) const;

        /**
         * @brief Get the length of the path.
         * @return The number of points of the particle's path.
         */
        unsigned int GetPathLength() const;

        /**
         * @brief Get displacement vector of the particle at step i.
         * @param i Step of the particle.
         * @return Displacement vector of the particle at given step.
         */
        GradientDirection GetVectorAtStep(unsigned int i) const;

        /**
         * @brief Normalize the last weight with a normalization coefficient.
         * @param normalizationCoefficient Normalization coefficient.
         */
        void NormalizeLastWeightWith(double normalizationCoefficient);

        /**
         * @brief Resample the particle at the specified point and with the specified weight.
         * @param p New position of the particle.
         * @param w New weight of the particle.
         */
        void Resample(PhysicalPoint p, double w);

        /**
         * @brief Add the likelihood to the current step of particle's path.
         * @param likelihood Likelihood density value.
         */
        void AddLikelihood(double likelihood);

        /**
         * @brief Get the likelihood log at given step.
         * @param i Step of the particle
         * @return The likelihood value at given step.
         */
        double GetLikelihoodAtStep(unsigned int i) const;

    private:
        /**
         * @brief Points of the particle's path.
         */
        std::vector< PhysicalPoint > m_Points;

        /**
         * @brief Displacements vectors of the particle's path.
         */
        std::vector< GradientDirection > m_Directions;

        /**
         * @brief Weights of the particle's path.
         */
        std::vector< double > m_Weights;

        /**
         * @brief True if the particle is active, false otherwise.
         */
        bool m_IsActive;

        /**
         * @brief Likelihood log of particle.
         */
        std::vector< double > m_LikelihoodLog;
};

} // namespace btk

#endif // BTK_PARTICLE_H
