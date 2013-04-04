/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 04/01/2013
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


#ifndef BTK_PARTICLE_FILTERING_TRACTOGRAPHY_ALGORITHM_H
#define BTK_PARTICLE_FILTERING_TRACTOGRAPHY_ALGORITHM_H

// ITK includes
#include "itkMacro.h"
#include "itkSmartPointer.h"

// Local includes
#include "btkMacro.h"
#include "btkTractographyAlgorithm.h"
#include "tmp/btkParticle.h"
#include "tmp/btkImportanceDensity.h"
#include "tmp/btkLikelihoodDensity.h"
#include "tmp/btkPriorDensity.h"

namespace btk
{

/**
 * @brief Define a probabilistic tractography algorithm within a particle filtering framework.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class ParticleFilteringTractographyAlgorithm : public TractographyAlgorithm
{
    public:
        typedef ParticleFilteringTractographyAlgorithm Self;
        typedef btk::TractographyAlgorithm             Superclass;
        typedef itk::SmartPointer< Self >              Pointer;
        typedef itk::SmartPointer< const Self >        ConstPointer;

        typedef Superclass::PhysicalPoint PhysicalPoint;
        typedef itk::Image< double,3 >    ProbabilityMap;

        itkNewMacro(Self);

        itkTypeMacro(ParticleFilteringTractographyAlgorithm, TractographyAlgorithm);

        btkSetMacro(NumberOfParticles, unsigned int);
        btkGetMacro(NumberOfParticles, unsigned int);

        btkSetMacro(ParticleStepSize, double);
        btkGetMacro(ParticleStepSize, double);

        btkSetMacro(ResamplingThreshold, double);
        btkGetMacro(ResamplingThreshold, double);

        btkSetMacro(CurveConstraint, double);
        btkGetMacro(CurveConstraint, double);

        btkSetMacro(ThresholdAngle, double);
        btkGetMacro(ThresholdAngle, double);

    protected:
        /**
         * @brief Constructor.
         */
        ParticleFilteringTractographyAlgorithm();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /**
         * @brief Initialize the filter.
         */
        virtual void Initialize();

        /**
         * @brief Propagate using the tractography algorithm at a seed point.
         * @param point Seed point.
         */
        virtual vtkSmartPointer< vtkPolyData > PropagateSeed(Self::PhysicalPoint point);

    private:
        /**
         * @brief Propagate a seed using a particle filter.
         * @param points Vector of points initilized with the coordinates of the seed.
         * @param nextDirection First direction of propagation.
         */
        void PropagateSeed(std::vector< Self::PhysicalPoint > &points, btk::GradientDirection nextDirection);

        /**
         * @brief Resampling function of the particles' cloud.
         * @param cloud Cloud of particles to resample.
         * @param numberOfActiveParticles Number of currently active particles
         */
        void ResampleParticlesCloud(std::vector< Particle > &cloud, unsigned int numberOfActiveParticles);

        /**
         * @brief Save the particles' cloud in file (filename si 'cloud-[world coordinates]-pathlength.vtk').
         * @param cloud Particles' cloud to save in file.
         */
        void SaveCloud(std::vector< Particle > &cloud);

        /**
         * @brief Compute a probability map from a cloud of particles.
         * @param cloud Cloud of particles defining the probability.
         */
        void ComputeProbabilityMap(std::vector< Particle > &cloud);

        /**
         * @brief Computethe maximum a posteriori of the particles' cloud.
         * @param map Points of maximum a posteriori path.
         * @param cloud Cloud of particles.
         * @param numberOfIterations Number of iterations of the particles' cloud.
         */
        void ComputeMaximumAPosteriori(std::vector< Self::PhysicalPoint > &map, std::vector< Particle > &cloud, unsigned int numberOfIterations);

    private:
        /**
         * @brief Number of particles used for the particle filter.
         */
        unsigned int m_NumberOfParticles;

        /**
         * @brief Displacement step of a particle (in mm).
         */
        double m_ParticleStepSize;

        /**
         * @brief Initial weight of the particles.
         */
        double m_InitialWeight;

        /**
         * @brief Resampling threshold of the process.
         */
        double m_ResamplingThreshold;

        /**
         * @brief Curve constraint (concentration of the prior von Mises-Fisher density).
         */
        double m_CurveConstraint;

        /**
         * @brief Angle threshold (the next mean direction of importance density is searched in the solid angle around the previous mean).
         */
        double m_ThresholdAngle;

        /**
         * @brief Importance probability density.
         */
        ImportanceDensity m_ImportanceDensity;

        /**
         * @brief Likelihood probability density.
         */
        LikelihoodDensity m_LikelihoodDensity;

        /**
         * @brief Prior probability density.
         */
        PriorDensity m_PriorDensity;

        /**
         * @brief Probability map of clouds of particles.
         */
        ProbabilityMap::Pointer m_ProbabilityMap;
};

} // namespace btk

#endif // BTK_PARTICLE_FILTERING_TRACTOGRAPHY_ALGORITHM_H
