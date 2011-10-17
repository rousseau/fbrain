/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

31 march 2010
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

#ifndef BTK_PARTICLEFILTER_H
#define BTK_PARTICLEFILTER_H

    // STL includes
    #include "vector"
    #include "string"

    // VTK includes
    #include "vtkSmartPointer.h"
    #include "vtkPolyData.h"

    // Local includes
    #include "btkTypes.h"
    #include "btkPoint.h"
    #include "btkParticle.h"
    #include "btkAPrioriDensity.h"
    #include "btkLikelihoodDensity.h"
    #include "btkImportanceDensity.h"
    #include "btkSHModel.h"
    #include "btkSignal.h"


    namespace btk
    {

        /**
         * @class ParticleFilter
         * @brief Particle filter
         * @author Julien Pontabry
         * @ingroup Tractography
         * Particle filter to approach the probability density presence of fibers in image
         */
        class ParticleFilter
        {
            public:
                /**
                 * @brief Constructor
                 * @param model Diffusion model
                 * @param aPriori A priori density
                 * @param likelihood Likelihood density
                 * @param importance Importance density
                 * @param maskFileName Image mask's filename
                 * @param size Size of the image
                 * @param origin Origin of the image
                 * @param spacing Spacing of the image
                 * @param M Amount of particle
                 * @param x0 Seed where all particles starts
                 * @param epsilon Resampling threshold
                 * @param stepSize Size of each step (norm of vectors)
                 * @param maxLength Maximal length of particles
                 */
                ParticleFilter(Signal *signal, SHModel *model, APrioriDensity aPriori, LikelihoodDensity likelihood, ImportanceDensity importance,
                               const std::string maskFileName, Image::SizeType size, Image::PointType origin, Image::SpacingType spacing,
                               unsigned int M, Point x0, Real epsilon, Real stepSize, unsigned int maxLength);

                /**
                 * @brief Constructor
                 * @param model Diffusion model
                 * @param aPriori A priori density
                 * @param likelihood Likelihood density
                 * @param importance Importance density
                 * @param mask Mask image
                 * @param size Size of the image
                 * @param origin Origin of the image
                 * @param spacing Spacing of the image
                 * @param M Amount of particle
                 * @param x0 Seed where all particles starts
                 * @param epsilon Resampling threshold
                 * @param stepSize Size of each step (norm of vectors)
                 */
                ParticleFilter(Signal *signal, SHModel *model, APrioriDensity aPriori, LikelihoodDensity likelihood, ImportanceDensity importance,
                               Mask::Pointer mask, Image::SizeType size, Image::PointType origin, Image::SpacingType spacing,
                               unsigned int M, Point x0, Real epsilon, Real stepSize, char displaMode);

                /**
                 * @brief Launch particle filter
                 * @param label Label of seed
                 */
                void run(int label);

                /**
                 * @brief Save cloud in VTK file
                 * @param label Label number (of label image)
                 * @param step Step number
                 */
                void saveCloudInVTK(int label, unsigned int step, Point begin);


                /**
                 * @brief Save fiber in VTK file
                 * @param label Label number (of label image)
                 * @param step Step number
                 */
                void saveFiber(int label/*, unsigned int step, Point begin*/);

                /**
                 * @brief Save Connection map in 3D image
                 * @param label Label number (of label image)
                 * @param begin Starting point
                 */
                void saveConnectionMap(int label, Point begin);

                /**
                 * @brief Get the computed connection map
                 * @return Connection map
                 */
                Image::Pointer GetConnectionMap();

                /**
                 * @brief Get the computed fibers
                 * @return Fibers' vtk data structure
                 */
                vtkSmartPointer<vtkPolyData> GetFiber();

                /**
                 * @brief Set on coordinates expressed in LPS (Left-Posterior-Superior)
                 */
                void SetLPSOn();


            private:
                /**
                 * @brief Run the filter
                 * @param label Label number
                 * @param dir Current starting direction
                 */
                void run(int label, Direction dir);

                /**
                 * @brief Get the MAP of the current cloud
                 * @return The particle with the maximal importance weight
                 */
                Particle GetMAP();

                /**
                 * @brief Compute the connection map
                 */
                void ComputeMap();

                /**
                 * @brief Compute the fiber given two MAP
                 * @param map1 MAP in the first direction
                 * @param map2 MAP in the second direction
                 */
                void ComputeFiber(Particle map1, Particle map2);

                /**
                 * @brief Resample all particles
                 * @param nbInsPart Number of particles inside mask
                 */
                 void ResampleCloud(unsigned int nbInsPart);


                double ComputeGFA(Point p);

            private:
                APrioriDensity    m_aPriori;       /**< A priori density */
                LikelihoodDensity m_likelihood;    /**< Likelihood density */
                ImportanceDensity m_importance;    /**< Importance density */

                Mask::Pointer m_mask;   /**< Image mask */

                Point        m_x0;          /**< Start point */
                unsigned int m_k;           /**< Current step */
                unsigned int m_M;           /**< Number of particles */
                Real         m_epsilon;     /**< Resampling threshold */
                Real         m_stepSize;    /**< Size of steps */
                Real         m_kx;          /**< Step size for x axis */
                Real         m_ky;          /**< Step size for y axis */
                Real         m_kz;          /**< Step size for z axis */
                unsigned int m_maxLength;   /**< Maximal length of particles */

                Image::PointType   m_origin;
                Image::SpacingType m_spacing;

                std::vector<Particle> m_cloud;  /**< Cloud of particles */

                Image::Pointer m_map;                   /**< Density map */
                vtkSmartPointer<vtkPolyData> m_fiber;   /**< Fibers' vtk structure */

                SHModel *m_model;   /**< Diffusion model */
                Signal  *m_signal;

                unsigned int m_dirNum; /**< Current direction number */

                char m_displayMode;

                std::vector<Real> m_vStepSize;

                unsigned int m_maxStep;
                unsigned int m_maxSize;

                bool m_lps;
        };

    } // namespace btk

#endif // BTK_PARTICLEFILTER_H

