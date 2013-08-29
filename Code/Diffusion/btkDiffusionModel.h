/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 09/07/2012
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

#ifndef BTK_DIFFUSION_MODEL_H
#define BTK_DIFFUSION_MODEL_H

// STL includes
#include "vector"

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkContinuousIndex.h"
#include "itkPoint.h"

// Local includes
#include "btkGradientDirection.h"

namespace btk
{

/**
 * @brief Diffusion MRI modelization.
 * @author Julien Pontabry
 * @ingroup Diffusion
 */
class DiffusionModel : public itk::ProcessObject
{
    public:
        typedef DiffusionModel                  Self;
        typedef itk::ProcessObject              Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef itk::ContinuousIndex< float,3 > ContinuousIndex;
        typedef itk::Point< float,3 >           PhysicalPoint;

        itkTypeMacro(DiffusionModel,itk::ProcessObject);

        btkGetMacro(Directions, std::vector< btk::GradientDirection >);

        btkSetMacro(SphericalResolution, unsigned int);
        btkGetMacro(SphericalResolution, unsigned int);

        /**
         * @brief Get modeling at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Model response in direction direction at cindex in image space.
         */
        virtual float ModelAt(ContinuousIndex cindex, btk::GradientDirection direction) = 0;

        /**
         * @brief Get modeling at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param directiosn Gradient directions were the model response is wanted.
         * @return Model response in direction direction at cindex in image space.
         */
        virtual std::vector< float > ModelAt(ContinuousIndex cindex, std::vector< btk::GradientDirection > &directions) = 0;

        /**
         * @brief Get modeling at continuous index.
         * @param cindex Location in the image space.
         * @return Model response at cindex in image space.
         */
        virtual std::vector< float > ModelAt(ContinuousIndex cindex) = 0;

        /**
         * @brief Get modeling at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Model response in direction direction at point in physical space.
         */
        virtual float ModelAt(PhysicalPoint point, btk::GradientDirection direction) = 0;

        /**
         * @brief Get modeling at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Model response in direction direction at point in physical space.
         */
        virtual std::vector< float > ModelAt(PhysicalPoint point, std::vector< btk::GradientDirection > &directions) = 0;

        /**
         * @brief Get modeling at physical point.
         * @param point Point in the physical space.
         * @return Model response at point in physical space.
         */
        virtual std::vector< float > ModelAt(PhysicalPoint point) = 0;

        /**
         * @brief Get signal at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Signal response in direction direction at cindex in image space.
         */
        virtual float SignalAt(ContinuousIndex cindex, btk::GradientDirection direction) = 0;

        /**
         * @brief Get signal at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Signal response in direction direction at cindex in image space.
         */
        virtual std::vector< float > SignalAt(ContinuousIndex cindex, std::vector< btk::GradientDirection > &directions) = 0;

        /**
         * @brief Get signal at continuous index.
         * @param cindex Location in the image space.
         * @return Signal response at cindex in image space.
         */
        virtual std::vector< float > SignalAt(ContinuousIndex cindex) = 0;

        /**
         * @brief Get signal at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Signal response in direction direction at point in physical space.
         */
        virtual float SignalAt(PhysicalPoint point, btk::GradientDirection direction) = 0;

        /**
         * @brief Get signal at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Signal response in direction direction at point in physical space.
         */
        virtual std::vector< float > SignalAt(PhysicalPoint point, std::vector< GradientDirection > &directions) = 0;

        /**
         * @brief Get signal at physical point.
         * @param point Point in the physical space.
         * @return Signal response at point in physical space.
         */
        virtual std::vector< float > SignalAt(PhysicalPoint point) = 0;

        /**
         * @brief Get mean directions at a location in the physical space.
         * @param point Point in the physical space.
         * @return Vector of mean directions of local model at a physical location point.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(PhysicalPoint point) = 0;

        /**
         * @brief Get mean directions in a solid angle formed by previous vector and search angle at a location in the physical space.
         * @param point Point in the physical space.
         * @param vector Previous vector
         * @param angle Angle of search.
         * @return Vector of mean directions of local model in a solid angle formed by previous vector and search angle at a physical location point.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(PhysicalPoint point, btk::GradientDirection vector, float angle) = 0;

        /**
         * @brief Get mean directions at a location in the image space.
         * @param cindex Continuous index in the image space.
         * @return Vector of mean directions of local model at a location in image space.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(ContinuousIndex cindex) = 0;

        /**
         * @brief Get mean directions in a solid angle formed by previous vector and search angle at a location in the physical space.
         * @param cindex Continuous index in the image space.
         * @param vector Previous vector
         * @param angle Angle of search.
         * @return Vector of mean directions of local model in a solid angle formed by previous vector and search angle at a physical location point.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(ContinuousIndex cindex, btk::GradientDirection vector, float angle) = 0;

        /**
         * @brief Update the parameters of the function (precompute some variables).
         */
        virtual void Update();

    protected:
        /**
         * @brief Constructor.
         */
        DiffusionModel();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /**
         * @brief Sampled directions on the unit sphere used by modeling reconstruction.
         */
        std::vector< btk::GradientDirection > m_Directions;

        /**
         * @brief Spherical resolution of the modeling reconstruction (number of point on the unit sphere).
         */
        unsigned int m_SphericalResolution;
};

} // namespace btk

#endif // BTK_DIFFUSION_MODEL_H
