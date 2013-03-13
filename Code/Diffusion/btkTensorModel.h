/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 13/08/2012
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

#ifndef BTK_TENSOR_MODEL_H
#define BTK_TENSOR_MODEL_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLinearInterpolateImageFunction.h"

// Local includes
#include "btkDiffusionModel.h"
#include "btkDiffusionTensorReconstructionFilter.h"

namespace btk
{

/**
 * @brief Tensor diffusion model computed by least square fitting.
 * @author Julien Pontabry
 * @ingroup Diffusion
 */
class TensorModel : public btk::DiffusionModel
{
    public:
        typedef TensorModel                     Self;
        typedef btk::DiffusionModel             Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef btk::DiffusionTensorReconstructionFilter::OutputImageType ModelImage;
        typedef itk::LinearInterpolateImageFunction< ModelImage,float >   InterpolateModelFunction;
        typedef Superclass::PhysicalPoint                                 PhysicalPoint;
        typedef Superclass::ContinuousIndex                               ContinuousIndex;

        itkNewMacro(Self);

        itkTypeMacro(TensorModel,btk::DiffusionModel);

        btkSetMacro(BValue, unsigned int);
        btkGetMacro(BValue, unsigned int);

        btkSetMacro(InputModelImage, ModelImage::Pointer);
        btkGetMacro(InputModelImage, ModelImage::Pointer);

        /**
         * @brief Update the process.
         */
        virtual void Update();

        /**
         * @brief Get modeling at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Model response in direction direction at cindex in image space.
         */
        virtual float ModelAt(ContinuousIndex cindex, btk::GradientDirection direction);

        /**
         * @brief Get modeling at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Model response in direction direction at cindex in image space.
         */
        virtual std::vector< float > ModelAt(ContinuousIndex cindex, std::vector< btk::GradientDirection > &directions);

        /**
         * @brief Get modeling at continuous index.
         * @param cindex Location in the image space.
         * @return Model response at cindex in image space.
         */
        virtual std::vector< float > ModelAt(ContinuousIndex cindex);

        /**
         * @brief Get modeling at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Model response in direction direction at point in physical space.
         */
        virtual float ModelAt(PhysicalPoint point, btk::GradientDirection direction);

        /**
         * @brief Get modeling at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Model response in direction direction at point in physical space.
         */
        virtual std::vector< float > ModelAt(PhysicalPoint point, std::vector< btk::GradientDirection > &directions);

        /**
         * @brief Get modeling at physical point.
         * @param point Point in the physical space.
         * @return Model response at point in physical space.
         */
        virtual std::vector< float > ModelAt(PhysicalPoint point);

        /**
         * @brief Get signal at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Signal response in direction direction at cindex in image space.
         */
        virtual float SignalAt(ContinuousIndex cindex, btk::GradientDirection direction);

        /**
         * @brief Get signal at continuous index and gradient direction.
         * @param cindex Location in the image space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Signal response in direction direction at cindex in image space.
         */
        virtual std::vector< float > SignalAt(ContinuousIndex cindex, std::vector< btk::GradientDirection > &directions);

        /**
         * @brief Get signal at continuous index.
         * @param cindex Location in the image space.
         * @return Signal response at cindex in image space.
         */
        virtual std::vector< float > SignalAt(ContinuousIndex cindex);

        /**
         * @brief Get signal at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param direction Gradient direction were the model response is wanted.
         * @return Signal response in direction direction at point in physical space.
         */
        virtual float SignalAt(PhysicalPoint point, btk::GradientDirection direction);

        /**
         * @brief Get signal at physical point and gradient direction.
         * @param point Point in the physical space.
         * @param directions Gradient directions were the model response is wanted.
         * @return Signal response in direction direction at point in physical space.
         */
        virtual std::vector< float > SignalAt(PhysicalPoint point, std::vector< btk::GradientDirection > &directions);

        /**
         * @brief Get signal at physical point.
         * @param point Point in the physical space.
         * @return Signal response at point in physical space.
         */
        virtual std::vector< float > SignalAt(PhysicalPoint point);

        /**
         * @brief Get mean directions at a location in the physical space.
         * @param point Point in the physical space.
         * @return Vector of mean directions of local model at a physical location point.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(ContinuousIndex cindex);

        /**
         * @brief Get mean directions in a solid angle formed by previous vector and search angle at a location in the physical space.
         * @param point Point in the physical space.
         * @param vector Previous vector
         * @param angle Angle of search.
         * @return Vector of mean directions of local model in a solid angle formed by previous vector and search angle at a physical location point.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(PhysicalPoint point, btk::GradientDirection vector, float angle);

        /**
         * @brief Get mean directions at a location in the image space.
         * @param cindex Continuous index in the image space.
         * @return Vector of mean directions of local model at a location in image space.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(PhysicalPoint point);

        /**
         * @brief Get mean directions in a solid angle formed by previous vector and search angle at a location in the physical space.
         * @param cindex Continuous index in the image space.
         * @param vector Previous vector
         * @param angle Angle of search.
         * @return Vector of mean directions of local model in a solid angle formed by previous vector and search angle at a physical location point.
         */
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(ContinuousIndex cindex, btk::GradientDirection vector, float angle);

    protected:
        /**
         * @brief Constructor.
         */
        TensorModel();

        /**
         * @brief Destructor.
         */
        virtual ~TensorModel();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /**
         * @brief Compute model from a diffusion tensor in a particular direction.
         * @param tensor Diffusion tensor.
         * @param direction Gradient direction.
         * @return Model response in direction direction computed from diffusion tensor tensor.
         */
        virtual float ModelAt(ModelImage::PixelType tensor, btk::GradientDirection direction);

        /**
         * @brief Compute model from a diffusion tensor in a particular direction.
         * @param tensor Diffusion tensor.
         * @param directions Gradient directions.
         * @return Model response in direction direction computed from diffusion tensor tensor.
         */
        virtual std::vector< float > ModelAt(ModelImage::PixelType tensor, std::vector< btk::GradientDirection > &directions);

        /**
         * @brief Compute signal from a diffusion tensor in a particular direction.
         * @param tensor Diffusion tensor.
         * @param direction Gradient direction.
         * @return Signal response in direction direction computed from diffusion tensor tensor.
         */
        virtual float SignalAt(ModelImage::PixelType tensor, btk::GradientDirection direction);

        /**
         * @brief Compute signal from a diffusion tensor in a particular direction.
         * @param tensor Diffusion tensor.
         * @param directions Gradient directions.
         * @return Signal response in direction direction computed from diffusion tensor tensor.
         */
        virtual std::vector< float > SignalAt(ModelImage::PixelType tensor, std::vector< btk::GradientDirection > &directions);

        /**
         * @brief Convert a point in physical space of input model image to a continuous index.
         * @param point Point in physical space.
         * @return Point in image space (continuous index).
         */
        virtual ContinuousIndex TransformPhysicalPointToContinuousIndex(PhysicalPoint point);

    private:
        /**
         * @brief B-value used during the data acquisition.
         */
        unsigned int m_BValue;

        /**
         * @brief Tensor image estimated by reconstruction filter.
         */
        ModelImage::Pointer m_InputModelImage;

        /**
         * @brief Interpolation function (linear).
         */
        InterpolateModelFunction::Pointer m_ModelImageFunction;
};

} // namespace btk

#endif // BTK_TENSOR_MODEL_H
