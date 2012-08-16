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

#ifndef BTK_ORIENTATION_DIFFUSION_FUNCTION_MODEL_H
#define BTK_ORIENTATION_DIFFUSION_FUNCTION_MODEL_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkLinearInterpolateImageFunction.h"

// Local includes
#include "btkDiffusionModel.h"
#include "btkSphericalHarmonicsDiffusionDecompositionFilter.h"

namespace btk
{

/**
 * @brief Orientation diffusion function computed by using spherical harmonics coefficients.
 * @author Julien Pontabry
 * @ingroup Diffusion
 */
class OrientationDiffusionFunctionModel : public btk::DiffusionModel
{
    public:
        typedef OrientationDiffusionFunctionModel Self;
        typedef btk::DiffusionModel               Superclass;
        typedef itk::SmartPointer< Self >         Pointer;
        typedef itk::SmartPointer< const Self >   ConstPointer;

        typedef btk::SphericalHarmonicsDiffusionDecompositionFilter::OutputImageType ModelImage;
        typedef itk::LinearInterpolateImageFunction< ModelImage,float >              InterpolateModelFunction;
        typedef Superclass::PhysicalPoint                                            PhysicalPoint;
        typedef Superclass::ContinuousIndex                                          ContinuousIndex;

        itkNewMacro(Self);

        itkTypeMacro(OrientationDiffusionFunctionModel,btk::DiffusionModel);

        btkSetMacro(InputModelImage, ModelImage::Pointer);
        btkGetMacro(InputModelImage, ModelImage::Pointer);

        virtual void Update();

        virtual float ModelAt(ContinuousIndex cindex, btk::GradientDirection direction);
        virtual std::vector< float > ModelAt(ContinuousIndex cindex);
        virtual float ModelAt(PhysicalPoint point, btk::GradientDirection direction);
        virtual std::vector< float > ModelAt(PhysicalPoint point);
        virtual float SignalAt(ContinuousIndex cindex, btk::GradientDirection direction);
        virtual std::vector< float > SignalAt(ContinuousIndex cindex);
        virtual float SignalAt(PhysicalPoint point, btk::GradientDirection direction);
        virtual std::vector< float > SignalAt(PhysicalPoint point);
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(ContinuousIndex cindex);
        virtual std::vector< btk::GradientDirection > MeanDirectionsAt(PhysicalPoint point);

    protected:
        /**
         * @brief Constructor.
         */
        OrientationDiffusionFunctionModel();

        /**
         * @brief Destructor.
         */
        virtual ~OrientationDiffusionFunctionModel();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        virtual float ModelAt(ModelImage::PixelType tensor, btk::GradientDirection direction);

        virtual float SignalAt(ModelImage::PixelType tensor, btk::GradientDirection direction);

        virtual ContinuousIndex TransformPhysicalPointToContinuousIndex(PhysicalPoint point);

    private:
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

#endif // BTK_ORIENTATION_DIFFUSION_FUNCTION_MODEL_H
