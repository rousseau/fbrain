/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 11/03/2013
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

#ifndef BTK_DIFFUSION_SIGNAL_H
#define BTK_DIFFUSION_SIGNAL_H

// STL includes
#include "vector"

// ITK includes
#include "itkMacro.h"
#include "itkSmartPointer.h"
#include "itkVectorImage.h"
#include "itkPoint.h"
#include "itkContinuousIndex.h"
#include "itkLinearInterpolateImageFunction.h"

// Local includes
#include "btkMacro.h"
#include "btkGradientDirection.h"

namespace btk
{

class DiffusionSignal : public itk::VectorImage< double,3 >
{
    public:
        typedef DiffusionSignal                 Self;
        typedef itk::VectorImage< double,3 >    Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef itk::LinearInterpolateImageFunction< Self > InterpolatedDiffusionSignal;

        typedef InterpolatedDiffusionSignal::ContinuousIndexType ContinuousIndex;
        typedef Superclass::PointType                            PhysicalPoint;


        itkNewMacro(Self);

        itkTypeMacro(DiffusionSignal,itk::VectorImage);


        btkSetMacro(GradientTable, std::vector< btk::GradientDirection >);
        btkGetMacro(GradientTable, std::vector< btk::GradientDirection >);

        btkSetMacro(PseudoResidualsStdDeviation, std::vector< double >);
        btkGetMacro(PseudoResidualsStdDeviation, std::vector< double >);


        /**
         * @brief Get signal at continuous index.
         * @param cindex Location in the image space.
         * @return Signal response at cindex in image space.
         */
        virtual PixelType SignalAt(ContinuousIndex cindex);

        /**
         * @brief Get signal at physical point.
         * @param point Point in the physical space.
         * @return Signal response at point in physical space.
         */
        virtual PixelType SignalAt(PhysicalPoint point);

    protected:
        /**
         * @brief Constructor.
         */
        DiffusionSignal();

        /**
         * @brief Desctructor.
         */
        virtual ~DiffusionSignal();

    private:
        /**
         * @brief Gradient directions.
         */
        std::vector< btk::GradientDirection > m_GradientTable;

        /**
         * @brief Standard deviation of the pseudo residuals of each gradient image.
         */
        std::vector< double > m_PseudoResidualsStdDeviation;

        /**
         * @brief Linear interpolation of the diffusion signal.
         */
        InterpolatedDiffusionSignal::Pointer m_Interpolation;
};

} // namespace btk

#endif // BTK_DIFFUSION_SIGNAL_H
