/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 18/07/2012
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

#ifndef BTK_SPHERICAL_HARMONICS_DIFFUSION_DECOMPOSITION_FILTER_H
#define BTK_SPHERICAL_HARMONICS_DIFFUSION_DECOMPOSITION_FILTER_H

// ITK includes
#include "itkImageToImageFilter.h"
#include "itkVariableSizeMatrix.h"
#include "itkVectorImage.h"

// Local includes
#include "btkMacro.h"
#include "btkDiffusionSequence.h"

namespace btk
{

/**
 * @brief Compute the spherical harmonics decomposition model to diffusion sequence.
 * @author Julien Pontabry
 * @ingroup Diffusion
 */
class SphericalHarmonicsDiffusionDecompositionFilter : public itk::ImageToImageFilter< btk::DiffusionSequence,itk::VectorImage< float,3 > >
{
    public:
        typedef SphericalHarmonicsDiffusionDecompositionFilter                                Self;
        typedef itk::ImageToImageFilter< btk::DiffusionSequence,itk::VectorImage< float,3 > > Superclass;
        typedef itk::SmartPointer< Self >                                                     Pointer;
        typedef itk::SmartPointer< const Self >                                               ConstPointer;

        typedef itk::VariableSizeMatrix< float > Matrix;

        itkNewMacro(Self);

        itkTypeMacro(SphericalHarmonicsDiffusionDecompositionFilter,itk::ImageToImageFilter);

        /**
         * @brief Set input diffusion sequence.
         * @param input Diffusion sequence.
         */
        void SetInput(const btk::DiffusionSequence::Pointer input);

        btkSetMacro(SphericalHarmonicsOrder, unsigned short);
        btkGetMacro(SphericalHarmonicsOrder, unsigned short);

        btkSetMacro(RegularizationParameter, float);
        btkGetMacro(RegularizationParameter, float);

    protected:
        /**
         * @brief Constructor.
         */
        SphericalHarmonicsDiffusionDecompositionFilter();

        /**
         * @brief Destructor.
         */
        virtual ~SphericalHarmonicsDiffusionDecompositionFilter();

        /**
         * @brief Run the filter and generate output data.
         */
        virtual void GenerateData();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:
        /**
         * @brief Compute the matrix of SH basis corresponding to the gradient table.
         */
        void ComputeSphericalHarmonicsBasisMatrix();

        /**
         * @brief Compute the regularization matrix for linear regression.
         */
        void ComputeRegularizationMatrix();

        /**
         * @brief Compute the linear regression matrix.
         */
        void ComputeTransitionMatrix();

    private:
        /** Diffusion sequence. */
        btk::DiffusionSequence::Pointer m_InputDiffusionSequence;

        /** Order of the spherical harmonics decomposition. */
        unsigned short m_SphericalHarmonicsOrder;

        /** Regularization parameter of the linear regression. */
        float m_RegularizationParameter;

        /** Number of spherical harmonics coefficients. */
        unsigned short m_NumberOfSHCoefficients;

        /** Spherical harmonics basis matrix. */
        Self::Matrix m_SphericalHarmonicsBasisMatrix;

        /** Regularization matrix (Laplace-Beltrami matrix). */
        Self::Matrix m_RegularizationMatrix;

        /** Transition matrix for linear regression. */
        Self::Matrix m_TransitionMatrix;
};

} // namespace btk

#endif // BTK_SPHERICAL_HARMONICS_DIFFUSION_DECOMPOSITION_FILTER_H
