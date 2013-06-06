    /*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/06/2013
  Author(s): Fr�deric Champ (champ(at)unistra.fr)
  
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

#ifndef BTKWEIGHTEDESTIMATIONBASE_H
#define BTKWEIGHTEDESTIMATIONBASE_H

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
 * @brief Compute the spherical harmonics decomposition model to diffusion sequence and signal estimation.
 * Use modified par of btkSphericalHarmonicsDiffusionDecompositionFilter and btkOrientationDiffusionFunctionModel
 * @author Frederic CHAMP
 * @ingroup Diffusion
 */
class WeightedEstimationBase : public itk::ProcessObject
{
    public:
        typedef WeightedEstimationBase                  Self;
        typedef ProcessObject                           Superclass;
        typedef itk::SmartPointer< Self >               Pointer;
        typedef itk::SmartPointer< const Self >         ConstPointer;

        itkNewMacro(Self);

        /** Diffusion Sequence type definitions */
        typedef btk::DiffusionSequence                  SequenceType;
        typedef SequenceType::IndexType                 SequenceIndexType;

        /** Containers type definition */
        typedef btk::DiffusionSequence::GradientTable   GradientTableType;
        typedef std::vector<SequenceIndexType>          IndexVectorType;
        typedef std::vector<float>                      VectorType;
        typedef itk::VariableSizeMatrix< float >        Matrix;


        /**
         * @brief Set/Get Order the Spherical Harmonics decomposition (default = 4 ).
         * @param Spherical Harmonics Order
         */
        btkSetMacro(SphericalHarmonicsOrder, unsigned short);
        btkGetMacro(SphericalHarmonicsOrder, unsigned short);

        /**
         * @brief Set/Get Regularization parameter of the linear regression ( default = 0.006 ).
         * @param Regularization Parameter
         */
        btkSetMacro(RegularizationParameter, float);
        btkGetMacro(RegularizationParameter, float);

        /**
         * @brief Set/Get the the value of the point to be estimated
         *        If the value is not set, weights will be only function of distance
         * @param Value of the query point.
         */
        btkSetMacro(QueryValue,unsigned short);
        btkGetMacro(QueryValue,unsigned short);

        /**
         * @brief Set/Get Neighbors distances from the query point (mm).
         * @param Vector of Neighbors distances
         */
        btkSetMacro(NeighborsDistances,VectorType);
        btkGetMacro(NeighborsDistances,VectorType);

        /**
         * @brief Set/Get Signal Values of neighbors
         * @param Vector of Signal Values (Must have the same size as "Neighbors" vector)
         */
        btkSetMacro(SignalValues,VectorType);
        btkGetMacro(SignalValues,VectorType);

        /**
         * @brief Set/Get Gradient direction of neighbors
         * @param Vector of gradient directions (Must have the same size as "Neighbors" vector)
         */
        btkSetMacro(GradientTable,GradientTableType);
        btkGetMacro(GradientTable,GradientTableType);

        /**
         * @brief Set/Get Spherical Resolution of the modeling reconstruction (number of point on the unit sphere, default 100).
         * @param Spherical Resolution
         */
        btkSetMacro(SphericalResolution, unsigned int);
        btkGetMacro(SphericalResolution, unsigned int);

        /**
         * @brief Set/Get the sigma parameters of the kernel used to calculate weights (default = 1.0)
         * @param Sigma
         */
        btkSetMacro(Sigma, float);
        btkGetMacro(Sigma, float);

        /**
         * @brief The intialize method
         */
        void Initialize();


        /**
         * @brief The execute method .
         */
        void Update();

        /**
         * @brief Compute signal from a diffusion tensor in a particular direction.
         * @param direction Gradient direction.
         * @return Signal response in direction direction computed from diffusion tensor tensor.
         */
        virtual float SignalAt( btk::GradientDirection direction);


    protected:
        /**
         * @brief Constructor.
         */
        WeightedEstimationBase();

        /**
         * @brief Destructor.
         */
        virtual ~WeightedEstimationBase();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:

        /**
         * @brief Compute the weighted Matrix
         */
        void ComputeWeightedMatrix();

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

        /** sigma of the kernel used to calculate wieghts */
        float m_Sigma;

        /** Minimum value */
        float m_VMin;

        /** Maximum value */
        float m_VMax;

        /** Minimum distance */
        float m_DMin;

        /** Maximum distance */
        float m_DMax; 

        /** Regularization parameter of the linear regression. */
        float m_RegularizationParameter;

        /** Query pont value */
        unsigned short      m_QueryValue;

        /** Order of the spherical harmonics decomposition. */
        unsigned short m_SphericalHarmonicsOrder;

        /** Number of spherical harmonics coefficients. */
        unsigned short m_NumberOfSHCoefficients;

        /** Number of Neighbors */
        unsigned int        m_numberOfNeighbors;

        /** Spherical resolution of the modeling reconstruction (number of point on the unit sphere) */
        unsigned int m_SphericalResolution;

        /** Vector of neighbors distances form the query point */
        VectorType    m_NeighborsDistances;

        /** Vector of neigignal values */
        VectorType    m_SignalValues;

        /** Vecotr of neighbors gradient direction */
        GradientTableType   m_GradientTable;

        /** Vecotr shperical harmonic coefficients */
        VectorType m_Coefficients;

         /** Sampled directions on the unit sphere used by modeling reconstruction. */
        GradientTableType m_Directions;

        /** Spherical harmonics basis matrix. */
        Self::Matrix m_CovarianceMatrix;

        /** Spherical harmonics basis matrix. */
        Self::Matrix m_SphericalHarmonicsBasisMatrix;

        /** Regularization matrix (Laplace-Beltrami matrix). */
        Self::Matrix m_RegularizationMatrix;

        /** Transition matrix for linear regression. */
        Self::Matrix m_TransitionMatrix;

        /** Matrix of weights */
        Self::Matrix m_WeightsMatrix;

        Self::Matrix m_SignalMatrix;


};

} // namespace btk

#endif // BTK_WEIGHTED_SPHERICAL_HARMONICS_DIFFUSION_DECOMPOSITION_FILTER_H
