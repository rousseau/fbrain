/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 09/01/2013
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

#ifndef BTK_KERNEL_COST_FUNCTION_H
#define BTK_KERNEL_COST_FUNCTION_H

// VNL includes
#include "vnl/vnl_diag_matrix.h"


// Local includes
#include "btkFeatureSelectionCostFunction.h"

namespace btk
{

/**
 * @brief Cost function for feature selection using kernels.
 * @author Julien Pontabry
 * @ingroup FeatureSelection
 */
class NadarayaWatsonReconstructionErrorFunction : public FeatureSelectionCostFunction
{
    public:
        typedef NadarayaWatsonReconstructionErrorFunction Self;
        typedef FeatureSelectionCostFunction              Superclass;
        typedef itk::SmartPointer< Self >                 Pointer;
        typedef itk::SmartPointer< const Self >           ConstPointer;

        itkTypeMacro(NadarayaWatsonReconstructionErrorFunction,FeatureSelectionCostFunction);
        itkNewMacro(Self);

        btkGetMacro(BandwidthMatrix, vnl_diag_matrix< double > *);
        btkSetMacro(BandwidthMatrix, vnl_diag_matrix< double > *);


        /**
         * @brief Evaluate the cost function with activated vectors.
         * @return Value of the cost function with activated vectors.
         */
        virtual double Evaluate();

        /**
         * @brief Evaluate cost function with parameter activation.
         * @param activatedIndex Index of the parameter to test activation.
         * @return The value of the cost function with the activation of the parameter with index activatedIndex.
         */
        double EvaluateActivation(unsigned int activatedIndex);

        /**
         * @brief Evaluate cost function with parameter desactivation.
         * @param desactivatedIndex Index of the parameter to test desactivation.
         * @return The value of the cost function with the desactivation of the parameter with index activatedIndex.
         */
        double EvaluateDesactivation(unsigned int desactivatedIndex);

        /**
         * @brief Activate parameter with index.
         * @param index Index of the parameter to activate.
         */
        void ActivateParameters(unsigned int index);

        /**
         * @brief Desactivate parameter with index.
         * @param index Index of the parameter to desactivate.
         */
        void DesactivateParameters(unsigned int index);

        /**
         * @brief Initialize the cost function (compute some matrix and determinants).
         */
        virtual void Initialize();

    protected:
        /**
         * @brief Constructor.
         */
        NadarayaWatsonReconstructionErrorFunction();

        /**
         * @brief Compute the nadaraya watson multivariate kernel estimator.
         * @param data Reduced parameters matrix.
         * @param v Input vector for kernel estimate.
         * @return A vector estimated by kernel regression.
         */
        vnl_vector< double > NadarayaWatsonMultivariateKernelEstimator(vnl_matrix< double > &data, vnl_vector< double > &v);

        /**
         * @brief Compute the Nadaraya-Watson multivariate kernel estimator with parameter activation.
         * @param activatedIndex Index of parameter to test activation
         * @param v Input vector of the estimator
         * @return The vector estimated with parameter activation
         */
        vnl_vector< double > NadarayaWatsonMultivariateKernelEstimatorActivation(unsigned int activatedIndex, vnl_vector< double > &v);

        /**
         * @brief Compute the Nadaraya-Watson multivariate kernel estimator with parameter desactivation.
         * @param desactivatedIndex Index of parameter to test desactivation
         * @param v Input vector of the estimator
         * @return The vector estimated with parameter desactivation
         */
        vnl_vector< double > NadarayaWatsonMultivariateKernelEstimatorDesactivation(unsigned int desactivatedIndex, vnl_vector< double > &v);

        /**
         * @brief Compute the gaussian multivariate kernel.
         * @param v Multivariate data (vector).
         * @return Value of the kernel with the vector v.
         */
        double MultivariateGaussianKernel(vnl_vector< double > v);

        /**
         * @brief Compute the gaussian kernel for univariate data.
         * @param v Univariate value.
         * @return Value of the kernel with the univariate data v.
         */
        double GaussianKernel(double v);

    private:
        /**
         * @brief Bandwidth matrix.
         */
        vnl_diag_matrix< double >* m_BandwidthMatrix;

        /**
         * @brief Inverse of the bandwidth matrix.
         */
        vnl_diag_matrix< double > m_BandwidthMatrixInverse;

        /**
         * @brief Precomputed coefficient of the gaussian kernel.
         */
        double m_GaussianCoefficient;
};

} // namespace btk

#endif // BTK_KERNEL_COST_FUNCTION_H
