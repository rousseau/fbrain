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

#ifndef BTK_SPACE_REDUCTION_COST_FUNCTION_H
#define BTK_SPACE_REDUCTION_COST_FUNCTION_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

// Local includes
#include "btkMacro.h"


namespace btk
{

/**
 * @brief Base class for cost function used by feature selection algorithm.
 * @author Julien Pontabry
 * @date 09/01/2013
 * @ingroup FeatureSelection
 */
class FeatureSelectionCostFunction : public itk::ProcessObject
{
    public:
        typedef FeatureSelectionCostFunction    Self;
        typedef itk::ProcessObject              Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        itkTypeMacro(FeatureSelectionCostFunction,itk::ProcessObject);


        btkSetMacro(InputParameters, vnl_matrix< double > *);
        btkGetMacro(InputParameters, vnl_matrix< double > *);

        btkSetMacro(ImagesWeightVector, vnl_vector< double > *);
        btkGetMacro(ImagesWeightVector, vnl_vector< double > *);

        btkSetMacro(CurrentParameters, vnl_matrix< double > *);
        btkGetMacro(CurrentParameters, vnl_matrix< double > *);

        /**
         * @brief Evaluate the cost function with activated vectors.
         * @return Value of the cost function with activated vectors.
         */
        virtual double Evaluate() = 0;

        /**
         * @brief Evaluate cost function with parameter activation.
         * @param activatedIndex Index of the parameter to test activation.
         * @return The value of the cost function with the activation of the parameter with index activatedIndex.
         */
        virtual double EvaluateActivation(unsigned int activatedIndex) = 0;

        /**
         * @brief Evaluate cost function with parameter desactivation.
         * @param desactivatedIndex Index of the parameter to test desactivation.
         * @return The value of the cost function with the desactivation of the parameter with index activatedIndex.
         */
        virtual double EvaluateDesactivation(unsigned int desactivatedIndex) = 0;

        /**
         * @brief Activate parameter with index.
         * @param index Index of the parameter to activate.
         */
        virtual void ActivateParameters(unsigned int index) = 0;

        /**
         * @brief Desactivate parameter with index.
         * @param index Index of the parameter to desactivate.
         */
        virtual void DesactivateParameters(unsigned int index) = 0;

        /**
         * @brief Get the residuals of the current parameters.
         * @return A matrix containing the residuals of each feature (rows) and each image (columns).
         */
        virtual vnl_matrix< double > GetResiduals() = 0;

        /**
         * @brief Get the reconstruction given the current parameters.
         * This method actually allocate space memory for the matrix and return it as a reference.
         * @param parameters Matrix of parameters used to reconstruct the data.
         * @return A reconstructed matrix from the given matrix of parameters.
         */
        virtual vnl_matrix< double > &GetReconstruction(vnl_matrix< double > &parameters) = 0;

        /**
         * @brief Initialize the cost function.
         */
        virtual void Initialize();

    protected:
        /**
         * @brief Constructor.
         */
        FeatureSelectionCostFunction();

        /**
         * @brief Destructor.
         */
        ~FeatureSelectionCostFunction();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    protected:
        /**
         * @brief Input parameters of the feature selection algorithm.
         */
        vnl_matrix< double > *m_InputParameters;

        /**
         * @brief Weights for each images
         */
        vnl_vector< double > *m_ImagesWeightVector;

        /**
         * @brief Number of parameters used.
         */
        unsigned int m_NumberOfParameters;

        /**
         * @brief Current parameters activated.
         */
        vnl_matrix< double > *m_CurrentParameters;

        /**
         * @brief Define if the current parameters activated have been set manualy.
         */
        bool m_ExternalCurrentParameters;

        /**
         * @brief Number of activated parameters.
         */
        unsigned int m_numberOfActivatedParameters;

        /**
         * @brief Two times the number of parameters.
         */
        unsigned int m_TwoTimeNumberOfParameters;

        /**
         * @brief Number of samples vectors.
         */
        unsigned int m_NumberOfVectors;
};

} // namespace btk

#endif // BTK_SPACE_REDUCTION_COST_FUNCTION_H
