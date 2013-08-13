/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 07/01/2012
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

#ifndef BTK_FEATURE_SELECTION_ALGORITHM_H
#define BTK_FEATURE_SELECTION_ALGORITHM_H

// STL includes
#include "vector"

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "vnl/vnl_matrix.h"

// Local includes
#include "btkMacro.h"
#include "btkFeatureSelectionCostFunction.h"


namespace btk
{

/**
 * @brief Algorithm base class for feature selection of parameters.
 * @author Julien Pontabry
 * @ingroup FeatureSelection
 */
class FeatureSelectionAlgorithm : public itk::ProcessObject
{
    public:
        typedef FeatureSelectionAlgorithm       Self;
        typedef itk::ProcessObject              Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        itkTypeMacro(FeatureSelectionAlgorithm,itk::ProcessObject);
        itkNewMacro(Self);

        btkSetMacro(InputParameters, vnl_matrix< double > *);
        btkGetMacro(InputParameters, vnl_matrix< double > *);

        btkGetMacro(WeightsVector, vnl_vector< short > *);

        btkGetMacro(EnergyVector, vnl_vector< double > *);

        btkSetMacro(CostFunction, FeatureSelectionCostFunction::Pointer);
        btkGetMacro(CostFunction, FeatureSelectionCostFunction::Pointer);

        btkGetMacro(MaxNumberOfParameters, unsigned int);
        btkSetMacro(MaxNumberOfParameters, unsigned int);

        btkGetMacro(UnexplainedVariance, double);

    protected:
        /**
         * @brief Constructor.
         */
        FeatureSelectionAlgorithm();

        /**
         * @brief Destructor.
         */
        ~FeatureSelectionAlgorithm();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /**
         * @brief Initialize the filter (called by update method).
         */
        virtual void Initialize();

    protected:
        /**
         * @brief Input parameters of the feature selection algorithm.
         */
        vnl_matrix< double > *m_InputParameters;

        /**
         * @brief Weights vector corresponding to the reduced parameter space.
         */
        vnl_vector< short > *m_WeightsVector;

        /**
         * @brief Energy vector corresponding to the reduced parameter space.
         */
        vnl_vector< double > *m_EnergyVector;

        /**
         * @brief Cost function used in optimization.
         */
        FeatureSelectionCostFunction::Pointer m_CostFunction;

        /**
         * @brief Maximal number of parameters.
         */
        unsigned int m_MaxNumberOfParameters;

        /**
         * @brief Unexplained variance.
         */
        double m_UnexplainedVariance;
};

} // namespace btk

#endif // BTK_FEATURE_SELECTION_ALGORITHM_H
