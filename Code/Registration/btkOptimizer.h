/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 13/12/2012
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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
#ifndef _btkOptimizer_h
#define _btkOptimizer_h

#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkSmartPointer.h"
#include "itkSingleValuedCostFunction.h"

namespace btk
{
class Optimizer : public itk::SingleValuedNonLinearOptimizer
{
public:
        typedef Optimizer   Self;
        typedef itk::SingleValuedNonLinearOptimizer Superclass;
        typedef itk::SmartPointer< Self > Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(Optimizer, itk::SingleValuedNonLinearOptimizer);

        /**  Parameters type.
         *  It defines a position in the optimization search space. */
        typedef Superclass::ParametersType ParametersType;

        /** InternalParameters typedef. */
        typedef   vnl_vector< double > InternalParametersType;

        /** Type of the Cost Function   */
        typedef  itk::SingleValuedCostFunction  CostFunctionType;
        typedef  CostFunctionType::Pointer CostFunctionPointer;

        /** Return Current Value */
        itkGetConstReferenceMacro(CurrentCost, MeasureType);
        MeasureType GetValue() const { return this->GetCurrentCost(); }

        virtual void StartOptimization(void);

        /** When users call StartOptimization, this value will be set false.
         * By calling StopOptimization, this flag will be set true, and
         * optimization will stop at the next iteration. */
        void StopOptimization()
        { m_Stop = true; }

        /** Plug in a Cost Function into the optimizer  */
        virtual void SetCostFunction(CostFunctionType *costFunction);

        const std::string GetStopConditionDescription() const;

        itkSetMacro( MaximumNumberOfIterations, unsigned int );
        itkGetConstMacro( MaximumNumberOfIterations, unsigned int );

    protected:
        Optimizer();
        virtual ~Optimizer();
        void PrintSelf(std::ostream &os, itk::Indent indent) const;

    protected:

        MeasureType m_CurrentValue;
        MeasureType m_CurrentCost;

        bool m_Stop;

        unsigned int m_MaximumNumberOfIterations;
        unsigned int m_CurrentIteration;

        std::ostringstream m_StopConditionDescription;

};
}// namespace

#endif
