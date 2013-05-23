/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 
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

#ifndef BTKSMARTSTEPGRADIENTDESCENTOPTIMIZER_H
#define BTKSMARTSTEPGRADIENTDESCENTOPTIMIZER_H

#include "btkOptimizer.h"
#include "btkMacro.h"

namespace btk
{
/**
 * This class is a gradient descent optimizer, the step at each iteration is computed with a smart search
 * allong gradient direction.
 * This optimizer can be used with an itk::SingleValueCostFunction
 *
 * @author Marc Schweitzer
 * \ingroup Registration
 */


class SmartStepGradientDescentOptimizer : public btk::Optimizer
{
    public:
        typedef SmartStepGradientDescentOptimizer   Self;
        typedef Optimizer Superclass;
        typedef itk::SmartPointer< Self > Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(SmartStepGradientDescentOptimizer, Optimizer);

        /**  Parameters type.
             *  It defines a position in the optimization search space. */
        typedef Superclass::ParametersType ParametersType;

        /** InternalParameters typedef. */
        typedef   vnl_vector< double > InternalParametersType;

        /** Type of the Cost Function   */
        typedef  itk::SingleValuedCostFunction  CostFunctionType;
        typedef  CostFunctionType::Pointer CostFunctionPointer;

        /** Type of Measure */

        typedef CostFunctionType::MeasureType MeasureType;

        /** Compute Optimization */
        virtual void StartOptimization(void);

        /** Activate of Desactivate optimization of few parameters (by a unit vector),
         * by default 1 for all parameters */
        void SetOptimizedParameters(ParametersType _Op)
        {
            m_OptimizedParameters = _Op;

            if(m_OptimizedParameters.min_value() == 0)
            {
                m_OptimizeAllParameters = false;
            }
        }

        /** Get current cost */
        MeasureType GetValue() const
        {
            return this->m_CurrentValue;
        }

        /** Set/Get NumberOfIterations*/
        btkSetMacro(NumberOfIterations, unsigned int);
        btkGetMacro(NumberOfIterations, unsigned int);

        /** Set/Get MaxStep length */
        btkSetMacro(MaxStep,double);
        btkGetMacro(MaxStep,double);

        /** Set/Get MinStep length */
        btkSetMacro(MinStep, double);
        btkGetMacro(MinStep, double);

        /** Set/Get Epsilon value between two consecutive cost */
        btkSetMacro(Epsilon, double);
        btkGetMacro(Epsilon, double);

        /** Set/Get Min Bounds */
        btkSetMacro(MinBounds, ParametersType);
        btkGetMacro(MinBounds, ParametersType);

        /** Set/Get max Bounds */
        btkSetMacro(MaxBounds, ParametersType);
        btkGetMacro(MaxBounds, ParametersType);

        /** Set/Get verbose mode */
        btkSetMacro(VerboseMode, bool);
        btkGetMacro(VerboseMode, bool);

        /** Set/Get Use bounds mode */
        btkSetMacro(UseBounds, bool);
        btkGetMacro(UseBounds, bool);

        /** Set/Get Samples for step searching */
        btkSetMacro(Samples, double);
        btkGetMacro(Samples, double);




    protected:

        SmartStepGradientDescentOptimizer();
        virtual ~SmartStepGradientDescentOptimizer(){}

        void PrintSelf(std::ostream &os, itk::Indent indent) const
        {
            Superclass::PrintSelf(os, indent);
        }

    private :

        bool IsInsideBounds(ParametersType _xIn);

        ParametersType ReverseBounds(ParametersType _xIn);

        double SearchStep(ParametersType _x , DerivativeType _gx);

        inline DerivativeType CheckParameters(DerivativeType _gx);

        ParametersType UpdateParameters(ParametersType _x);

        double SmartSearch(ParametersType _x , DerivativeType _gx);

        unsigned int m_NumberOfIterations;

        double m_Epsilon;

        ParametersType m_X;

        DerivativeType m_gX;

        double m_MaxStep;
        double m_MinStep;

        MeasureType m_CurrentValue;

        ParametersType m_OptimizedParameters;

        bool m_OptimizeAllParameters;

        ParametersType m_MinBounds;
        ParametersType m_MaxBounds;

        bool m_VerboseMode;
        bool m_UseBounds;

        double m_Samples;

};

}

#endif // BTKSMARTSTEPGRADIENTDESCENTOPTIMIZER_H
