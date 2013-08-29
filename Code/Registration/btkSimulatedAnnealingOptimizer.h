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
#ifndef BTKRANDOMOPTIMIZER_H
#define BTKRANDOMOPTIMIZER_H

#include "btkOptimizer.h"
#include "btkMacro.h"

/* OTHERS */
#include "cfloat"
#include "cmath"
#include "algorithm"
#include "ctime"

namespace btk
{
/**
 * @class SimulatedAnnealingOptimizer
 * @brief This is an optimizer using Simulated Anneal Strategy.
 * It can be used with any class inherited form itk::SingleValuedCostFunction.
 *
 * @author Marc Schweitzer
 * @ingroup Registration
 */
class SimulatedAnnealingOptimizer : public btk::Optimizer
{

    public:
        typedef SimulatedAnnealingOptimizer   Self;
        typedef btk::Optimizer Superclass;
        typedef itk::SmartPointer< Self > Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(SimulatedAnnealingOptimizer, Optimizer);

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

        /** Get the current cost value */
        MeasureType GetValue() const
        {
            return this->m_CurrentValue;
        }

        /** Set/Get temperature */
        btkSetMacro(Temperature,double);
        btkGetMacro(Temperature, double);

        /** Set/Get Number of Iterations */
        btkSetMacro(Iteration,unsigned int);
        btkGetMacro(Iteration,unsigned int);


    protected:

        SimulatedAnnealingOptimizer();
        virtual ~SimulatedAnnealingOptimizer(){}
        void PrintSelf(std::ostream &os, itk::Indent indent) const
        {
            Superclass::PrintSelf(os, indent);
        }

    private :


        ParametersType GenerateRandomNeighboor(ParametersType _x);

        inline double Random()
        {
            return (double)rand()/(double)RAND_MAX;
        }

        ParametersType m_Range;
        ParametersType m_RangeNeighboor;

        double m_Temperature;
        unsigned int m_Iteration;
        MeasureType m_CurrentValue;


};

}

#endif // BTKRANDOMOPTIMIZER_H
