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

#ifndef BTKSUPERRESOLUTIONCOSTFUNCTIONITKWRAPPER_H
#define BTKSUPERRESOLUTIONCOSTFUNCTIONITKWRAPPER_H

#include "itkSingleValuedCostFunction.h"


#include "btkSuperResolutionCostFunction.hxx"
#include "btkMacro.h"

namespace btk
{
/**
 * @class SuperResolutionCostFunctionITKWrapper
 * @brief Wrapper ofr use the SuperResolutionCostFunction with an itk::Optimizer
 * @author Marc Schweitzer
 * @ingroup SuperResolution
 */
template < class TImage >
class SuperResolutionCostFunctionITKWrapper: public itk::SingleValuedCostFunction
{
    public:
        /** Standard class typedefs. */
        typedef SuperResolutionCostFunctionITKWrapper   Self;
        typedef itk::SingleValuedCostFunction               Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef TImage ImageType;
        typedef typename ImageType::PixelType VoxelType;

        typedef btk::SuperResolutionCostFunction< ImageType >   CostFunctionType;



        /** Run-time type information (and related methods). */
        itkTypeMacro(SuperResolutionCostFunctionITKWrapper, CostFunction);

        /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
        typedef double MeasureType;

        /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
        typedef Superclass::ParametersType      ParametersType;
        typedef Superclass::ParametersValueType ParametersValueType;

        /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
        typedef itk::Array< ParametersValueType > DerivativeType;

        /** This method returns the value of the cost function corresponding
      * to the specified parameters.    */
        virtual MeasureType GetValue(const ParametersType & parameters) const;

        /** This method returns the derivative of the cost function corresponding
      * to the specified parameters.   */
        virtual void GetDerivative(const ParametersType & parameters,
                                   DerivativeType & derivative) const;

        /** Return the number of parameters required to compute
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
        virtual unsigned int GetNumberOfParameters(void) const ;


        /** Return a pointer to a btk cost function */

        btkGetMacro(CostFunction,typename CostFunctionType::Pointer );

        btkSetMacro(NumberOfParameters, unsigned int);


    protected:
        SuperResolutionCostFunctionITKWrapper();
        virtual ~SuperResolutionCostFunctionITKWrapper();


    private:

        typename CostFunctionType::Pointer m_CostFunction;

        unsigned int m_NumberOfParameters;
};
}


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSuperResolutionCostFunctionITKWrapper.txx"
#endif

#endif // BTKSUPERRESOLUTIONCOSTFUNCTIONITKWRAPPER_H
