/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 16/05/2013
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
#include "vnl/vnl_math.h"
#include "itkSingleValuedCostFunction.h"

#include "btkSmartStepGradientDescentOptimizer.h"



/**
 *  The objectif function is the quadratic form:
 *
 *  1/2 x^T A x - b^T x
 *
 *  Where A is a matrix and b is a vector
 *  The system in this example is:
 *
 *     | 3  2 ||x|   | 2|   |0|
 *     | 2  6 ||y| + |-8| = |0|
 *
 *
 *   the solution is the vector | 2 -2 |
 *
 * \class gradientCostFunction
 */
class gradientCostFunction : public itk::SingleValuedCostFunction
{
public:

  typedef gradientCostFunction            Self;
  typedef itk::SingleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>         Pointer;
  typedef itk::SmartPointer<const Self>   ConstPointer;
  itkNewMacro( Self );
  itkTypeMacro( gradientCostFunction, SingleValuedCostFunction );

  enum { SpaceDimension=2 };

  typedef Superclass::ParametersType      ParametersType;
  typedef Superclass::DerivativeType      DerivativeType;
  typedef Superclass::MeasureType         MeasureType;

  gradientCostFunction()
  {
  }


  MeasureType  GetValue( const ParametersType & parameters ) const
  {

    double x = parameters[0];
    double y = parameters[1];

    std::cout << "GetValue( ";
    std::cout << x << " ";
    std::cout << y << ") = ";

    MeasureType measure = 0.5*(3*x*x+4*x*y+6*y*y) - 2*x + 8*y;

    std::cout << measure << std::endl;

    return measure;

  }

  void GetDerivative( const ParametersType & parameters,
                            DerivativeType & derivative ) const
  {

    double x = parameters[0];
    double y = parameters[1];

    std::cout << "GetDerivative( ";
    std::cout << x << " ";
    std::cout << y << ") = ";

    DerivativeType temp(SpaceDimension);
    temp.Fill( 0 );
    derivative = temp;
    derivative[0] = 3 * x + 2 * y -2;
    derivative[1] = 2 * x + 6 * y +8;

    std::cout << derivative << std::endl;

  }


  unsigned int GetNumberOfParameters(void) const
    {
    return SpaceDimension;
    }

private:


};

int main(int, char* [] )
{

  //---------------

    std::cout << "Gradient Descent Optimizer Test ";
    std::cout << std::endl << std::endl;

    typedef  btk::SmartStepGradientDescentOptimizer  OptimizerType;

    // Declaration of a itkOptimizer
    OptimizerType::Pointer  optimizer = OptimizerType::New();


    // Declaration of the CostFunction
    gradientCostFunction::Pointer costFunction = gradientCostFunction::New();


    optimizer->SetCostFunction( costFunction.GetPointer() );


    typedef gradientCostFunction::ParametersType    ParametersType;

    const unsigned int spaceDimension =
                        costFunction->GetNumberOfParameters();

    // We start not so far from  | 2 -2 |
    ParametersType  initialPosition( spaceDimension );

    initialPosition[0] =  100;
    initialPosition[1] = -100;

    optimizer->SetNumberOfIterations(1000);
    optimizer->SetMaxStep(2.0);
    optimizer->SetMinStep(0.01);
    optimizer->SetUseBounds(false);

    optimizer->SetInitialPosition( initialPosition );
    optimizer->SetVerboseMode(true);

    try
      {
      optimizer->StartOptimization();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << "Exception thrown ! " << std::endl;
      std::cout << "An error occurred during Optimization" << std::endl;
      std::cout << "Location    = " << e.GetLocation()    << std::endl;
      std::cout << "Description = " << e.GetDescription() << std::endl;
      return EXIT_FAILURE;
      }

    ParametersType finalPosition = optimizer->GetCurrentPosition();
    double error =  optimizer->GetValue();
    std::cout << "Solution        = (";
    std::cout << finalPosition[0] << ",";
    std::cout << finalPosition[1] << ")" << std::endl;


    std::cout<<"Error       = "<<error<<std::endl;
    //
    // check results to see if it is within range
    //
    bool pass = true;
    double trueParameters[2] = { 2, -2 };
    for( unsigned int j = 0; j < 2; j++ )
      {
      if( vnl_math_abs( finalPosition[j] - trueParameters[j] ) > 0.01 )
        pass = false;
      }

    if( !pass )
      {
      std::cout << "Test failed." << std::endl;
      return EXIT_FAILURE;
      }

    std::cout << "Test passed." << std::endl;
    return EXIT_SUCCESS;


}
