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

#include "btkSmartStepGradientDescentOptimizer.h"

namespace btk
{
//-------------------------------------------------------------------------------------------------
SmartStepGradientDescentOptimizer::SmartStepGradientDescentOptimizer(): m_MaxStep(2.0),
    m_MinStep(0.05),m_NumberOfIterations(100),m_Epsilon(0.00005),m_OptimizeAllParameters(true)
{
}
//-------------------------------------------------------------------------------------------------
void SmartStepGradientDescentOptimizer::StartOptimization()
{
    //TODO: Take scales into account !

    if(this->m_CostFunction.IsNull())
    {
        return;
    }

    unsigned int NumberOfParameters = this->m_CostFunction->GetNumberOfParameters();

    if(m_OptimizeAllParameters)
    {
        m_OptimizedParameters.SetSize(NumberOfParameters);
        m_OptimizedParameters.Fill(1);
    }


    m_StopConditionDescription.str("");
    m_StopConditionDescription<<this->GetNameOfClass()<<": ";

    this->InvokeEvent(itk::StartEvent());
    m_Stop = false;



    ParametersType x,newX;
    DerivativeType gx;
    x.SetSize(NumberOfParameters);
    gx.SetSize(NumberOfParameters);

    x = this->GetInitialPosition();

    this->m_CostFunction->GetDerivative(x,gx);

    gx.normalize();
    gx = -gx;
    unsigned int i =0;
    double cost = 0.0;
    cost = this->m_CostFunction->GetValue(x);
    double mu =0.0;

    double epsilon = 100.0;
    this->SetCurrentPosition(x);
    this->m_CurrentValue = cost;

    //this->SmartSearch(x,gx);
    gx = this->CheckParameters(gx);

    while(i < m_NumberOfIterations &&  epsilon > m_Epsilon)
    {
        //std::cout<<"Iteration : "<<i<<std::endl;
        //mu = this->LinMin(x, gx);
        //mu=0.01;
        mu = this->SearchStep(x,gx);
        newX = x + (mu * gx);
        //std::cout<<"New optimized parameters : "<<newX<<std::endl;

        double newCost = this->m_CostFunction->GetValue(newX);
        //std::cout<<"New Cost : "<<newCost<<std::endl;

        epsilon = newCost-cost;

        //std::cout<<"Epsilon : "<<epsilon<<std::endl;

        if(newCost < cost)
        {

            this->SetCurrentPosition(newX);
            this->m_CurrentValue = newCost;
            cost = newCost;
            //std::cout<<"NEW BEST COST :"<<cost<<std::endl;
            //std::cout<<"NEW BEST PARAMETERS :"<<newX<<std::endl;
        }
        //x = this->GetCurrentPosition();
        x = newX;
        this->m_CostFunction->GetDerivative(x,gx);
        gx.normalize();
        gx = -gx;
        gx = this->CheckParameters(gx);

        //std::cout<<"Derivative : "<<gx<<std::endl;

        i++;
    }


    this->InvokeEvent(itk::EndEvent());
}
//-------------------------------------------------------------------------------------------------
double SmartStepGradientDescentOptimizer::SearchStep(ParametersType _x , DerivativeType _gx)
{
    DerivativeType newGx;
    ParametersType newX;

    _gx = this->CheckParameters(_gx); // if we don't want to optimize all parameters
    double max = std::fabs(_gx.max_value());

    if(max == 0)
    {
        return  m_MinStep;
    }


    double maxmu = m_MaxStep/max;

    double minmu = m_MinStep/max;

    double mu = minmu;
    //double mu = maxmu;

    double step = (maxmu - minmu)/50.0;

    //double initialCost = this->m_CostFunction->GetValue(_x);


    newGx = mu * _gx;
    newX = _x + newGx;
    double cost = this->m_CostFunction->GetValue(newX);
    mu+=step;

    //std::cout<<"    *Gradient : "<<_gx<<std::endl;
    //std::cout<<"    *New Parameters : "<<newX<<std::endl;

    while(mu <= maxmu)
    {
        newGx = mu * _gx;
        newX = _x + newGx;
        double newCost = this->m_CostFunction->GetValue(newX);

        if(newCost > cost)
            break;

        _x = newX;
        this->m_CostFunction->GetDerivative(_x,_gx);
        _gx.normalize();
        _gx = -_gx;

        mu+=step;
    }

    //std::cout<<"    *µ : "<<mu<<std::endl;
    return mu;







}
//-------------------------------------------------------------------------------------------------
SmartStepGradientDescentOptimizer::DerivativeType
SmartStepGradientDescentOptimizer::CheckParameters(DerivativeType _gx)
{
    unsigned int size = m_OptimizedParameters.Size();

    for(unsigned int i=0; i<size; i++)
    {
        _gx[i] = _gx[i]* m_OptimizedParameters[i];
    }
    return _gx;
}
//-------------------------------------------------------------------------------------------------
}
