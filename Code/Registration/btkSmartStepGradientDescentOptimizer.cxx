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
    m_MinStep(0.05),m_NumberOfIterations(1000),m_Epsilon(10e-6),m_OptimizeAllParameters(true),
    m_VerboseMode(false),m_UseBounds(false),m_Samples(100.0)

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

    //std::cout<<"Optimizer number of parameters : "<<NumberOfParameters<<std::endl;

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

    m_X.SetSize(NumberOfParameters);
    m_gX.SetSize(NumberOfParameters);

    m_X = this->GetInitialPosition();



    this->m_CostFunction->GetDerivative(m_X,gx);


    gx.normalize();
    m_gX = -gx; // inverse the gradient !
    unsigned int i = 0;
    double cost = 0.0;
    cost = this->m_CostFunction->GetValue(m_X);
    double mu = 0.0;

    double epsilon = 100.0; // to avoid random initialization
    this->SetCurrentPosition(m_X);
    this->m_CurrentValue = cost;

    m_gX = this->CheckParameters(m_gX);

    while(i < m_NumberOfIterations &&  epsilon > m_Epsilon)
    {

        mu = this->SearchStep(m_X,m_gX); // look for the step
        x = m_X + (mu * m_gX); // new parameters

        if(m_UseBounds)
        {
            if(!this->IsInsideBounds(x)) // we looking if parameters are inside bounds
            {
                //if not we initialize parameter which is outside with the opposite bounds (for looping)
                newX = this->ReverseBounds(x);
                this->m_CostFunction->GetDerivative(newX,gx);
                gx.normalize();
                m_gX = -gx;
                mu = this->SearchStep(newX,m_gX);
                x = newX + (mu * m_gX);
            }
        }



        double newCost = this->m_CostFunction->GetValue(x); // cost of the new parameters


        //epsilon = std::fabs(newCost-cost); /**** BUG ****/
        epsilon = cost - newCost;

        //FIXME : epsilon should not be negative

        if(m_VerboseMode)
        {
            std::cout<<"******************************"<<std::endl;
            std::cout<<"    old X : "<<m_X<<std::endl;
            std::cout<<"    gradient : "<<m_gX<<std::endl;
            std::cout<<"    µ : "<<mu<<std::endl;
            std::cout<<"    new X : "<<x<<std::endl;
            std::cout<<"    old Cost : "<<cost<<std::endl;
            std::cout<<"    new Cost : "<<newCost<<std::endl;
            std::cout<<"    epsilon : "<<epsilon<<std::endl;
            std::cout<<"******************************"<<std::endl;
        }

        if(newCost < cost)
        {
            this->SetCurrentPosition(x);
            this->m_CurrentValue = newCost;
            cost = newCost;

        }
        //----------------------------------
        else
        {
            this->SetCurrentPosition(m_X);
            this->m_CurrentValue = cost;
        }
        //--------------------------------
        
		//Should be in the previous brackets? -----------------
        m_X = x;
        this->m_CostFunction->GetDerivative(m_X,gx);
        gx.normalize();
        m_gX = -gx;
        m_gX = this->CheckParameters(m_gX);

        i++;
    }
    if(m_VerboseMode)
    {
    	std::cout<<"Stopping optimization. Number of iterations : "<<i<<std::endl;
    	std::cout<<"--------------------------------------------------"<<std::endl;
    }


    this->InvokeEvent(itk::EndEvent());
}
//-------------------------------------------------------------------------------------------------
double SmartStepGradientDescentOptimizer::SearchStep(ParametersType _x , DerivativeType _gx)
{
    DerivativeType newGx;
    ParametersType newX;

    _gx = this->CheckParameters(_gx); // if we don't want to optimize all parameters
    double max = 0.0;

    for(unsigned int i = 0; i<_gx.size(); i++)
    {
        if(std::fabs(_gx[i]) > max)
        {
            max = std::fabs(_gx[i]); //looking for the direction
        }
    }


    if(max == 0)
    {
        return  0.0; // avoid division by 0
    }

    //Should we divide by max ? (if normalized max can't be > 1)
    double maxmu = m_MaxStep/max; // compute a max step

    double minmu = m_MinStep/max; // compute a min step

    double mu = minmu; // use the min step at the begining
    //double mu = maxmu;

    double step = (maxmu - minmu)/m_Samples; // sample between minmu and maxmu

    //double initialCost = this->m_CostFunction->GetValue(_x);

    double lastmu = 0.0;
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
        {
            //if newCost is > cost we stop (mu is then the previous value)
            mu = lastmu;
            break;
        }

        _x = newX;
        this->m_CostFunction->GetDerivative(_x,_gx);
        _gx.normalize();
        _gx = -_gx;
        lastmu = mu;
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
bool SmartStepGradientDescentOptimizer::IsInsideBounds(ParametersType _xIn)
{
    unsigned int size = _xIn.Size();
    bool isInside = true;

    for(unsigned int i = 0; i <size; i++)
    {
        if(_xIn[i] > m_MaxBounds[i])
        {
            isInside = false;
        }
        if(_xIn[i] < m_MinBounds[i])
        {
            isInside = false;
        }

    }

    return isInside;
}

//-------------------------------------------------------------------------------------------------
SmartStepGradientDescentOptimizer::ParametersType
SmartStepGradientDescentOptimizer::ReverseBounds(ParametersType _xIn)
{
    unsigned int size = _xIn.Size();
    ParametersType xOut(size);

    for(unsigned int i = 0; i <size; i++)
    {
        xOut[i] = _xIn[i];

        if(_xIn[i] > m_MaxBounds[i])
        {
            xOut[i] = m_MinBounds[i];
        }
        if(_xIn[i] < m_MinBounds[i])
        {
            xOut[i] = m_MaxBounds[i];
        }

    }

    return xOut;
}
}
