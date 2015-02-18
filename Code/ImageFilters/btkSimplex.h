/*
 Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

 14 january 2015
 rousseau@unistra.fr

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
 */

#ifndef BTK_SIMPLEX_H
#define BTK_SIMPLEX_H

#include <cmath>

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

/* Itk includes */
//#include "itkVector.h"

#include "btkMacro.h"

namespace btk
{
class Simplex
{
public:

    vnl_vector< double > m_FunctionValues;
    vnl_matrix< double > m_CurrentSimplex;
    vnl_vector< double > m_ConvergenceToleranceForEachParameter;

    btkSetMacro(ConvergenceTolerance,double);
    btkGetMacro(ConvergenceTolerance,double);
    btkSetMacro(ConvergenceToleranceForEachParameter,vnl_vector<double>);
    btkGetMacro(ConvergenceToleranceForEachParameter,vnl_vector<double>);
    btkSetMacro(Dimension,int);
    btkGetMacro(Dimension,int);
    btkSetMacro(CurrentSimplex,vnl_matrix<double>);
    btkSetMacro(NumberOfFunctionsEvaluations,int);
    btkGetMacro(NumberOfFunctionsEvaluations,int);
    btkSetMacro(MinimumValue,double);

    Simplex()
    {
        this->SetNumberOfFunctionsEvaluations(1000);
        this->SetConvergenceTolerance(1e-3);
    }

    template <class T> vnl_vector< double > minimize(vnl_vector< double > & point, double displacement, T & f)
    {
        vnl_vector< double > displacements(point.size());
        displacements.fill(displacement);
        return minimize(point, displacements, f);
    }

    template <class T> vnl_vector< double > minimize(vnl_vector< double > & point, vnl_vector< double > & displacements, T & f)
    {
        //std::cout<<"Run Simplex Optimization"<<std::endl;
        //std::cout<<"Initial displacements : "<<displacements<<std::endl;

        this->SetDimension(point.size());
        vnl_matrix<double> initialSimplex(this->GetDimension()+1,this->GetDimension());
        //Initialization of the simplex
        for(unsigned int i=0; i < this->GetDimension()+1; i++)
        {
            for(unsigned int j=0; j < this->GetDimension(); j++)
                initialSimplex(i,j) = point(j);
            if(i > 0)
                initialSimplex(i,i-1) += displacements(i-1);
        }
        return minimize(initialSimplex,f);
    }

    template <class T> vnl_vector< double > minimize(vnl_matrix<double> & initialSimplex, T & f)
    {
        /*
        std::cout<<"Initial point : "<<initialSimplex.get_row(0)<<std::endl;
        if(m_ConvergenceToleranceForEachParameter.size() != this->GetDimension())
            std::cout<<"Global tolerance (on cost function value) : "<<this->GetConvergenceTolerance()<<std::endl;
        else
            std::cout<<"Tolerance (on each parameter value) : "<<this->GetConvergenceToleranceForEachParameter()<<std::endl;
        */

        vnl_vector< double > sum( this->GetDimension() );
        vnl_vector< double > min( this->GetDimension() );
        vnl_vector< double > x( this->GetDimension() );

        this->SetCurrentSimplex(initialSimplex);
        m_FunctionValues.set_size( this->GetDimension() + 1);

        //Compute function values at the vertices of the simplex
        for(unsigned int i=0; i < this->GetDimension()+1; i++)
        {
            x = m_CurrentSimplex.get_row(i);
            m_FunctionValues(i) = f(x);
        }
        //Compute sum of each column of the current simplex
        for(unsigned int i=0; i < this->GetDimension(); i++)
            sum(i) = m_CurrentSimplex.get_column(i).sum();

        //main loop
        unsigned int currentIteration = 0;
        while( currentIteration < this->GetNumberOfFunctionsEvaluations() ) //c est con, il vaut mieux le nombre d'iterations
        {
            //Determine highest, next highest and lowest function values.
            int highestValueIndex = 1;
            int nextHighestValueIndex = 0;
            int lowestValueIndex = 0;
            if(m_FunctionValues(1) < m_FunctionValues(0))
            {
                highestValueIndex = 0;
                nextHighestValueIndex = 1;
            }
            for(unsigned int i=0; i < this->GetDimension()+1; i++)
            {
                if(m_FunctionValues(i) < m_FunctionValues(lowestValueIndex))
                    lowestValueIndex = i;
                if(m_FunctionValues(i) > m_FunctionValues(highestValueIndex))
                {
                    nextHighestValueIndex = highestValueIndex;
                    highestValueIndex = i;
                }
                else if (m_FunctionValues(i) > m_FunctionValues(nextHighestValueIndex) && i != highestValueIndex)
                    nextHighestValueIndex = i;
            }

            //Compute the current fractional range
            bool convergenceReached = false;
            //If no tolerance for each parameter has been properly set, then we use a global tolerance
            if(m_ConvergenceToleranceForEachParameter.size() != this->GetDimension())
            {
                double currentTolerance = 2.0 * std::abs(m_FunctionValues(highestValueIndex) - m_FunctionValues(lowestValueIndex)) / ( std::abs(m_FunctionValues(highestValueIndex)) +std::abs(m_FunctionValues(lowestValueIndex)) + 1.0e-10 );
                if (currentTolerance < this->GetConvergenceTolerance() )
                    convergenceReached = true;
            }
            else
            {
                vnl_vector< double > best( this->GetDimension() );
                vnl_vector< double > worst( this->GetDimension() );
                best = m_CurrentSimplex.get_row(lowestValueIndex);
                worst= m_CurrentSimplex.get_row(highestValueIndex);

                //std::cout<<"difference between best and worst : "<<best-worst<<std::endl;

                convergenceReached = true;
                for(unsigned int i = 0; i < this->GetDimension(); i++)
                {
                    double currentTolerance = std::abs(best(i) - worst(i));
                    if (currentTolerance > this->m_ConvergenceToleranceForEachParameter(i) )
                    {
                        convergenceReached = false;
                        break;
                    }
                }
            }
            //std::cout<<"current best estimate : ";
            //std::cout<<m_CurrentSimplex.get_row(lowestValueIndex)<<std::endl;
            //std::cout<<"current cost function value : "<<m_FunctionValues(lowestValueIndex)<<std::endl;
            //std::cout<<"currentTolerance : "<<currentTolerance<<std::endl;
            //if the minimum has been reached ... (should add also a test on the maximum iterations, and it there is no change of simplex shape)
            if ( convergenceReached == true )
            {
                this->SetMinimumValue( m_FunctionValues(lowestValueIndex) );
                min = m_CurrentSimplex.get_row(lowestValueIndex);
                return min;
            }

            currentIteration += 2;
            //Reflect the simplex (factor -1.0)
            double newValue = this->TryAndReplace(sum,highestValueIndex,-1.0,f);
            if(newValue <= m_FunctionValues(lowestValueIndex) )
                //Better estimate, so continue into this direction (use a factor 2.0)
                this->TryAndReplace(sum,highestValueIndex,2.0,f);
            else if(newValue >= m_FunctionValues(nextHighestValueIndex) )
            {
                double previousValue = m_FunctionValues(highestValueIndex);
                //Worse estimate, so try another factor (0.5)
                newValue = this->TryAndReplace(sum,highestValueIndex,0.5,f);
                if(newValue >= previousValue)
                {
                    //Worse estimater again, build another simplex around the best point
                    for(unsigned int i=0; i < this->GetDimension()+1; i++)
                    {
                        if( i != lowestValueIndex)
                        {
                            for(unsigned int j=0; j < this->GetDimension(); j++)
                                m_CurrentSimplex(i,j) = 0.5*(m_CurrentSimplex(i,j) + m_CurrentSimplex(lowestValueIndex,j));
                            m_FunctionValues(i) = f(m_CurrentSimplex.get_row(i));
                        }
                    }
                    for(unsigned int j=0; j < this->GetDimension(); j++)
                        sum(j) = m_CurrentSimplex.get_column(j).sum();
                }
            }
        }
    }

    template <class T> double TryAndReplace(vnl_vector<double> & sum, int index, double factor, T&f)
    {
        vnl_vector<double> newPoint(this->GetDimension());
        double factor1 = (1.0 - factor) / this->GetDimension();
        double factor2 = factor1-factor;
        for(unsigned int i=0; i < this->GetDimension(); i++)
            newPoint(i) = sum(i)*factor1 - m_CurrentSimplex(index,i)*factor2;
        double newValue = f(newPoint);
        //If better, update simplex
        if(newValue < m_FunctionValues(index))
        {
            m_FunctionValues(index) = newValue;
            for(unsigned int i=0; i < this->GetDimension(); i++)
            {
                m_CurrentSimplex(index,i) = newPoint(i);
                sum(i) = m_CurrentSimplex.get_column(i).sum();
            }
        }
        return newValue;
    }


private:
    double               m_ConvergenceTolerance;
    int                  m_NumberOfFunctionsEvaluations;
    int                  m_Dimension;
    double               m_MinimumValue;


};


} // namespace btk

#endif
