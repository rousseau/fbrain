/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

14/08/2012
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

#ifndef __BTKHISTOGRAM_CXX__
#define __BTKHISTOGRAM_CXX__

#include "btkHistogram.h"

namespace btk
{    
  Histogram::Histogram()
  {
    m_numberOfSamples = 0;
    m_numberOfBins = 0;
  }
  
  Histogram::~Histogram()
  {
  
  }
  
  void Histogram::SetNumberOfBins(unsigned int numberOfBins)
  {
    m_data.resize(numberOfBins);
    m_normalizedData.resize(numberOfBins);
    m_numberOfBins = m_data.size();
  }
  
  void Histogram::SetSampleQuantification(unsigned int sampleQuantification)
  {
    m_sampleQuantification = sampleQuantification;
  }  
    
  void Histogram::SetLowerBound(float lowerBound)
  {
    m_lowerBound = lowerBound;
  }
  
  void Histogram::SetUpperBound(float upperBound)
  {
    m_upperBound = upperBound;
  }
  
  //This function computes the linear coefficients used to convert bin to values and m_widthOfBin
  void Histogram::Setup()
  {
    m_aCoefficient = (m_numberOfBins-1) * 1.0 / (m_upperBound - m_lowerBound);
    m_bCoefficient = - m_aCoefficient * m_lowerBound;
    m_widthOfBin   = m_upperBound / (m_numberOfBins - 1) ;
  }
  
  
  void Histogram::AddSample(float sample)
  {
    //No boundary checking !
    unsigned int bin = (int) (sample * m_aCoefficient + m_bCoefficient);
    m_data[bin]++;
    m_numberOfSamples++;
  }
  
  
  void Histogram::RemoveSample(float sample)
  {
    //No boundary checking !
    unsigned int bin = (int) (sample * m_aCoefficient + m_bCoefficient);
    m_data[bin]--;
    m_numberOfSamples--;    
  }
  
  
  float Histogram::BinToValue(unsigned int bin)
  {
    return ( bin*(m_upperBound - m_lowerBound)/m_numberOfBins + m_lowerBound ) + m_widthOfBin / 2.0;
  }
  
  
  void Histogram::ComputeCumulativeDistributionFunction()
  {
    m_cumulativeDistributionFunction.resize(m_numberOfBins);
    m_cumulativeDistributionFunction[0] = m_data[0];
    for(unsigned int i=1; i<m_numberOfBins; i++)
      m_cumulativeDistributionFunction[i] = m_data[i] + m_cumulativeDistributionFunction[i-1];
  }
  
  
  void Histogram::ComputeInverseCumulativeDistributionFunction()
  {
    m_inverseCumulativeDistributionFunction.resize(m_numberOfSamples);
    m_inverseCumulativeDistributionFunction[0] = 0;
    unsigned int currentBin = 0;
    
    for(unsigned int i=1; i<m_numberOfSamples+1; i++)
    {
      m_inverseCumulativeDistributionFunction[i] = currentBin;
      for(unsigned int bin=currentBin; bin<m_numberOfBins; bin++)
        if( m_cumulativeDistributionFunction[bin] >= i )
        {
          currentBin = bin;
          m_inverseCumulativeDistributionFunction[i] = bin;
          break;          
        }    
    }

  }
  
  //NORMALIZED VERSIONS -----------------------------------------------------------
  void Histogram::NormalizeData()
  {
    for(unsigned int i=0; i<m_numberOfBins; i++)
      m_normalizedData[i] = m_data[i] / m_numberOfSamples * (m_sampleQuantification-1);
  }

  void Histogram::ComputeNormalizedCumulativeDistributionFunction()
  {
    m_normalizedCumulativeDistributionFunction.resize(m_numberOfBins);
    m_normalizedCumulativeDistributionFunction[0] = m_normalizedData[0];
    for(unsigned int i=1; i<m_numberOfBins; i++)
      m_normalizedCumulativeDistributionFunction[i] = m_normalizedData[i] + m_normalizedCumulativeDistributionFunction[i-1];
  }
  
  
  void Histogram::ComputeNormalizedInverseCumulativeDistributionFunction()
  {
    m_normalizedInverseCumulativeDistributionFunction.resize(m_sampleQuantification);
    m_normalizedInverseCumulativeDistributionFunction[0] = 0;
    unsigned int currentBin = 0;
    
    for(unsigned int i=1; i<m_sampleQuantification; i++)
    {
      m_normalizedInverseCumulativeDistributionFunction[i] = currentBin;
      for(unsigned int bin=currentBin; bin<m_numberOfBins; bin++)
        if( m_normalizedCumulativeDistributionFunction[bin] >= i )
        {
          currentBin = bin;
          m_normalizedInverseCumulativeDistributionFunction[i] = bin;
          break;          
        }    
    }

  }

}
#endif // btkHistogram_CXX