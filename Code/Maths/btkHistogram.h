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


#ifndef btkHistogram_H
#define btkHistogram_H


#include "string"
#include "iomanip"
#include "sstream"
#include "fstream"
#include "iostream"

namespace btk
{
class Histogram
{
  public:
  Histogram();
  ~Histogram();
  
  void SetNumberOfBins(unsigned int numberOfBins);
  void SetSampleQuantification(unsigned int sampleQuantification);
  void SetLowerBound(float lowerBound);
  void SetUpperBound(float upperBound);
  
  void Setup();
  
  void AddSample(float sample);
  void RemoveSample(float sample);
  float BinToValue(unsigned int bin);
  void AddWeightedSample(float sample, float weight);
  void RemoveWeightedSample(float sample, float weight);
  
  void ComputeCumulativeDistributionFunction();
  void ComputeInverseCumulativeDistributionFunction();
  
  void NormalizeData();
  void ComputeNormalizedCumulativeDistributionFunction();
  void ComputeNormalizedInverseCumulativeDistributionFunction();
  
  void SaveHistogram(const std::string & filename);
  void SaveNormalizedHistogram(const std::string & filename);
  void SaveCumulativeDistributionFunction(const std::string & filename);
  void SaveNormalizedCumulativeDistributionFunction(const std::string & filename);
  void SaveInverseCumulativeDistributionFunction(const std::string & filename);
  void SaveNormalizedInverseCumulativeDistributionFunction(const std::string & filename);  
    
  std::vector<float> m_data;  
  std::vector<float> m_cumulativeDistributionFunction; //not normalized
  std::vector<float> m_inverseCumulativeDistributionFunction; //not normalized
  
  std::vector<float> m_normalizedData;  
  std::vector<float> m_normalizedCumulativeDistributionFunction;
  std::vector<float> m_normalizedInverseCumulativeDistributionFunction; 
  
  
  unsigned int m_numberOfBins;
  unsigned int m_sampleQuantification;
  float m_numberOfSamples;
  float m_lowerBound;
  float m_upperBound;
  float m_aCoefficient;
  float m_bCoefficient;
  float m_widthOfBin;
  
  
  
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkHistogram.cxx"
#endif


#endif

