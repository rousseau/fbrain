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
#include "vector"

/* BTK */
#include "btkMacro.h"

namespace btk
{
/**
 * @class Histogram
 * @brief The Histogram class
 * @author François Rousseau
 * @ingroup Maths
 */
class Histogram
{

public:
    /**
   * @brief Histogram
   */
  Histogram();
  /**
    * @brief destructor
    */
  ~Histogram();
  

  /**
   * @brief Setup
   */
  void Setup();
  /**
   * @brief AddSample
   * @param sample
   */
  void AddSample(float sample);
  /**
   * @brief RemoveSample
   * @param sample
   */
  void RemoveSample(float sample);
  /**
   * @brief BinToValue
   * @param bin
   * @return
   */
  float BinToValue(unsigned int bin);
  /**
   * @brief AddWeightedSample
   * @param sample
   * @param weight
   */
  void AddWeightedSample(float sample, float weight);
  /**
   * @brief RemoveWeightedSample
   * @param sample
   * @param weight
   */
  void RemoveWeightedSample(float sample, float weight);
  /**
   * @brief ComputeCumulativeDistributionFunction
   */
  void ComputeCumulativeDistributionFunction();
  /**
   * @brief ComputeInverseCumulativeDistributionFunction
   */
  void ComputeInverseCumulativeDistributionFunction();
  /**
   * @brief NormalizeData
   */
  void NormalizeData();
  /**
   * @brief ComputeNormalizedCumulativeDistributionFunction
   */
  void ComputeNormalizedCumulativeDistributionFunction();
  /**
   * @brief ComputeNormalizedInverseCumulativeDistributionFunction
   */
  void ComputeNormalizedInverseCumulativeDistributionFunction();
  /**
   * @brief SaveHistogram
   * @param filename
   */
  void SaveHistogram(const std::string & filename);
  /**
   * @brief SaveNormalizedHistogram
   * @param filename
   */
  void SaveNormalizedHistogram(const std::string & filename);
  /**
   * @brief SaveCumulativeDistributionFunction
   * @param filename
   */
  void SaveCumulativeDistributionFunction(const std::string & filename);
  /**
   * @brief SaveNormalizedCumulativeDistributionFunction
   * @param filename
   */
  void SaveNormalizedCumulativeDistributionFunction(const std::string & filename);
  /**
   * @brief SaveInverseCumulativeDistributionFunction
   * @param filename
   */
  void SaveInverseCumulativeDistributionFunction(const std::string & filename);/**
   * @brief SaveNormalizedInverseCumulativeDistributionFunction
   * @param filename
   */
  void SaveNormalizedInverseCumulativeDistributionFunction(const std::string & filename);  

  // GET / SET

  /**
   * @brief btkGetMacro
   */
  btkGetMacro(NormalizedInverseCumulativeDistributionFunction,std::vector<float>);
  /**
   * @brief btkGetMacro
   */
  btkGetMacro(NormalizedCumulativeDistributionFunction,std::vector<float>);
  /**
   * @brief btkGetMacro
   */
  btkGetMacro(ACoefficient,float);
  /**
   * @brief btkGetMacro
   */
  btkGetMacro(BCoefficient,float);
  /**
   * @brief SetNumberOfBins
   * @param numberOfBins
   */
  void SetNumberOfBins(unsigned int numberOfBins);
  /**
   * @brief btkSetMacro
   */
  btkSetMacro(SampleQuantification,unsigned int);
  /**
   * @brief btkSetMacro
   */
  btkSetMacro(LowerBound,float);
  /**
   * @brief btkSetMacro
   */
  btkSetMacro(UpperBound,float);


private:
    
  std::vector<unsigned int> m_Data;
  std::vector<float> m_CumulativeDistributionFunction; //not normalized
  std::vector<float> m_InverseCumulativeDistributionFunction; //not normalized
  
  std::vector<float> m_NormalizedData;
  std::vector<float> m_NormalizedCumulativeDistributionFunction;
  std::vector<float> m_NormalizedInverseCumulativeDistributionFunction;
  
  
  unsigned int m_NumberOfBins;
  unsigned int m_SampleQuantification;
  unsigned int m_NumberOfSamples;
  float m_LowerBound;
  float m_UpperBound;
  float m_ACoefficient;
  float m_BCoefficient;
  float m_WidthOfBin;
  
  
  
};
}



#endif

