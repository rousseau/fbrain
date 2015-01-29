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

#ifndef btkJointHistogram_H
#define btkJointHistogram_H

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

/* BTK */
#include "btkMacro.h"

namespace btk
{

class JointHistogram
{

public:

  vnl_matrix< double >  m_Data;

  JointHistogram();
  ~JointHistogram();

  btkSetMacro(Ax,double);
  btkGetMacro(Ax,double);
  btkSetMacro(Bx,double);
  btkGetMacro(Bx,double);
  btkSetMacro(Ay,double);
  btkGetMacro(Ay,double);
  btkSetMacro(By,double);
  btkGetMacro(By,double);
  btkSetMacro(NumberOfSamples,double);
  btkGetMacro(NumberOfSamples,double);
  btkSetMacro(Data,vnl_matrix<double>);
  btkGetMacro(Data,vnl_matrix<double>);

  void SetNumberOfBins(unsigned int nx, unsigned int ny);
  void ClearJointHistogram();
  unsigned int GetNumberOfBinsX();
  unsigned int GetNumberOfBinsY();
  void AddSample(double x, double y, double w);

  double MutualInformation();
  double NormalizedMutualInformation();

  double EntropyX();
  double EntropyY();
  double JointEntropy();


private:

  //linear mapping parameters to convert intensities to bin
  double m_Ax, m_Bx, m_Ay, m_By;
  double m_NumberOfSamples;

  
};
}



#endif

