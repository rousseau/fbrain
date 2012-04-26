/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 22/03/2012
  Author(s): Schweitzer Marc (marc.schweitzer@unistra.fr)

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
#ifndef __BTK_HIGHRESOLUTIONIBPFILTER_H__
#define __BTK_HIGHRESOLUTIONIBPFILTER_H__




/* BTK */
#include "btkMacro.h"
#include "btkHighResolutionReconstructionFilter.h"
#include "btkCreateHRMaskFilter.h"
#include "btkSimulateLRImageFilter.h"



/* OTHERS */
#include "iostream"

namespace btk
{
class HighResolutionIBPFilter : public HighResolutionReconstructionFilter
{


public:



    HighResolutionIBPFilter();
    ~HighResolutionIBPFilter();

    void Initialize();

    virtual void Update();

    btkGetMacro(Nlm,int);
    btkSetMacro(Nlm,int);

    btkGetMacro(Beta,float);
    btkSetMacro(Beta,float);

    btkGetMacro(MedianIBP,int);
    btkSetMacro(MedianIBP,int);


protected:

    void HComputation();
    void InitializePSF();
    void UpdateX();
    double ComputeIterativeBackProjection( int & nlm, float & beta, int & medianIBP);

private:

    int m_Nlm;
    float m_Beta;
    int m_MedianIBP;


    CreateHRMaskFilter* m_HRMaskFilter;
    SimulateLRImageFilter* m_SimuLRImagesFilter;



};
}

#endif
