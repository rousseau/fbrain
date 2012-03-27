/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 09/12/2011
  Author(s): François Rousseau (rousseau@unistra.fr)

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

#ifndef __btkSuperResolutionManager_h
#define __btkSuperResolutionManager_h

#include "itkImage.h"
#include "itkImageDuplicator.h"

#include "btkSuperResolutionDataManager.h"
#include "btkSuperResolutionTools.h"


class SuperResolutionManager
{
public:
  typedef float PixelType;
  typedef itk::Image<PixelType, 3>           itkImage;
  typedef itkImage::Pointer                  itkPointer;  
  typedef itk::ImageDuplicator< itkImage >   itkDuplicator;
  
  SuperResolutionDataManager    data;
  SuperResolutionTools          tool;
  
  void Initialize();
  void SimulateLRImages();
  void IteratedBackProjection(int & loops, int & nlm, float & beta, int & medianIBP);
  
};

void SuperResolutionManager::Initialize()
{
  tool.InitializePSF(data);
  tool.CreateMaskHRImage(data);
  tool.HComputation(data);
}

void SuperResolutionManager::SimulateLRImages()
{
  tool.SimulateLRImages(data);
}

void SuperResolutionManager::IteratedBackProjection(int & loops, int & nlm, float & beta, int & medianIBP)
{
  double e = 0.1; //current change between two consecutive estimate.
  int i = 0;
  
  while( (e>0) && (i<loops) ){
    e = tool.IteratedBackProjection(data, nlm, beta, medianIBP);
    i++;
    std::cout<<"Loop "<<i+1<<", current e: "<<e<<"\n";
  }
}
#endif
