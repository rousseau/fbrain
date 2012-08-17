/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

16 august 2012
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

  //Principle of the algorithm: apply a midway histogram equalization to all input images 
  //See these papers: 
  //Cox et al. Dynamic histogram warping of image pairs for constant image brightness, ICIP 1995
  //Delon, Midway Image Equalization, Journal of Mathematical Imaging and Vision (JMIV), vol.21, no.2, pp.119-134, 2004.
  
  //Algorithm:
  //Rescale all images between 0 and 1
  //Compute the histograms of each image
  //Compute the cumulative distribution functions (cdf) H_i
  //Compute the inverse cdf
  //Average all the inverse cdf 
  //Inverse this mean inverse cdf (called H_new)
  //Inverse H_new
  //Apply to all images H_new^-1 o H_i
  
  
#ifndef btkMidwayImageEqualization_H
#define btkMidwayImageEqualization_H


#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"

#include "btkHistogram.h" 

#include "string"
#include "iomanip"
#include "sstream"
#include "fstream"

 
namespace btk
{

template <typename TPixelType>
class MidwayImageEqualization
{
 public:
  
  typedef typename itk::Image< TPixelType, 3> itkTImage;
  typedef typename itkTImage::Pointer itkTPointer;
  typedef typename itk::ImageDuplicator< itkTImage > itkTDuplicator;
  typedef typename itk::ImageRegionIterator< itkTImage > itkTIterator;

  unsigned int m_numberOfBins;
  unsigned int m_sampleQuantification;
  
  void SetNumberOfBins(unsigned int n);
  void SetSampleQuantification(unsigned int n);
  
  void Do(std::vector<itkTPointer> inputImages, std::vector<itkTPointer> maskImages, std::vector<itkTPointer> outputImages);

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkMidwayImageEqualization.txx"
#endif


#endif

