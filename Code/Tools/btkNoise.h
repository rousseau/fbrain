/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 january 2013
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

#ifndef BTK_NOISE_H
#define BTK_NOISE_H

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkLaplacianImageFilter.h"
#include "itkCastImageFilter.h"

#include "btkMacro.h"

namespace btk
{

  /**
  * @class Noise
  * @brief Noise estimation
  * @author François Rousseau
  * @ingroup Tools
  */
  template<typename T>
  class Noise
  {
    public:

      typedef typename itk::Image< T, 3>                          itkTImage;
      typedef typename itkTImage::Pointer                         itkTImagePointer;
      typedef typename itk::ImageRegionIterator< itkTImage >      itkTIterator;
      typedef typename itk::ImageRegionConstIterator< itkTImage > itkConstIterator;
	  typedef typename itkTImage::SpacingType                     itkTSpacing;  
	      
  	  typedef typename itk::Image< float, 3>                               itkFloatImage;
      typedef typename itkFloatImage::Pointer                              itkFloatImagePointer;
      typedef typename itk::ImageRegionIterator< itkFloatImage >           itkFloatIterator;

	  typedef typename itk::CastImageFilter< itkTImage, itkFloatImage >        itkCastFilterType;
	  typedef typename itk::LaplacianImageFilter<itkFloatImage,itkFloatImage>  itkLaplacianFilter;
  
      btkGetMacro(Sigma2Image, itkFloatImagePointer);
      btkSetMacro(Sigma2Image, itkFloatImagePointer);

      void ComputeGlobalSigma(itkTImagePointer & image, itkTImagePointer & mask);

    private:
    
      itkFloatImagePointer m_Sigma2Image; //variance of the noise
  
  };

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkNoise.txx"
#endif

#endif // BTK_NOISE_H
