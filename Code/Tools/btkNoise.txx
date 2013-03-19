/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 january 2013
< rousseau at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include "btkNoise.h"


namespace btk
{

template<typename T>
void Noise<T>::ComputeGlobalSigma(itkTImagePointer & image, itkTImagePointer & mask)
{
  //std::cout<<"Compute global sigma over the masked image\n";
  //std::cout<<"sigma is estimated using pseudo-residual technique (which involves Laplacian filtering)\n";
  typename itkCastFilterType::Pointer castFilter = itkCastFilterType::New();
  castFilter->SetInput(image);
  castFilter->Update();
  
  typename itkLaplacianFilter::Pointer laplacian = itkLaplacianFilter::New();  
  laplacian->SetInput( castFilter->GetOutput() ); // NOTE: input image type must be double or float
  laplacian->Update();
    
  itkFloatImagePointer laplacianImage = laplacian->GetOutput();
  std::vector<float> vecei;
  
  itkFloatIterator itLaplacian( laplacianImage, laplacianImage->GetLargestPossibleRegion() );
  itkTIterator itMask( mask, mask->GetLargestPossibleRegion());
  
  for(itLaplacian.GoToBegin(), itMask.GoToBegin(); !itLaplacian.IsAtEnd(); ++itLaplacian, ++itMask){
    if( itMask.Get() > 0 )
      vecei.push_back( itLaplacian.Get() / sqrt(42.0) ); // 42 = 6*6+1+1+1+1+1+1 !      
  }
   
  //Estimation of sigma with MAD
  std::sort(vecei.begin(), vecei.end());
  float med = vecei[(int)(vecei.size()/2)];
  for(unsigned int i=0; i<vecei.size(); i++)
    vecei[i] = fabs(vecei[i] - med);
  std::sort(vecei.begin(), vecei.end());
    
  double sigma = 1.4826 * vecei[(int)(vecei.size()/2)];
  //std::cout<<"Estimated sigma:"<<sigma<<"\n";
  m_Sigma2Image->FillBuffer(sigma*sigma);  
}



} // namespace btk

