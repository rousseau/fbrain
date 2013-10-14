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
  std::cout<<"Compute global sigma over the masked image\n";
  std::cout<<"sigma is estimated using pseudo-residual technique (which involves Laplacian filtering)\n";

  std::vector<double> vecei;

  //IMPORTANT NOTE:
  //Using ITK Laplacian does not provide correct result !!!!
  //Manual computation or ITK convolution is correct.


  //------------------------------------------------------------------------
  //             ITK LAPLACIAN (not correct)
  /*
  typename itkCastFilterType::Pointer castFilter = itkCastFilterType::New();
  castFilter->SetInput(image);
  castFilter->Update();

  typename itkLaplacianFilter::Pointer laplacian = itkLaplacianFilter::New();
  laplacian->SetInput( castFilter->GetOutput() ); // NOTE: input image type must be double or float
  laplacian->Update();

  itkFloatImagePointer laplacianImage = laplacian->GetOutput();
  */
  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //             ITK CONVOLUTION

  typename itkCastFilterType::Pointer castFilter = itkCastFilterType::New();
  castFilter->SetInput(image);
  castFilter->Update();

  itkFloatImagePointer kernel = itkFloatImage::New();
  itkFloatImage::IndexType kernelStart;
  kernelStart.Fill(0);

  itkFloatImage::SizeType kernelSize;
  kernelSize.Fill(3);

  itkFloatImage::RegionType kernelRegion;
  kernelRegion.SetSize(kernelSize);
  kernelRegion.SetIndex(kernelStart);

  kernel->SetRegions(kernelRegion);
  kernel->Allocate();
  kernel->FillBuffer(0.0);

  typename itkFloatImage::IndexType p;
  p[0] = 1; p[1] = 1; p[2] = 1; kernel->SetPixel(p,6);
  p[0] = 0; p[1] = 1; p[2] = 1; kernel->SetPixel(p,-1);
  p[0] = 2; p[1] = 1; p[2] = 1; kernel->SetPixel(p,-1);
  p[0] = 1; p[1] = 0; p[2] = 1; kernel->SetPixel(p,-1);
  p[0] = 1; p[1] = 2; p[2] = 1; kernel->SetPixel(p,-1);
  p[0] = 1; p[1] = 1; p[2] = 0; kernel->SetPixel(p,-1);
  p[0] = 1; p[1] = 1; p[2] = 2; kernel->SetPixel(p,-1);


  itk::ConvolutionImageFilter<itkFloatImage>::Pointer convolutionFilter = itk::ConvolutionImageFilter<itkFloatImage>::New();
  convolutionFilter->SetInput( castFilter->GetOutput() );
  convolutionFilter->SetKernelImage(kernel);
  convolutionFilter->Update();
  itkFloatImagePointer laplacianImage = convolutionFilter->GetOutput();

  itkFloatIterator itLaplacian( laplacianImage, laplacianImage->GetLargestPossibleRegion() );
  itkTIterator itMask( mask, mask->GetLargestPossibleRegion());

  for(itLaplacian.GoToBegin(), itMask.GoToBegin(); !itLaplacian.IsAtEnd(); ++itLaplacian, ++itMask){
    if( itMask.Get() > 0 )
      vecei.push_back( itLaplacian.Get() / (double)sqrt(42.0) ); // 42 = 6*6+1+1+1+1+1+1 !
  }

  //------------------------------------------------------------------------
  //------------------------------------------------------------------------
  //             MANUAL COMPUTATION

  /*
  typename itkTImage::SizeType    size    = image->GetLargestPossibleRegion().GetSize();

  for(int z=1; z < (int)size[2]-1; z++)
  for(int y=1; y < (int)size[1]-1; y++)
  for(int x=1; x < (int)size[0]-1; x++)
  {
    typename itkTImage::IndexType p;
    p[0] = x;
    p[1] = y;
    p[2] = z;
    if( mask->GetPixel(p) > 0 )
    {

        double value = 6*image->GetPixel(p);

        p[0] = x+1;	  p[1] = y;	  p[2] = z;
        value = value - image->GetPixel(p);
        p[0] = x-1;	  p[1] = y;	  p[2] = z;
        value = value - image->GetPixel(p);
        p[0] = x;	    p[1] = y+1;	p[2] = z;
        value = value - image->GetPixel(p);
        p[0] = x;  	  p[1] = y-1;	p[2] = z;
        value = value - image->GetPixel(p);
        p[0] = x;  	  p[1] = y;	  p[2] = z+1;
        value = value - image->GetPixel(p);
        p[0] = x;  	  p[1] = y;	  p[2] = z-1;
        value = value - image->GetPixel(p);

        vecei.push_back( value / (double)sqrt(42.0) );

    }
  }
  */

  //------------------------------------------------------------------------

  //MEAN AND STANDARD DEVIATION USING STANDARD FORMULA
  /*
  double stddev = 0;
  double mean = 0;
  for(unsigned int i=0; i<vecei.size(); i++)
      mean += vecei[i];
  std::cout<<"mean = "<<mean/vecei.size()<<std::endl;
  mean /= vecei.size();

  for(unsigned int i=0; i<vecei.size(); i++)
    stddev +=  (vecei[i]-mean)*(vecei[i]-mean);
  std::cout<<"std dev = "<<stddev/vecei.size()<<std::endl;
  */

  //Estimation of sigma with MAD
  std::sort(vecei.begin(), vecei.end());
  double med = vecei[(int)(vecei.size()/2)];
  for(unsigned int i=0; i<vecei.size(); i++)
    vecei[i] = fabs(vecei[i] - med);
  std::sort(vecei.begin(), vecei.end());
    
  double sigma = 1.4826 * vecei[(int)(vecei.size()/2)];
  std::cout<<"Estimated sigma:"<<sigma<<"\n";
  m_Sigma2Image->FillBuffer(sigma*sigma);
  SetGlobalSigma2(sigma*sigma);
}



} // namespace btk

