/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 15/10/2013
  Author(s): François Rousseau
  
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

#include "btkPandoraBoxImageFilters.h"

namespace btk
{

void PandoraBoxImageFilters::ProbabilityImageNormalization(std::vector< itkFloatImagePointer > & inputImages, std::vector< itkFloatImagePointer > & outputImages)
{
  int x,y,z;
  itkFloatImage::SizeType    size    = inputImages[0]->GetLargestPossibleRegion().GetSize();

  #pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)size[2]; z++)
  for(y=0; y < (int)size[1]; y++)
  for(x=0; x < (int)size[0]; x++)
  {
    itkFloatImage::IndexType index;
    index[0] = x;
    index[1] = y;
    index[2] = z;

    double sumProba = 0;
    for(unsigned int l = 0; l < inputImages.size(); l++)
      sumProba += inputImages[l]->GetPixel( index );

    if(sumProba > 0)
      for(unsigned int l = 0; l < inputImages.size(); l++)
        outputImages[l]->SetPixel(index, inputImages[l]->GetPixel(index) / sumProba);
    else
        for(unsigned int l = 0; l < inputImages.size(); l++)
            outputImages[l]->SetPixel(index, 0);
    }
}

void PandoraBoxImageFilters::GetLabelWithMaxProbabilityImage(std::vector< itkFloatImagePointer > & inputImages, itkShortImagePointer & outputImage)
{
  int x,y,z;
  itkFloatImage::SizeType    size    = inputImages[0]->GetLargestPossibleRegion().GetSize();

  #pragma omp parallel for private(x,y,z) schedule(dynamic)
  for(z=0; z < (int)size[2]; z++)
  for(y=0; y < (int)size[1]; y++)
  for(x=0; x < (int)size[0]; x++)
  {
    itkShortImage::IndexType index;
    index[0] = x;
    index[1] = y;
    index[2] = z;

    float maxProba = 0;
    short maxLabel = 0;
    for(unsigned int l = 0; l < inputImages.size(); l++)
      if( maxProba < inputImages[l]->GetPixel( index ) )
      {
        maxProba = inputImages[l]->GetPixel( index );
        maxLabel = l+1;
      }

    outputImage->SetPixel(index, maxLabel);
  }
}

void PandoraBoxImageFilters::DiscreteGaussianFiltering(itkFloatImagePointer & inputImage, itkFloatImagePointer & outputImage, float variance)
{
  itkGaussianFilter::Pointer filter = itkGaussianFilter::New();
  filter->SetInput(inputImage);

  itkGaussianFilter::ArrayType var;
  var[0] = variance;
  var[1] = variance;
  var[2] = variance;
  filter->SetVariance(var);

  filter->SetUseImageSpacingOff();
  filter->Update();

  outputImage = filter->GetOutput();
}




} // namespace btk
