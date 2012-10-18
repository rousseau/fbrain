/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 15/10/2012
  Author(s): Julien Pontabry (pontabry@unistra.fr)

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

#include "btkResampleImagesToBiggestImageFilter.h"


// ITK includes
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"

// BTK includes
#include "btkMacro.h"
#include "btkImageHelper.h"


namespace btk
{

template< class TImage >
ResampleImagesToBiggestImageFilter< TImage >::ResampleImagesToBiggestImageFilter() : Superclass()
{
    Self::SetNumberOfRequiredInputs(0);
    Self::SetNumberOfRequiredOutputs(0);

    m_Interpolator = Self::ResampleImageFilter::LinearInterpolatorType::New();
}

//------------------------------------------------------------------------------------------------

template< class TImage >
ResampleImagesToBiggestImageFilter< TImage >::~ResampleImagesToBiggestImageFilter()
{
    // ----
}

//------------------------------------------------------------------------------------------------

template< class TImage >
void ResampleImagesToBiggestImageFilter< TImage >::GenerateData()
{
    // FIXME : the image with the maximal volume may not
    // containt all the other images. Instead, we should
    // search the bounding box of the whole images and
    // create a new image with this bounding box.
    // Note that the images are supposed to be closed enough.
    // Search the image with the maximum volume.
    unsigned int  j = 0;
    float maxVolume = 0;

    for(int i = 0; i < m_Images.size(); i++)
    {
        typename TImage::SpacingType spacing = m_Images[i]->GetSpacing();
        typename TImage::SizeType       size = m_Images[i]->GetLargestPossibleRegion().GetSize();
        float volume = size[0]*spacing[0]*size[1]*spacing[1]*size[2]*spacing[2];

        if(volume > maxVolume)
        {
            j = i;
            maxVolume = volume;
        }
    }

    // Define a resampler with the reference image with the maximal volume.
    typename Self::ResampleImageFilter::Pointer filter = ResampleImageFilter::New();

    filter->SetReferenceImage(m_Images[j]);
    filter->UseReferenceImageOn();
    filter->SetSize(m_Images[j]->GetLargestPossibleRegion().GetSize());
    filter->SetInterpolator(m_Interpolator);

    // Resample all input images and set them as output.
    for(int i = 0; i < m_Images.size(); i++)
    {
        if(i != j)
        {
            filter->SetInput(m_Images[i]);
            filter->Update();
            m_Images[i] = filter->GetOutput();
        }
    }
}

//------------------------------------------------------------------------------------------------

template< class TImage >
void ResampleImagesToBiggestImageFilter< TImage >::SetInputs(const std::vector< typename TImage::Pointer > &inputs)
{
    // Exceptions
    if(inputs.size() < 2)
    {
        btkException("ResampleImagesToBiggestImageFilter: There is not enough input images (it should be at least 2) !");
    }


    m_Images = inputs;
}

//------------------------------------------------------------------------------------------------

template< class TImage >
std::vector< typename TImage::Pointer > ResampleImagesToBiggestImageFilter< TImage >::GetOutputs()
{
    return m_Images;
}

} // namespace btk
