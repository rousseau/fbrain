/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 02/04/2012
  Author(s): Marc Schweitzer (marc.schweitzer@unistra.fr)

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

#ifndef __BTK_APPLYTRANSFORMTOIMAGEFILTER_TXX__
#define __BTK_APPLYTRANSFORMTOIMAGEFILTER_TXX__

#include "btkApplyTransformToImageFilter.h"


namespace btk
{
template<typename TImageIn, typename TImageOut>
ApplyTransformToImageFilter<TImageIn,TImageOut>::ApplyTransformToImageFilter()
{
    m_OutputImage = itkImageOut::New();
}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
ApplyTransformToImageFilter<TImageIn,TImageOut>::~ApplyTransformToImageFilter()
{
    m_Transform = NULL;
    m_OutputImage = NULL;
    m_InputImage = NULL;

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Initialize()
{
    if(!m_Transform)
    {
        btkException("Missing input transform !");
    }
    if(!m_InputImage)
    {
        btkException("Missing input image !");
    }
    if(!m_ReferenceImage)
    {
        btkException("Missing reference image !");
    }

}
//-------------------------------------------------------------------------------------------
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Update()
{
    try
    {
        this->Resample();
    }
    catch(itk::ExceptionObject & excp)
    {
        std::cerr << "Error when apply transformation" << std::endl;
        std::cerr << excp << std::endl;
        std::cout << "[FAILED]" << std::endl;
        throw excp;
    }



}
template<typename TImageIn, typename TImageOut>
void ApplyTransformToImageFilter<TImageIn,TImageOut>::Resample() throw(itk::ExceptionObject &)
{
    m_OutputImage = btk::ImageHelper<itkImageOut>::DeepCopy(m_ReferenceImage.GetPointer());

    IteratorIn it(m_InputImage, m_InputImage->GetLargestPossibleRegion());
    IteratorOut itO(m_OutputImage, m_OutputImage->GetLargestPossibleRegion());

    typename Interpolator::Pointer interpolator = Interpolator::New();
    interpolator->SetInputImage(m_InputImage);

    for(it.GoToBegin(), itO.GoToBegin(); !it.IsAtEnd() && !itO.IsAtEnd(); ++it, ++itO)
    {
        typename itkImage::IndexType InputIndex = it.GetIndex();
        typename itkImage::IndexType OutputIndex = itO.GetIndex();
        typename itkImage::PointType point, Tpoint;
        typename Interpolator::ContinuousIndexType Tindex;

        m_OutputImage->TransformIndexToPhysicalPoint(OutputIndex, point);

        Tpoint = m_Transform->TransformPoint(point);

         m_InputImage->TransformPhysicalPointToContinuousIndex(Tpoint, Tindex);

        if(interpolator->IsInsideBuffer(Tindex))
        {
            typename itkImage::PixelType value = interpolator->EvaluateAtContinuousIndex(Tindex);
            itO.Set(value);
        }
        else
        {
            itO.Set(0);
        }
    }


}
}

#endif
