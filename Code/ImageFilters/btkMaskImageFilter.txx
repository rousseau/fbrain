/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 10/01/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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
#ifndef BTK_MASKIMAGEFILTER_TXX
#define BTK_MASKIMAGEFILTER_TXX

#include "btkMaskImageFilter.h"

namespace btk
{
//-------------------------------------------------------------------------------------------------
template < typename TImageIn, typename TMaskIn, typename TImageOut>
MaskImageFilter<  TImageIn,  TMaskIn,  TImageOut>::MaskImageFilter()
{
    m_Threshold = 0.5;
    m_Output = NULL;
    m_Input =  NULL;
    m_Mask = NULL;
}
//-------------------------------------------------------------------------------------------------
template < typename TImageIn, typename TMaskIn, typename TImageOut>
MaskImageFilter<  TImageIn,  TMaskIn,  TImageOut>::~MaskImageFilter()
{

}
//-------------------------------------------------------------------------------------------------
template < typename TImageIn, typename TMaskIn, typename TImageOut>
void
MaskImageFilter<  TImageIn,  TMaskIn,  TImageOut>::Update() throw(itk::ExceptionObject &)
{
    if(!m_Input || !m_Mask)
    {
        btkException("Missing input image, or mask.");
    }


    m_Output = ImageTypeOut::New();
    m_Output = ImageHelper<ImageTypeIn, ImageTypeOut>::CreateNewImageFromPhysicalSpaceOf(m_Input.GetPointer());

    // iterate over mask

    ConstIteratorMask it_mask(m_Mask, m_Mask->GetLargestPossibleRegion());
    ConstIteratorImageIn it_input(m_Input, m_Input->GetLargestPossibleRegion());
    IteratorImageOut it_output(m_Output, m_Output->GetLargestPossibleRegion());

    for(it_mask.GoToBegin(), it_input.GoToBegin(), it_output.GoToBegin();
        !it_mask.IsAtEnd(), !it_input.IsAtEnd(), !it_output.IsAtEnd() ;
        ++it_mask, ++it_input, ++it_output)
    {
        if(it_mask.Get() > m_Threshold)
        {
            it_output.Set(it_input.Get());
        }
    }


}
//-------------------------------------------------------------------------------------------------
}


#endif
