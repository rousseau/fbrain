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

#ifndef BTK_MASKIMAGEFILTER_H
#define BTK_MASKIMAGEFILTER_H

/* ITK */
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"

/* BTK */
#include "btkMacro.h"
#include "btkImageHelper.h"
namespace btk
{

/**
 * Class for apply a mask on an image, templated on input image, mask image and output image.
 * A Threshold can be applied if mask is a probability map, default value is 0.5 (all pixel in mask strictly greater
 * than threshold are considerated).
 * @author Marc Schweitzer
 * \ingroup ImageFilters
 */
template < typename TImageIn, typename TMaskIn, typename TImageOut = TImageIn >
class MaskImageFilter
{

        //TODO: Do an ITK::Filter for using it into a pipeline.
    public:
        /** Typedefs */
        typedef TImageIn ImageTypeIn;
        typedef TImageOut ImageTypeOut;
        typedef TMaskIn   MaskType;
        typedef itk::ImageRegionConstIteratorWithIndex<ImageTypeIn> ConstIteratorImageIn;
        typedef itk::ImageRegionConstIteratorWithIndex<MaskType> ConstIteratorMask;
        typedef itk::ImageRegionIteratorWithIndex<ImageTypeOut> IteratorImageOut;

        /** Constructor */
        MaskImageFilter();
        /** Destructor */
        virtual ~MaskImageFilter();

       /** Set/Get Input image */
        btkSetMacro(Input,typename ImageTypeIn::Pointer );
        btkGetMacro(Input,typename ImageTypeIn::Pointer );

        /** Set/Get Mask image */
        btkSetMacro(Mask,typename MaskType::Pointer );
        btkGetMacro(Mask,typename MaskType::Pointer );

        /** Get Output image (result of the filter) */
        btkGetMacro(Output,typename ImageTypeOut::Pointer );

        /** Set/Get Threshold (if mask is a probability map) default value is 0.5 */
        btkSetMacro(Threshold, float);
        btkGetMacro(Threshold, float);

        /** Update Method */
        void Update() throw(itk::ExceptionObject &);

    protected :

    private :
        /** SmartPointer on input image */
        typename ImageTypeIn::Pointer m_Input;
        /** SmartPointer on mask image */
        typename MaskType::Pointer m_Mask;
        /** SmartPointer on ouput image */
        typename ImageTypeOut::Pointer m_Output;
        /** Threshold for mask */
        float m_Threshold;

};

}
#ifndef ITK_MANUAL_INSTANTIATION
#include "btkMaskImageFilter.txx"
#endif

#endif // BTKMASKIMAGEFILTER_H
