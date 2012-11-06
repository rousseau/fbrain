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

#ifndef __BTK_APPLYTRANSFORMTOIMAGEFILTER_H__
#define __BTK_APPLYTRANSFORMTOIMAGEFILTER_H__

/* ITK */
#include "itkObject.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkSmartPointer.h"

/* BTK */

#include "btkMacro.h"

namespace btk
{
template <typename TImageIn, typename TImageOut>
class ApplyTransformToImageFilter : public itk::Object
{

public:
    typedef TImageIn itkImage;
    typedef typename itkImage::Pointer itkImagePointer;
    typedef TImageOut itkImageOut;
    typedef typename itkImageOut::Pointer itkImageOutPointer;
    typedef itk::Transform<double, 3,3> itkTransform;
    typedef itk::ResampleImageFilter<itkImage,itkImageOut> Resampler;

    typedef ApplyTransformToImageFilter Self;
    typedef itk::SmartPointer<Self>         Pointer;
    typedef itk::SmartPointer<const Self>     ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    virtual void Update();
    virtual void Initialize();

    btkSetMacro(InputImage,itkImagePointer);
    btkGetMacro(InputImage,itkImagePointer);

    btkGetMacro(OutputImage,itkImageOut*);


    btkSetMacro(Transform, itkTransform*);


protected:
    ApplyTransformToImageFilter();
    virtual ~ApplyTransformToImageFilter();

private:
    itkImagePointer m_InputImage;
    itkImageOut* m_OutputImage;
    itkTransform* m_Transform;
    typename Resampler::Pointer m_Resampler;




};

}

#include "btkApplyTransformToImageFilter.txx"

#endif // BTKWarpTRANSFORMTOIMAGEFILTER_H
