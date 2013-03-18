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
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkTransform.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkSmartPointer.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"


/* BTK */

#include "btkMacro.h"
#include "btkImageHelper.h"

namespace btk
{
/**
 * Class for apply a transformation to an image, templated on input and output image
 * @author Marc Schweitzer
 * \ingroup ImageFilters
 */
template <typename TImageIn, typename TImageOut>
class ApplyTransformToImageFilter : public itk::Object
{

    public:
        /** Typedefs */
        typedef TImageIn itkImage;
        typedef typename itkImage::Pointer itkImagePointer;
        typedef TImageOut itkImageOut;
        typedef typename itkImageOut::Pointer itkImageOutPointer;
        typedef itk::Transform<double, 3,3> itkTransform;
        typedef itk::LinearInterpolateImageFunction<itkImage> Interpolator;
        typedef itk::ImageRegionIteratorWithIndex<itkImage> IteratorIn;
        typedef itk::ImageRegionIteratorWithIndex<itkImageOut> IteratorOut;
        typedef itk::ContinuousIndex<itkImage> ContinuousIndexIn;
        typedef itk::ContinuousIndex<itkImageOut> ContinuousIndexOut;

        typedef ApplyTransformToImageFilter Self;
        typedef itk::SmartPointer<Self>         Pointer;
        typedef itk::SmartPointer<const Self>     ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        /** Update Method */
        virtual void Update();
        /** Initialize Method to call before Update() */
        virtual void Initialize();

        /** Set/Get Input Image */
        btkSetMacro(InputImage,itkImagePointer);
        btkGetMacro(InputImage,itkImagePointer);

        /** Set/Get Reference Image */
        btkSetMacro(ReferenceImage, itkImagePointer);
        btkGetMacro(ReferenceImage, itkImagePointer);

        /** Get Output Image (result of the filter) */
        btkGetMacro(OutputImage,itkImageOutPointer);

        /** Set Transform to Apply */
        btkSetMacro(Transform, itkTransform*);


    protected:
        /** Constructor */
        ApplyTransformToImageFilter();
        /** Destructor */
        virtual ~ApplyTransformToImageFilter();
        /** Resample method, called by Update(),  it can throw an exception */
        virtual void Resample() throw(itk::ExceptionObject &);

    private:
        /** SmartPointer on Input Image */
        itkImagePointer m_InputImage;
        /** SmartPointer on Reference Image */
        itkImagePointer m_ReferenceImage;
        /** SmartPointer on Output Image */
        itkImageOutPointer m_OutputImage;
        /** Pointer to an itk::Transform */
        itkTransform* m_Transform;




};

}

#include "btkApplyTransformToImageFilter.txx"

#endif // BTKWarpTRANSFORMTOIMAGEFILTER_H
