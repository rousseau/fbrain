/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 09/05/2012
  Author(s): Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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

#ifndef __BTK_REGISTRATION_H__
#define __BTK_REGISTRATION_H__

#include "itkImageRegistrationMethod.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkTransform.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

#include "itkNumericTraits.h"
#include "btkMacro.h"
#include "btkCommandIterationUpdate.h"



namespace btk
{

using namespace itk;

/** \class Registration
 * \brief Describe the class briefly here.
 *
 * Full class description
 * Full class description
 * Full class description

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

template <typename TImage>
class Registration :public ImageRegistrationMethod <TImage,TImage>
{
    public:
        /** Standard class typedefs. */
        typedef Registration  Self;
        typedef ImageRegistrationMethod<TImage,TImage>       Superclass;
        typedef SmartPointer<Self>                           Pointer;
        typedef SmartPointer<const Self>                     ConstPointer;

        /** Method for creation through the object factory. */
        //itkNewMacro(Self); // pure virtual

        /** Run-time type information (and related methods). */
        itkTypeMacro(Registration, ImageRegistrationMethod);

        /**  Type of the Fixed image. */
        typedef          TImage                               ImageType;
        typedef typename ImageType::Pointer                   ImagePointer;

        /**  Type of the image mask. */
        typedef Image< unsigned char, 3 >                     ImageMaskType;
        typedef typename ImageMaskType::Pointer               ImageMaskPointer;

        /**  Type of the image mask spatial object. */
        typedef ImageMaskSpatialObject< 3 >                   MaskType;
        typedef typename MaskType::Pointer                    MaskPointer;

        typedef typename ImageType::PointType                 PointType;


        /**  Type of the metric. */
        /**  Type of the metric.
   * --------------------- /!\ Warning /!\ ------------------------------
   *
   * - itk::MattesMutualInformationImageToImageMetric is NOT THREAD SAFE....
   *
   * ---------------------------------------------------------------------
   *
   * typedef MattesMutualInformationImageToImageMetric<
   */
        //typedef MattesMutualInformationImageToImageMetric<
        //                                      ImageType,
        //                                    ImageType >   MetricType;


        typedef NormalizedCorrelationImageToImageMetric<
        ImageType,
        ImageType >   MetricType;


        typedef typename MetricType::Pointer                  MetricPointer;



        /**  Type of the Transform . */
        typedef itk::Transform< double, ImageType::ImageDimension,
                                        ImageType::ImageDimension > TransformType;
        typedef typename TransformType::Pointer               TransformPointer;
        typedef typename TransformType::ParametersType        ParametersType;

        /**  Type of the optimizer. */
        typedef RegularStepGradientDescentOptimizer           OptimizerType;

        /**  Type of the interpolator. */
        typedef LinearInterpolateImageFunction<
                                          ImageType,
                                          double>             InterpolatorType;
        typedef typename InterpolatorType::Pointer            InterpolatorPointer;


        //virtual TransformType* GetTransform(){};

        //  virtual void SetFixedImageMask(ImageMaskType *imMask)
        //  {
        //      m_FixedImageMask = imMask;
        //  };
        //  virtual ImageMaskType* GetFixedImageMask()
        //  {
        //      return m_FixedImageMask.GetPointer();
        //  };

        //  virtual void SetIterations(const unsigned int it){};
        //  virtual unsigned int GetIterations(){};

        //  virtual void SetEnableObserver(const bool arg){};
        //  virtual bool GetEnableObserver(){};

        //  virtual OptimizerType* GetOptimizer(){};

        //  virtual void SetMovingImageMask(ImageMaskType* imMask){};
        //  virtual ImageMaskType* GetMovingImageMask(){};

        /** Method that initiates the registration. */

        /** Set/Get transform center. */
        itkSetMacro(TransformCenter, PointType);

        virtual PointType GetTransformCenter() const
        {
            return m_TransformCenter;
        }


          virtual TransformType* GetTransform()
          {
            return m_Transform.GetPointer();
          }

        itkSetObjectMacro(FixedImageMask, ImageMaskType);
        itkGetObjectMacro(FixedImageMask, ImageMaskType);

        /** Set/Get image mask for moving image. */
        itkSetObjectMacro(MovingImageMask, ImageMaskType);
        itkGetObjectMacro(MovingImageMask, ImageMaskType);

        /** Set/Get iterations. */
        itkSetMacro(Iterations, unsigned int);
        itkGetMacro(Iterations, unsigned int);

        /** Set/Get observer status. */
        itkSetMacro(EnableObserver, bool);
        itkGetMacro(EnableObserver, bool);

        itkSetObjectMacro(FixedMask,MaskType);
        itkGetObjectMacro(FixedMask,MaskType);

        itkGetObjectMacro(Optimizer, OptimizerType);

        /** Initialization is performed with the provided transform. */
        virtual void InitializeWithTransform() = 0;


        /** Initialization is performed with the provided image masks. */
        virtual void InitializeWithMask() = 0;


    protected:
        Registration(){};
        virtual ~Registration() {};

        virtual void PrintSelf(std::ostream& os, Indent indent) const
        {
            Superclass::PrintSelf( os, indent );
        }


        /** Initialize by setting the interconnects between the components.
   */
        virtual void Initialize() throw (ExceptionObject)
        {
            Superclass::Initialize();
        }

        typename MetricType::Pointer        m_Metric;
        typename OptimizerType::Pointer     m_Optimizer;
        typename TransformType::Pointer         m_Transform;
        typename InterpolatorType::Pointer  m_Interpolator;

        unsigned int m_Iterations;
        bool m_EnableObserver;

        ImageMaskPointer     m_FixedImageMask;
        ImageMaskPointer     m_MovingImageMask;

        MaskPointer          m_FixedMask;

        CommandIterationUpdate::Pointer  m_Observer;

        PointType                 m_TransformCenter;


    private:
        Registration(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented




};


} // end namespace itk


#endif
