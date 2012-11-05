/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/04/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr)

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

#ifndef __BTK_AFFINEREGISTRATION_H__
#define __BTK_AFFINEREGISTRATION_H__

#include "btkRegistration.h"

#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkImageMaskSpatialObject.h"

#include "itkNumericTraits.h"
#include "btkUserMacro.h"


namespace btk
{

using namespace itk;

/** \class AffineRegistration
 * \brief Describe the class briefly here.
 *
 * Full class description
 * Full class description
 * Full class description

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

template <typename TImage>
class AffineRegistration :public btk::Registration <TImage >
{
public:
  /** Standard class typedefs. */
  typedef AffineRegistration  Self;
  typedef btk::Registration<TImage>   Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AffineRegistration, ImageRegistrationMethod);

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
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;

  typedef typename ImageType::RegionType                RegionType;



  /**  Type of the metric. */
  typedef MattesMutualInformationImageToImageMetric<
                                          ImageType,
                                          ImageType >   MetricType;
  typedef typename MetricType::Pointer                  MetricPointer;

  typedef typename MetricType::FixedImageMaskType   FixedImageMaskType;
  typedef typename FixedImageMaskType::Pointer     FixedImageMaskPointer;

  /**  Type of the Transform . */
  typedef AffineTransform<
                          double,
                          ImageType::ImageDimension >   TransformType;
  typedef typename TransformType::Pointer               TransformPointer;
  typedef typename TransformType::ParametersType        ParametersType;

  /**  Type of the Interpolator. */
  typedef LinearInterpolateImageFunction<
                                    ImageType,
                                    double>             InterpolatorType;
  typedef typename InterpolatorType::Pointer            InterpolatorPointer;

  /**  Type of the optimizer. */
  typedef RegularStepGradientDescentOptimizer           OptimizerType;
  typedef typename OptimizerType::Pointer               OptimizerPointer;
  typedef OptimizerType::ScalesType 					OptimizerScalesType;

  /** Method that initiates the registration. */

  /** Set/Get transform center. */
  itkSetMacro(TransformCenter, PointType);

  virtual PointType GetTransformCenter() const
    { return m_Transform -> GetCenter(); }

  //itkGetObjectMacro(Transform, TransformType);

  virtual TransformType* GetTransform()
  {
      return static_cast<TransformType*>(m_Transform);
  }


  itkSetObjectMacro(FixedImageMask, ImageMaskType);
  itkGetObjectMacro(FixedImageMask, ImageMaskType);

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
  virtual void InitializeWithTransform(){};


  /** Initialization is performed with the provided image masks. */
  virtual void InitializeWithMask(){};



protected:
  AffineRegistration();
  virtual ~AffineRegistration() {};
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */
  virtual void Initialize() throw (ExceptionObject);

private:
  AffineRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MetricPointer        m_Metric;
  OptimizerPointer     m_Optimizer;
  TransformPointer     m_Transform;
  InterpolatorPointer  m_Interpolator;

  unsigned int m_Iterations;
  bool m_EnableObserver;

  ImageMaskPointer m_FixedImageMask;

  MaskPointer          m_FixedMask;

  CommandIterationUpdate::Pointer  m_Observer;

  PointType                 m_TransformCenter;
};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkAffineRegistration.txx"
#endif

#endif
