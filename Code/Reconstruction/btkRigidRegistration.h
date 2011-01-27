/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 12/11/2010
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

#ifndef __btkRigidRegistration_h
#define __btkRigidRegistration_h

#include "itkImageRegistrationMethod.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "itkCenteredTransformInitializer.h"

#include "itkImageMaskSpatialObject.h"

#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"

#include "itkNumericTraits.h"
#include "btkUserMacro.h"

#include "itkCommand.h"

#ifndef __CommandIterationUpdate
#define __CommandIterationUpdate

class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer =
        dynamic_cast< OptimizerPointer >( object );
      if( ! itk::IterationEvent().CheckEvent( &event ) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer -> GetCurrentStepLength() << "   ";
      std::cout << optimizer->GetValue() << "   ";

      std::cout << optimizer->GetCurrentPosition() << std::endl;
    }
};

#endif

namespace btk
{
/** \class btkRigidRegistration
 * \brief Initialized registration method.
 *
 * This class aims to encapsulate registration code common to many
 * application/classes of the library. The used components provided good
 * results in fetal brain imaging, and therefore can remain fixed.

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

using namespace itk;

template <typename TImage>
class RigidRegistration :
public ImageRegistrationMethod <TImage,TImage>
{
public:
  /** Standard class typedefs. */
  typedef RigidRegistration  Self;
  typedef ImageRegistrationMethod<TImage,TImage>       Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RigidRegistration, ImageRegistrationMethod);

  /**  Type of the Fixed image. */
  typedef          TImage                               ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;

  typedef Image< unsigned char, 3 >                     ImageMaskType;
  typedef typename ImageMaskType::Pointer               ImageMaskPointer;

  typedef ImageMaskSpatialObject< 3 >                   MaskType;
  typedef typename MaskType::Pointer                    MaskPointer;

  typedef typename ImageType::PointType                 PointType;
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;

  /**  Type of the metric. */
  // typedef MattesMutualInformationImageToImageMetric<
  //typedef MeanSquaresImageToImageMetric<
  // FIXME Change metric to MI after testing NC
  typedef NormalizedCorrelationImageToImageMetric<
                                          ImageType,
                                          ImageType >   MetricType;



  typedef typename MetricType::Pointer                  MetricPointer;

  /* TODO: Check if better results can be obtained with other transformation types
     like VersorRigid3DTransform. The ITK user's guide says this transformation
     could have problems for large rotations */

  /**  Type of the Transform . */
  typedef Euler3DTransform< double >                    TransformType;
  typedef typename TransformType::Pointer               TransformPointer;
  typedef  TransformType::ParametersType                ParametersType;

  typedef CenteredTransformInitializer<
                                    TransformType,
                                    ImageMaskType,
                                    ImageMaskType >  TransformInitializerType;


  /**  Type of the Interpolator. */
  typedef LinearInterpolateImageFunction<
                                    ImageType,
                                    double>             InterpolatorType;
  typedef typename InterpolatorType::Pointer            InterpolatorPointer;

  /**  Type of the optimizer. */
  typedef RegularStepGradientDescentOptimizer           OptimizerType;
  typedef typename OptimizerType::Pointer               OptimizerPointer;
  typedef OptimizerType::ScalesType 										OptimizerScalesType;


  /** Set/Get transform. */
  itkGetObjectMacro(Transform, TransformType);

  /** Set/Get transform center. */
  itkSetMacro(TransformCenter, PointType);
  itkGetMacro(TransformCenter, PointType);

    /** Set/Get iterations. */
  itkSetMacro(Iterations, unsigned int);
  itkGetMacro(Iterations, unsigned int);

  /** Set/Get observer status. */
  itkSetMacro(EnableObserver, bool);
  itkGetMacro(EnableObserver, bool);

  /** Set/Get image mask for fixed image. */
  itkSetObjectMacro(FixedImageMask, ImageMaskType);
  itkGetObjectMacro(FixedImageMask, ImageMaskType);

  /** Set/Get image mask for moving image. */
  itkSetObjectMacro(MovingImageMask, ImageMaskType);
  itkGetObjectMacro(MovingImageMask, ImageMaskType);

  void InitializeWithTransform()
  {
    m_InitializeWithTransform = true;
    m_InitializeWithMask = false;
    m_InitializeWithRegion = false;
  }

  void InitializeWithMask()
  {
    m_InitializeWithTransform = false;
    m_InitializeWithMask = true;
    m_InitializeWithRegion = false;
  }

  void InitializeWithRegion()
  {
    m_InitializeWithTransform = false;
    m_InitializeWithMask = false;
    m_InitializeWithRegion = true;
  }



protected:
  RigidRegistration();
  virtual ~RigidRegistration() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */
  void Initialize() throw (ExceptionObject);

private:
  RigidRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MetricPointer        m_Metric;
  OptimizerPointer     m_Optimizer;
  TransformPointer     m_Transform;
  InterpolatorPointer  m_Interpolator;

  ImageMaskPointer     m_FixedImageMask;
  ImageMaskPointer     m_MovingImageMask;

  MaskPointer          m_FixedMask;

  PointType            m_TransformCenter;

  bool m_InitializeWithMask;
  bool m_InitializeWithRegion;
  bool m_InitializeWithTransform;


  unsigned int m_Iterations;
  bool         m_EnableObserver;

  CommandIterationUpdate::Pointer  m_Observer;
};


} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkRigidRegistration.txx"
#endif

#endif
