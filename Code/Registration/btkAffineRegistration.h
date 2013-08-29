/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/04/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr)
             Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkImageMaskSpatialObject.h"
#include "itkNumericTraits.h"
#include "btkUserMacro.h"


namespace btk
{

/** @class AffineRegistration
 * @brief This class is usefull to do affine registration
 * @author Marc Schweitzer
 * @ingroup Registration
 */

template <typename TImage>
class AffineRegistration :public btk::Registration < TImage >
{
public:
  /** Standard class typedefs. */
  typedef AffineRegistration  Self;
  typedef btk::Registration<TImage>   Superclass;
  typedef itk::SmartPointer<Self>                           Pointer;
  typedef itk::SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AffineRegistration, ImageRegistrationMethod);

  /**  Type of the Fixed image. */
  typedef          TImage                               ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;

  /**  Type of the image mask. */
  typedef itk::Image< unsigned char, 3 >                     ImageMaskType;
  typedef typename ImageMaskType::Pointer               ImageMaskPointer;

  /**  Type of the image mask spatial object. */
  typedef itk::ImageMaskSpatialObject< 3 >                   MaskType;
  typedef typename MaskType::Pointer                    MaskPointer;

  typedef typename ImageType::PointType                 PointType;
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;

  typedef typename ImageType::RegionType                RegionType;


  /**  Type of the metric.
   * --------------------- /!\ Warning /!\ ------------------------------
   *
   * - itk::MattesMutualInformationImageToImageMetric is NOT THREAD SAFE....
   *
   * ---------------------------------------------------------------------
   *
   * typedef MattesMutualInformationImageToImageMetric<
   */

  typedef itk::NormalizedCorrelationImageToImageMetric<
                                          ImageType,
                                          ImageType >   MetricType;

  typedef typename MetricType::Pointer                  MetricPointer;


  typedef typename MetricType::FixedImageMaskType   FixedImageMaskType;
  typedef typename FixedImageMaskType::Pointer     FixedImageMaskPointer;

  /**  Type of the Transform . */
  typedef itk::AffineTransform<
                          double,
                          ImageType::ImageDimension >   AffineTransformType;
  typedef typename AffineTransformType::Pointer               TransformPointer;
  typedef typename AffineTransformType::ParametersType        ParametersType;

  /**  Type of the Interpolator. */
  typedef itk::LinearInterpolateImageFunction<
                                    ImageType,
                                    double>             InterpolatorType;
  typedef typename InterpolatorType::Pointer            InterpolatorPointer;

  /**  Type of the optimizer. */
  typedef itk::RegularStepGradientDescentOptimizer           OptimizerType;
  typedef typename OptimizerType::Pointer               OptimizerPointer;
  typedef OptimizerType::ScalesType 					OptimizerScalesType;

  virtual AffineTransformType* GetTransform()
  {
      return dynamic_cast<AffineTransformType*>(Superclass::m_Transform.GetPointer());
  }



  /** Initialization is performed with the provided transform. */
  virtual void InitializeWithTransform(){};

  /** Initialization is performed with the provided image masks. */
  virtual void InitializeWithMask(){};



protected:
  AffineRegistration();
  virtual ~AffineRegistration() {};
  virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */
  virtual void Initialize() throw (itk::ExceptionObject);

private:
  AffineRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};


} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkAffineRegistration.txx"
#endif

#endif
