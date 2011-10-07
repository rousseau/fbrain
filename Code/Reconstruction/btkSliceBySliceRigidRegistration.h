/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 24/01/2011
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

#ifndef __btkSliceBySliceRigidRegistration_h
#define __btkSliceBySliceRigidRegistration_h

#include "itkProcessObject.h"
#include "itkEuler3DTransform.h"
#include "btkSliceBySliceTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "btkRigidRegistration.h"

#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"

#include "itkImageMaskSpatialObject.h"

#include "itkNumericTraits.h"
#include "btkUserMacro.h"

namespace btk
{

using namespace itk;

/** \class SliceBySliceRigidRegistration
 * \brief This class registers two images by using a slice by slice rigid transform.
 *
 * Full class description
 * Full class description
 * Full class description

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

// TODO Change the template to provide two types of images (Low and High resolution)

template <typename TImage>
class SliceBySliceRigidRegistration : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SliceBySliceRigidRegistration  Self;
  typedef ProcessObject                                Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SliceBySliceRigidRegistration, ProcessObject);

  /**  Type of the Fixed image. */
  typedef          TImage                               ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;
  typedef          std::vector<ImagePointer>            ImageArrayPointer;

  itkStaticConstMacro(ImageDimension, unsigned int,ImageType::ImageDimension);

  typedef typename ImageType::RegionType                ImageRegionType;
  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::PointType                 PointType;
  typedef typename ImageType::SpacingType               SpacingType;
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;

  typedef typename ImageType::RegionType               	RegionType;
  typedef  std::vector<RegionType>                      RegionArray;

  typedef itk::Image< unsigned char, 3 >                ImageMaskType;
  typedef typename ImageMaskType::Pointer               ImageMaskPointer;
  typedef std::vector<ImageMaskType::Pointer>           ImageMaskArray;

  typedef itk::ImageMaskSpatialObject< ImageDimension > MaskType;

  typedef ImageFileWriter< ImageType >  ImageWriterType;

  typedef ImageRegionIteratorWithIndex< ImageType >  IteratorType;

  /**  Type of the Transform . */
  typedef Euler3DTransform< double >             				TransformType;
  typedef typename TransformType::Pointer               TransformPointer;
  typedef  std::vector<TransformPointer>                TransformPointerArray;

  typedef SliceBySliceTransform< double, ImageDimension >  SliceBySliceTransformType;
  typedef typename SliceBySliceTransformType::Pointer   SliceBySliceTransformPointer;

  /**  Type of the transform reader. */
  typedef itk::TransformFileReader                      TransformReaderType;
  typedef TransformReaderType::TransformListType*       TransformListType;

  /**  Type of the Interpolator. */
  typedef LinearInterpolateImageFunction<
                                    ImageType,
                                    double>             InterpolatorType;
  typedef typename InterpolatorType::Pointer            InterpolatorPointer;

   /**  Type of the registration method. */
  typedef RigidRegistration<ImageType>  RegistrationType;
  typedef typename RegistrationType::Pointer     RegistrationPointer;

  typedef ResampleImageFilter< ImageType, ImageType >    ResampleType;
  typedef typename ResampleType::Pointer 								 ResamplePointer;


  typedef  TransformType::ParametersType    ParametersType;
  typedef  std::vector<ParametersType>      ParametersArrayType;

  typedef TransformFileWriter TransformWriterType;

  /** Method that initiates the registration. */
  void StartRegistration();

  /** Set/Get the fixed image. */
  itkSetObjectMacro( FixedImage, ImageType );
  itkGetObjectMacro( FixedImage, ImageType );

  /** Set/Get the moving image. */
  itkSetObjectMacro( MovingImage, ImageType );
  itkGetObjectMacro( MovingImage, ImageType );

  /** Set/Get the image mask. */
  itkSetObjectMacro( ImageMask, ImageMaskType );
  itkGetObjectMacro( ImageMask, ImageMaskType );

  /** Set/Get the transform. */
  itkSetObjectMacro( Transform, SliceBySliceTransformType );
  itkGetObjectMacro( Transform, SliceBySliceTransformType );

  /** Set/Get the number of iterations. */
  itkSetMacro( Iterations, unsigned int );
  itkGetMacro( Iterations, unsigned int );


protected:
  SliceBySliceRigidRegistration();
  virtual ~SliceBySliceRigidRegistration() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */
  void Initialize() throw (ExceptionObject);



private:
  SliceBySliceRigidRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ImagePointer                     m_FixedImage;
  ImagePointer                     m_MovingImage;
  ImageMaskPointer                 m_ImageMask;

  SliceBySliceTransformPointer     m_Transform;

  RegionType                       m_ROI;
  RegistrationPointer							 m_Registration;

  InterpolatorPointer              m_Interpolator;

  ParametersType                   m_InitialTransformParameters;
  ParametersType                   m_LastTransformParameters;
  ImageRegionType                  m_FixedImageRegion;

  unsigned int 										 m_Iterations;

};


} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSliceBySliceRigidRegistration.txx"
#endif

#endif
