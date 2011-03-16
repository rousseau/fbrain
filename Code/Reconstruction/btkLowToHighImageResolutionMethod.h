/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/5/2010
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


/*=========================================================================

  Purpose: class for obtaining a high resolution image from different views
  by applying affine registration and averaging. First step of the method
  described in:

  Rousseau, F., Glenn, O.A., Iordanova, B., Rodriguez-Carranza, C.,
  Vigneron, D.B., Barkovich, J.A., Studholme, C.: Registration-based approach
  for reconstruction of high-resolution in utero fetal MR brain images. Acad
  Radiol 13(9) (Sep 2006) 1072–1081

=========================================================================*/


#ifndef __btkLowToHighImageResolutionMethod_h
#define __btkLowToHighImageResolutionMethod_h

#include "itkProcessObject.h"
#include "itkEuler3DTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "btkRigidRegistration.h"

#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"

#include "itkNumericTraits.h"
#include "btkUserMacro.h"

namespace btk
{

using namespace itk;

/** \class LowToHighImageResolutionMethod
 * \brief Describe the class briefly here.
 *
 * Full class description
 * Full class description
 * Full class description

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

// TODO Change the template to provide two types of images (Low and High resolution)

template <typename TImage>
class LowToHighImageResolutionMethod : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef LowToHighImageResolutionMethod  Self;
  typedef ProcessObject                                Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LowToHighImageResolutionMethod, ProcessObject);

  /**  Type of the Fixed image. */
  typedef          TImage                               ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;
  typedef          std::vector<ImagePointer>            ImageArrayPointer;

  typedef typename ImageType::RegionType                ImageRegionType;
  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::PointType                 PointType;
  typedef typename ImageType::SpacingType               SpacingType;
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;

  typedef typename ImageType::RegionType               	RegionType;
  typedef  std::vector<RegionType>                      RegionArray;

  typedef itk::Image< unsigned char, 3 >                ImageMaskType;
  typedef std::vector<ImageMaskType::Pointer>           ImageMaskArray;

  typedef ImageFileWriter< ImageType >  ImageWriterType;

  typedef ImageRegionIteratorWithIndex< ImageType >  IteratorType;

  /**  Type of the Transform . */
  typedef Euler3DTransform< double >             				TransformType;
  typedef typename TransformType::Pointer               TransformPointer;
  typedef  std::vector<TransformPointer>                TransformPointerArray;

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

  /** Write the final transforms to disk. */
  void WriteTransforms( const char *folder );

  /** Write a specific transform to disk. */
  void WriteTransforms( unsigned int i, const char *filename );

  /** Write the resampled images to disk. */
  void WriteResampledImages( const char *folder );

  /** Write a specific resampled image to disk. */
  void WriteResampledImages( unsigned int i, const char *filename );

  /** Get the resampled image array. */
  UserGetObjectMacro( ResampledImageArray, ImageType );

  /** Method to stop the registration. */
//  void StopRegistration();

  /** Get the high resolution Image. */
  itkGetObjectMacro( HighResolutionImage, ImageType );

  /** Set/Get the FixedImageRegion. */
  itkSetMacro( FixedImageRegion, ImageRegionType );
  itkGetMacro( FixedImageRegion, ImageRegionType );

  /** Set/Get the Initialization. */
  itkSetMacro( InitializeWithMask, bool );
  itkGetMacro( InitializeWithMask, bool );

  /** Set/Get the target image. */
  itkSetMacro( TargetImage, unsigned int );
  itkGetMacro( TargetImage, unsigned int );

  /** Set/Get the margin. */
  itkSetMacro( Margin, double );
  itkGetMacro( Margin, double );

  /** Set/Get the Transfrom. */
  UserSetObjectMacro( TransformArray, TransformType );
  UserGetObjectMacro( TransformArray, TransformType );

  /** Set/Get the image array. */
  UserSetObjectMacro( ImageArray, ImageType );
  UserGetObjectMacro( ImageArray, ImageType );

  /** Set/Get the image array. */
  UserSetObjectMacro( ImageMaskArray, ImageMaskType );
  UserGetObjectMacro( ImageMaskArray, ImageMaskType );

  /** Set/Get the image array. */
  UserSetMacro( RegionArray, RegionType );
//  UserGetMacro( RegionArray, RegionType );

  /** Set the number of images */
  void SetNumberOfImages(int N);

  /** write output images and extract their slices after each resolution */


protected:
  LowToHighImageResolutionMethod();
  virtual ~LowToHighImageResolutionMethod() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */
  void Initialize() throw (ExceptionObject);



private:
  LowToHighImageResolutionMethod(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  TransformPointerArray            m_TransformArray;
  ImageArrayPointer 			         m_ImageArray;
  RegionArray                      m_RegionArray;
  ImageMaskArray                   m_ImageMaskArray;
//  RegistrationPointer							 m_Registration;
//  ParametersArrayType              m_InitialRigidParameters;

  SpacingType											 m_ResampleSpacing;
  SizeType											   m_ResampleSize;
  ImageArrayPointer						 		 m_ResampledImageArray;
  ResamplePointer                  m_Resample;
  std::vector<bool>                m_ResamplingStatus;
  double                           m_Margin;

  InterpolatorPointer              m_Interpolator;

  ImagePointer										 m_HighResolutionImage;

  ParametersType                   m_InitialTransformParameters;
  ParametersType                   m_LastTransformParameters;
  ImageRegionType                  m_FixedImageRegion;

  unsigned int                     m_NumberOfImages;
  unsigned int                     m_TargetImage;

  bool                             m_InitializeWithMask;

};


} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkLowToHighImageResolutionMethod.txx"
#endif

#endif
