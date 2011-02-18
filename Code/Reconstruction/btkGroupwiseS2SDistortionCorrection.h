/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/03/2010
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

#ifndef __btkGroupwiseS2SDistortionCorrection_h
#define __btkGroupwiseS2SDistortionCorrection_h

#include "itkProcessObject.h"
#include "itkImageToImageMetric.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"

#include "btkAffineRegistration.h"
#include "itkTransformFileWriter.h"

#include "itkImageMaskSpatialObject.h"
#include "itkExtractImageFilter.h"

#include "itkDiscreteGaussianImageFilter.h"
#include "btkSliceBySliceRegistration.h"
#include "btkSliceToSliceInterpolateImageFunction.h"

#include "itkImageFileWriter.h"

#include "btkUserMacro.h"


namespace btk
{

using namespace itk;

/** \class GroupwiseS2SDistortionCorrection
 * \brief Describe the class briefly here.
 *
 * Full class description
 * Full class description
 * Full class description

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

template <typename TSequence >
class GroupwiseS2SDistortionCorrection : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef GroupwiseS2SDistortionCorrection  Self;
  typedef ProcessObject                                Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GroupwiseS2SDistortionCorrection, ProcessObject);

  typedef          TSequence                            SequenceType;
  typedef typename SequenceType::Pointer                SequencePointer;

  typedef typename SequenceType::PixelType              SequencePixelType;
  typedef typename SequenceType::RegionType             SequenceRegionType;
  typedef typename SequenceType::SizeType             	SequenceSizeType;
  typedef typename SequenceType::IndexType             	SequenceIndexType;

  typedef ImageRegionIteratorWithIndex< SequenceType > SequenceIteratorType;


  /**  Type of the Fixed image. */
  typedef          Image< SequencePixelType, 3 >        ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;
  typedef          std::vector<ImagePointer>            ImageArrayPointer;

  typedef          Image< unsigned char, 3 >       			ImageMaskType;
  typedef typename ImageMaskType::Pointer               ImageMaskPointer;

  typedef ImageMaskSpatialObject< 3 >   								MaskType;
  typedef typename MaskType::Pointer               			MaskPointer;

  typedef typename ImageType::PixelType                 PixelType;
  typedef typename ImageType::PointType                 PointType;
  typedef typename ImageType::SpacingType               SpacingType;
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;
  typedef typename ImageType::RegionType               	RegionType;

  typedef ImageRegionIteratorWithIndex< ImageType >  IteratorType;

  /**  Type of the Transform . */
  typedef AffineTransform< double, 3>           				TransformType;
  typedef typename TransformType::Pointer               TransformPointer;
  typedef  std::vector<TransformPointer>                S2STransformType;
  typedef  std::vector<S2STransformType>                S2STransformArray;

  typedef TransformFileWriter                           TransformWriterType;

  typedef typename TransformType::MatrixType            MatrixType;
  typedef typename TransformType::OffsetType            OffsetType;

  /**  Type of the Interpolator. */
  typedef SliceToSliceInterpolateImageFunction<
                                      ImageType,
                                      double >          RBFInterpolatorType;
  typedef typename RBFInterpolatorType::Pointer         RBFInterpolatorPointer;

    typedef LinearInterpolateImageFunction<
                                    ImageType,
                                    double>             InterpolatorType;
  typedef typename InterpolatorType::Pointer            InterpolatorPointer;

   /**  Type of the registration method. */
  typedef SliceBySliceRegistration< ImageType >    RegistrationType;
  typedef typename RegistrationType::Pointer            RegistrationPointer;

  /** Type of affine registration method */
  typedef AffineRegistration<ImageType>                 AffineRegistrationType;
  typedef typename AffineRegistrationType::Pointer      AffineRegistrationPointer;

  typedef ResampleImageFilter< ImageType, ImageType >   ResampleType;
  typedef typename ResampleType::Pointer 								ResamplePointer;

  typedef ExtractImageFilter< SequenceType, ImageType > ImageExtractorType;
  typedef typename ImageExtractorType::Pointer 					ImageExtractorPointer;

  typedef ImageFileWriter< ImageType >  ImageWriterType;

  typedef  TransformType::ParametersType    ParametersType;

  typedef DiscreteGaussianImageFilter< ImageType, ImageType > GaussianFilterType;
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

  typedef vnl_matrix<double> VnlMatrixType;
  typedef vnl_vector<double> VnlVectorType;

  void SetGradientTable( const char* input );
//  void WriteGradientTable( const char* output );
//  void RotateGradients( );

  /** Method that initiates the registration. */
  void StartRegistration();

    /** Set/Get the original sequence (to be corrected). */
  itkSetObjectMacro( Input, SequenceType );
  itkGetObjectMacro( Input, SequenceType );

   /** Get the corrected sequence. */
//  itkGetObjectMacro( Output, SequenceType );

  /** Set/Get the FixedImageRegion. */
  TSequence * GetOutput();

  /** Set/Get the maximum number of iterations. */
  itkSetMacro( Iterations, unsigned int );
  itkGetMacro( Iterations, unsigned int );

  /** Set/Get the FixedImageRegion. */
  void SetFixedImageRegion( const  RegionType & region );

  /** Set/Get the image mask. */
  itkSetObjectMacro( ImageMask, ImageMaskType );
  itkGetObjectMacro( ImageMask, ImageMaskType );

   /** Set/Get the mean dw-image. */
  virtual ImageType * GetMeanGradient()
  {
    return this -> m_Resample -> GetOutput();
  }

   /** Set/Get the image array. */
//  UserGetObjectMacro( TransformArray, TransformType );

   /** Write transforms ( DW -> T2). */
  void WriteTransforms( const char* folder );

  void Initialize();


protected:
  GroupwiseS2SDistortionCorrection();
  virtual ~GroupwiseS2SDistortionCorrection() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */




private:
  GroupwiseS2SDistortionCorrection(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  SequencePointer									 m_Input;
  SequencePointer									 m_Output;

  unsigned int                     m_GradientDirections;
  unsigned int                     m_Iterations;

  ImageMaskPointer								 m_ImageMask;
  ImagePointer										 m_FixedImage;
  ImagePointer                     m_T2epi;
  ImagePointer										 m_MeanGradient;

  float  *m_OriginalGradients;
  double *m_CorrectedGradients;

  S2STransformArray                m_TransformArray;
//  TransformPointerArray            m_NormalizedTransformArray;
  InterpolatorPointer			         m_Interpolator;
  ImageArrayPointer 			         m_ImageArray;
  RegistrationPointer							 m_Registration;
  ResamplePointer                  m_Resample;

  RegionType 											 m_FixedImageRegion;
  bool                             m_FixedImageRegionDefined;

  SequenceRegionType 							 m_SequenceRegion;
  SequenceSizeType 							 	 m_SequenceSize;
  SequenceIndexType 							 m_SequenceIndex;

  AffineRegistrationPointer m_AffineRegistration;

  vnl_matrix< double > m_GradientTable;

};


} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkGroupwiseS2SDistortionCorrection.txx"
#endif

#endif
