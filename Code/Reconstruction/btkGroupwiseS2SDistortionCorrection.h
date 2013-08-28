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

/** @class GroupwiseS2SDistortionCorrection
 * @brief It performs a groupwise registration by using slice by slice transforms.
 *
 * This class performs a groupwise registration of a set of images by using slice
 * by slice affine transforms.
 *
 * @author Estanislao Oubel
 * @ingroup Reconstruction
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

  /** Type of the sequence. */
  typedef          TSequence                            SequenceType;
  typedef typename SequenceType::Pointer                SequencePointer;

  /** Sequence pixel typedef support. */
  typedef typename SequenceType::PixelType              SequencePixelType;

  /** Sequence region typedef support. */
  typedef typename SequenceType::RegionType             SequenceRegionType;

  /** Sequence size typedef support. */
  typedef typename SequenceType::SizeType             	SequenceSizeType;

  /** Sequence index typedef support. */
  typedef typename SequenceType::IndexType             	SequenceIndexType;

  /** Type of the sequence iterator. */
  typedef ImageRegionIteratorWithIndex< SequenceType > SequenceIteratorType;

  /** Type of the image. */
  typedef          Image< SequencePixelType, 3 >        ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;

  /** Image pixel typedef support. */
  typedef typename ImageType::PixelType                 PixelType;

  /** Image point typedef support. */
  typedef typename ImageType::PointType                 PointType;

  /** Image spacing typedef support. */
  typedef typename ImageType::SpacingType               SpacingType;

  /** Image size typedef support. */
  typedef typename ImageType::SizeType               		SizeType;

  /** Image index typedef support. */
  typedef typename ImageType::IndexType               	IndexType;

  /** Image region typedef support. */
  typedef typename ImageType::RegionType               	RegionType;

  /** Type of the image pointer array. */
  typedef          std::vector<ImagePointer>            ImageArrayPointer;

  /** Type of the image mask. */
  typedef          Image< unsigned char, 3 >       			ImageMaskType;
  typedef typename ImageMaskType::Pointer               ImageMaskPointer;

  /** Type of the mask. */
  typedef ImageMaskSpatialObject< 3 >   								MaskType;
  typedef typename MaskType::Pointer               			MaskPointer;

  /** Type of the image iterator. */
  typedef ImageRegionIteratorWithIndex< ImageType >  IteratorType;

  /** Type of the Transform . */
  typedef AffineTransform< double, 3>           				TransformType;
  typedef typename TransformType::Pointer               TransformPointer;

  //TODO This should be replaced by the use of the affine version of
  // btkSliceBySliceTransform (to be created)
  /** Type of the slice by slice transform. */
  typedef  std::vector<TransformPointer>                S2STransformType;

  /** Type of the slice by slice transform array. */
  typedef  std::vector<S2STransformType>                S2STransformArray;

  /** Type of the transform writer. */
  typedef TransformFileWriter                           TransformWriterType;

  /** Type of the transform matrix. */
  typedef typename TransformType::MatrixType            MatrixType;

  /** Type of the transform offset . */
  typedef typename TransformType::OffsetType            OffsetType;

  /** Type of the transform parameters. */
  typedef  TransformType::ParametersType    ParametersType;

  /** Type of the RBF interpolator. */
  typedef SliceToSliceInterpolateImageFunction<
                                      ImageType,
                                      double >          RBFInterpolatorType;
  typedef typename RBFInterpolatorType::Pointer         RBFInterpolatorPointer;

  /** Type of the linear interpolator. */
  typedef LinearInterpolateImageFunction<
                                    ImageType,
                                    double>             InterpolatorType;
  typedef typename InterpolatorType::Pointer            InterpolatorPointer;

  /** Type of the registration method. */
  typedef SliceBySliceRegistration< ImageType >   			RegistrationType;
  typedef typename RegistrationType::Pointer            RegistrationPointer;

  /** Type of affine registration method. */
  typedef AffineRegistration<ImageType>                 AffineRegistrationType;
  typedef typename AffineRegistrationType::Pointer      AffineRegistrationPointer;

  /** Type of the image resampler. */
  typedef ResampleImageFilter< ImageType, ImageType >   ResampleType;
  typedef typename ResampleType::Pointer 								ResamplePointer;

  /** Type of the image extractor. */
  typedef ExtractImageFilter< SequenceType, ImageType > ImageExtractorType;
  typedef typename ImageExtractorType::Pointer 					ImageExtractorPointer;

  /** Type of the image writer. */
  typedef ImageFileWriter< ImageType >  ImageWriterType;

  /** Type of the Gaussian filter. */
  typedef DiscreteGaussianImageFilter< ImageType, ImageType > GaussianFilterType;
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;

  /** Type of the Vnl matrix. */
  typedef vnl_matrix<double> VnlMatrixType;

  /** Type of the Vnl vector. */
  typedef vnl_vector<double> VnlVectorType;

  //TODO The library has now a class called btkDiffusionGradientTable to manage
  // gradient tables. We could set the gradient information by setting an object
  // of this type, instead of providing the file name.
  /** Sets the gradient table (file name). */
  void SetGradientTable( const char* input );
//  void WriteGradientTable( const char* output );
//  void RotateGradients( );

  /** Method that initiates the registration. */
  void StartRegistration();

  /** Sets the input sequence. */
  itkSetObjectMacro( Input, SequenceType );

  /** Gets the input sequence. */
  itkGetObjectMacro( Input, SequenceType );

  /** Get the corrected sequence. */
  TSequence * GetOutput();

  /** Sets the maximum number of iterations. */
  itkSetMacro( Iterations, unsigned int );

  /** Gets the maximum number of iterations. */
  itkGetMacro( Iterations, unsigned int );

  /** Sets the FixedImageRegion. */
  void SetFixedImageRegion( const  RegionType & region );

  /** Sets the image mask. */
  itkSetObjectMacro( ImageMask, ImageMaskType );

  /** Gets the image mask. */
  itkGetObjectMacro( ImageMask, ImageMaskType );

   /** Gets the mean dw-image. */
  virtual ImageType * GetMeanGradient()
  {
    return m_MeanGradient;
  }

  /** Gets the mean gradient registered to B0. */
  itkGetObjectMacro( RegisteredMeanGradient, ImageType);

  /** Write transforms ( DW -> T2) to the specified location. */
  void WriteTransforms( const char* folder );

  /** Performs the class initialization.*/
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
  ImagePointer										 m_RegisteredMeanGradient;

  float  *m_OriginalGradients;
  double *m_CorrectedGradients;

  S2STransformArray                m_TransformArray;
  InterpolatorPointer			         m_Interpolator;
  ImageArrayPointer 			         m_ImageArray;
  RegistrationPointer							 m_Registration;

  RegionType 											 m_FixedImageRegion;
  bool                             m_FixedImageRegionDefined;

  SequenceRegionType 							 m_SequenceRegion;
  SequenceSizeType 							 	 m_SequenceSize;
  SequenceIndexType 							 m_SequenceIndex;

  AffineRegistrationPointer 			 m_AffineRegistration;

  vnl_matrix< double > m_GradientTable;

};


} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkGroupwiseS2SDistortionCorrection.txx"
#endif

#endif
