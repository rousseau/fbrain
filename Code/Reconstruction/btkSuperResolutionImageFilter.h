/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 02/12/2010
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

#ifndef __btkSuperResolutionImageFilter_h
#define __btkSuperResolutionImageFilter_h

#include "itkFixedArray.h"
#include "itkTransform.h"
#include "itkEuler3DTransform.h"
#include "itkImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageToImageFilter.h"
#include "itkSize.h"
#include "btkUserMacro.h"
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/algo/vnl_conjugate_gradient.h"
#include "btkLeastSquaresVnlCostFunction.h"
#include "itkImageDuplicator.h"
#include "itkContinuousIndex.h"

namespace btk
{

using namespace itk;

/** \class SuperResolutionImageFilter
 * \brief Super-resolution by using a set of low resolution images.and the
 * reconstructed image.
 *
 * SuperResolutionImageFilter allows to obtain a super-resolution image from a
 * set of low-resolution images and the reconstructed image. The class is templated
 * over the types of the input and output images.
 *
 * The implemented method is similar to the one described in:
 *
 * F. Rousseau,  K. Kim,  C. Studholme,  M. Koob,  J.-L. Dietemann
 * On Super-Resolution for Fetal Brain MRI, Medical Image Computing and Computer
 * Assisted Intervention Pékin, Chine, pp. 355--362, Springer-Verlag, Lecture
 * Notes in Computer Science, Vol. 6362, doi:10.1007/978-3-642-15745-5, 2010
 *
 * \ingroup Reconstruction
 */

template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType=double>
class SuperResolutionImageFilter:
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef enum {
    BOXCAR=0,
    GAUSSIAN=1,
  } PSF_type;

  /** Standard class typedefs. */
  typedef SuperResolutionImageFilter                Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                             InputImageType;
  typedef TOutputImage                            OutputImageType;
  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename OutputImageType::Pointer       OutputImagePointer;

  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef std::vector<InputImageRegionType>       InputImageRegionVectorType;

  typedef ImageMaskSpatialObject< TInputImage::ImageDimension > MaskType;
  typedef typename MaskType::Pointer   MaskPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SuperResolutionImageFilter, ImageToImageFilter);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Transform typedef. */
  typedef Euler3DTransform<TInterpolatorPrecisionType> TransformType;
  typedef typename TransformType::Pointer TransformPointerType;

  //TODO This should be replaced by a std::vector of btkSliceBySliceTransform.
  /** Type of the transform list. */
  typedef std::vector< std::vector<TransformPointerType> > TransformPointerArrayType;

  /** Image size typedef. */
  typedef Size<itkGetStaticConstMacro(ImageDimension)> SizeType;

  /** Image index typedef support. */
  typedef typename TOutputImage::IndexType IndexType;

  /** Image point typedef support. */
  typedef typename TOutputImage::PointType    PointType;

  /** Image pixel typedef support. */
  typedef typename TOutputImage::PixelType   PixelType;

  /** Input image pixel typedef support. */
  typedef typename TInputImage::PixelType    InputPixelType;

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Image spacing typedef support. */
  typedef typename TOutputImage::SpacingType   SpacingType;

  /** Image direction typedef support. */
  typedef typename TOutputImage::DirectionType DirectionType;

  /** base type for images of the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

  /**Const iterator typedef. */
  typedef ImageRegionConstIteratorWithIndex< OutputImageType >  OutputIteratorType;

  /**VnlVectorType typedef. */
  typedef vnl_vector<double> VnlVectorType;

  /** Overrides SetInput to resize the transform. */
  void AddInput(InputImageType* _arg);

  /** Set the transform array. */
  void SetTransform( int i, int j, TransformType* transform )
  {
    m_Transform[i][j] = transform;
  }

  /** Converts from a linear index (li = i+i*x_size+k*x_size*y_size) to an absolute
   * index (ai = [i j k]). */
  IndexType LinearToAbsoluteIndex( unsigned int linearIndex, InputImageRegionType region );

  /** Set the size of the output image. */
  itkSetMacro( Size, SizeType );

  /** Get the size of the output image. */
  itkGetConstReferenceMacro( Size, SizeType );

  /** Set the pixel value when a transformed pixel is outside of the
   * image.  The default default pixel value is 0. */
  itkSetMacro( DefaultPixelValue, PixelType );

  /** Get the pixel value when a transformed pixel is outside of the image */
  itkGetConstReferenceMacro( DefaultPixelValue, PixelType );

  /** Set the output image spacing. */
  itkSetMacro( OutputSpacing, SpacingType );

  /** Set the output image spacing as a const array of values. */
  virtual void SetOutputSpacing( const double* values );

  /** Get the output image spacing. */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  /** Set the output image origin. */
  itkSetMacro( OutputOrigin, PointType );

  /** Set the output image origin as a const array of values. */
  virtual void SetOutputOrigin( const double* values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro( OutputOrigin, PointType );

  /** Set the output direciton cosine matrix. */
  itkSetMacro( OutputDirection, DirectionType );
  itkGetConstReferenceMacro( OutputDirection, DirectionType );

  /** Helper method to set the output parameters based on this image */
  void SetOutputParametersFromImage ( const ImageBaseType * image );

  /** Set the start index of the output largest possible region.
   * The default is an index of all zeros. */
  itkSetMacro( OutputStartIndex, IndexType );

  /** Get the start index of the output largest possible region. */
  itkGetConstReferenceMacro( OutputStartIndex, IndexType );

  /** Copy the output information from another Image.  By default,
   *  the information is specified with the SetOutputSpacing, Origin,
   *  and Direction methods. UseReferenceImage must be On and a
   *  Reference image must be present to override the default behavior.
   *  NOTE: This function seems redundant with the
   *  SetOutputParametersFromImage( image ) function */
  void SetReferenceImage ( const TOutputImage *image );

  /** Gets the reference image. */
  const TOutputImage * GetReferenceImage( void ) const;

  /** Sets the use of a reference image to true/false. */
  itkSetMacro( UseReferenceImage, bool );

  /** Adds UseReferenceImageOff/On. */
  itkBooleanMacro( UseReferenceImage );

  /** Gets the status of the UseReferenceImage variable. */
  itkGetConstMacro( UseReferenceImage, bool );

  /** SuperResolutionImageFilter produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation( void );

  /** SuperResolutionImageFilter needs a different input requested region than
   * the output requested region.  As such, SuperResolutionImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion( void );

  /** Method Compute the Modified Time based on changed to the components. */
  unsigned long GetMTime( void ) const;

  /** Adds an image region. The regions must be added in the same order than the
   * input images.*/
  void AddRegion(InputImageRegionType _arg)
  {
    m_InputImageRegion.push_back(_arg);
  }

  /** Adds an image mask. Masks must be added in the same order than the
   * input images.*/
  void AddMask(MaskType *mask)
  {
    m_MaskArray.push_back( mask );
  }

  /** Sets the output image region.*/
  itkSetMacro(OutputImageRegion, OutputImageRegionType);

  /** Gets the output image region.*/
  itkGetMacro(OutputImageRegion, OutputImageRegionType);

  /** Sets the number of iterations.*/
  itkSetMacro(Iterations, unsigned int);

  /** Gets the number of iterations.*/
  itkGetMacro(Iterations, unsigned int);

  /** Sets the lambda value for regularization.*/
  itkSetMacro(Lambda, float);

  /** Gets the lambda value for regularization.*/
  itkGetMacro(Lambda, float);

  /** Sets the type of PSF (Boxcar/Gaussian).*/
  itkSetMacro(PSF, unsigned int);

  /** Gets the type of PSF (Boxcar/Gaussian).*/
  itkGetMacro(PSF, unsigned int);


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<PixelType>));
  /** End concept checking */
#endif

protected:
  SuperResolutionImageFilter( void );
  ~SuperResolutionImageFilter( void ) {};

  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** ResampleImageFilter cannot be implemented as a multithreaded filter.  Therefore,
   * this implementation only provides a GenerateData() routine that allocates output
   * image data.
   * */
  void GenerateData();

private:

  SuperResolutionImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void OptimizeByLeastSquares();

  SizeType                    m_Size;         /**< Size of the output image. */
  TransformPointerArrayType   m_Transform;
  InputImageRegionVectorType  m_InputImageRegion;

  OutputImageRegionType       m_OutputImageRegion;

  std::vector<InputImagePointer>  m_ImageArray;
  std::vector<MaskPointer>  			m_MaskArray;

  VnlVectorType     m_x;

  float m_Lambda;

  unsigned int 			m_Iterations;

  PixelType         m_DefaultPixelValue; /**< Default pixel value if the point
                                              falls outside the image. */
  SpacingType       m_OutputSpacing;     /**< Spacing of the output image. */
  PointType   			m_OutputOrigin;      /**< Origin of the output image. */
  DirectionType     m_OutputDirection;   /**< Direction of the output image. */
  IndexType         m_OutputStartIndex;  /**< Start index of the output image.*/
  bool              m_UseReferenceImage;

  unsigned int m_PSF;

};


} // end namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSuperResolutionImageFilter.txx"
#endif

#endif
