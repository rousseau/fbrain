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
#include "itkNeighborhoodAlgorithm.h"
#include "itkSize.h"
#include "btkUserMacro.h"
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/algo/vnl_conjugate_gradient.h"
#include "btkLeastSquaresVnlCostFunction.h"
#include "btkLinearInterpolateImageFunctionWithWeights.h"
#include "btkOrientedSpatialFunction.h"
#include "itkImageDuplicator.h"
#include "itkContinuousIndex.h"

namespace btk
{

using namespace itk;

/** \class SuperResolutionImageFilter
 * \brief Resample an image via a coordinate transform
 *
 * SuperResolutionImageFilter resamples an existing image through some coordinate
 * transform, interpolating via some image function.  The class is templated
 * over the types of the input and output images.
 *
 * Output information (spacing, size and direction) for the output
 * image should be set. This information has the normal defaults of
 * unit spacing, zero origin and identity direction. Optionally, the
 * output information can be obtained from a reference image. If the
 * reference image is provided and UseReferenceImage is On, then the
 * spacing, origin and direction of the reference image will be used.
 *
 * Since this filter produces an image which is a different size than
 * its input, it needs to override several of the methods defined
 * in ProcessObject in order to properly manage the pipeline execution model.
 * In particular, this filter overrides
 * ProcessObject::GenerateInputRequestedRegion() and
 * ProcessObject::GenerateOutputInformation().
 *: 105
 *
 * \ingroup GeometricTransforms
 */

template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType=double>
class SuperResolutionImageFilter:
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:

  typedef enum {
    MSE=0,
    BACKPROJECTION=1,
  } OptimizationMethods;

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

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SuperResolutionImageFilter, ImageToImageFilter);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);


  /** Transform typedef. */
  typedef Euler3DTransform<TInterpolatorPrecisionType> TransformType;
  typedef typename TransformType::Pointer TransformPointerType;

  typedef std::vector< std::vector<TransformPointerType> > TransformPointerArrayType;

  /** Image size typedef. */
  typedef Size<itkGetStaticConstMacro(ImageDimension)> SizeType;

  /** Image index typedef. */
  typedef typename TOutputImage::IndexType IndexType;

  /** Image continuous index typedef. */
  typedef ContinuousIndex<double, ImageDimension> ContinuousIndexType;

  /** Image point typedef. */
  typedef typename TOutputImage::PointType    PointType;

  /** Image pixel value typedef. */
  typedef typename TOutputImage::PixelType   PixelType;
  typedef typename TInputImage::PixelType    InputPixelType;

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Image spacing, origin and direction typedef */
  typedef typename TOutputImage::SpacingType   SpacingType;
  typedef typename TOutputImage::PointType     OriginPointType;
  typedef typename TOutputImage::DirectionType DirectionType;

  /** base type for images of the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

  /**Const iterator typedef. */
  typedef ImageRegionConstIteratorWithIndex< InputImageType >  ConstIteratorType;

  /**Const iterator typedef. */
  typedef ImageRegionConstIteratorWithIndex< OutputImageType >  OutputIteratorType;

  /**Neighborhood iterator typedef. */
  typedef NeighborhoodIterator< InputImageType >   NeighborhoodIteratorType;

  /**Blurring function typedef. */
  typedef OrientedSpatialFunction<double, 3, PointType> FunctionType;

  /**Interpolator typedef. */
  typedef LinearInterpolateImageFunctionWithWeights<TOutputImage, double> InterpolatorType;

  /**Interpolator typedef. */
//  typedef NearestNeighborInterpolateImageFunction<TInputImage, double> NNInterpolatorType;

  /** Duplicator typedef. */
  typedef ImageDuplicator< InputImageType > DuplicatorType;

  /**Face calculator typedef. */
  typedef NeighborhoodAlgorithm
  ::ImageBoundaryFacesCalculator< OutputImageType > FaceCalculatorType;
  typedef typename FaceCalculatorType::FaceListType FaceListType;

  typedef vnl_matrix<float> VnlMatrixType;
  typedef itk::Matrix<float,3,3> MatrixType;
  typedef vnl_vector<float> VnlVectorType;

  typedef vnl_sparse_matrix<float> VnlSparseMatrixType;

  // Overrides SetInput to resize the transform

  void AddInput(InputImageType* _arg)
  {
    m_ImageArray.push_back(_arg);

    this -> SetInput(_arg);

    m_Transform.resize( m_Transform.size() + 1 );
    SizeType _argSize = _arg -> GetLargestPossibleRegion().GetSize();
    m_Transform[m_Transform.size()-1].resize(_argSize[2]);

    typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
    duplicator -> SetInputImage (_arg);
    duplicator -> Update();
    m_SimulatedImages.push_back( duplicator -> GetOutput() );
    m_SimulatedImages[m_SimulatedImages.size()-1]->FillBuffer(0);

  }

  /** Set the transform array. */
  void SetTransform( int i, int j, TransformType* transform )
  {
    m_Transform[i][j] = transform;
  }

  /** Get the simulated image */
  InputImageType* GetSimulatedImage( int i )
  {
    return m_SimulatedImages[i];
  }

  IndexType LinearToAbsoluteIndex( unsigned int linearIndex, InputImageRegionType region );

  void UpdateSimulatedImages();


  /** Get the transform array. */
//  TransformType* GetTransform( int i)
//  {
//    return this -> m_Transform[i].GetPointer();
//  }

  /** Get the transform array. */
//  itkGetMacroNoDeb( Transform, TransformPointerArrayType );

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
  virtual void SetOutputSpacing( const double* values );

  /** Get the output image spacing. */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  /** Set the output image origin. */
  itkSetMacro( OutputOrigin, OriginPointType );
  virtual void SetOutputOrigin( const double* values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro( OutputOrigin, OriginPointType );

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
   *  Reference image must be present to override the defaul behavior.
   *  NOTE: This function seems redundant with the
   *  SetOutputParametersFromImage( image ) function */
  void SetReferenceImage ( const TOutputImage *image );
  const TOutputImage * GetReferenceImage( void ) const;

  itkSetMacro( UseReferenceImage, bool );
  itkBooleanMacro( UseReferenceImage );
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

  void AddRegion(InputImageRegionType _arg)
  {
    m_InputImageRegion.push_back(_arg);

  }

  // Set/Get output image region
  itkSetMacro(OutputImageRegion, OutputImageRegionType);
  itkGetMacro(OutputImageRegion, OutputImageRegionType);

  // Set number of iterations
  itkSetMacro(Iterations, unsigned int);

  // Set/Get output image region
  itkSetMacro(OptimizationMethod, unsigned int);
  itkGetMacro(OptimizationMethod, unsigned int);

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

/*  class vnl_my_cost_fun : public vnl_cost_function
  {
    public: vnl_my_cost_fun(): vnl_cost_function(1) {}

    double f(const vnl_vector<double>& x)
    {
      std::cout << m_x.size();
      return (x[0]-5)*(x[0]-5)+10;
    }

    void gradf(const vnl_vector<double>& x, vnl_vector<double>& g)
    {
      g[0] = 2 *x[0]-10;
    }
  };*/

  SuperResolutionImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  void CreateH();
  void OptimizeByBackprojection();
  void OptimizeByLeastSquares();

  SizeType                    m_Size;              // Size of the output image
  TransformPointerArrayType   m_Transform;         // Coordinate transform to use
  InputImageRegionVectorType  m_InputImageRegion;

  OutputImageRegionType       m_OutputImageRegion;

  std::vector<InputImagePointer>     m_ImageArray;
  std::vector<InputImagePointer>     m_SimulatedImages;
  bool    m_SimulatedImagesUpdated;

  std::vector< VnlSparseMatrixType > m_H;
  std::vector< VnlSparseMatrixType > m_Ht;
  std::vector< VnlSparseMatrixType > m_Hbp;
  std::vector< VnlVectorType >       m_y;
  std::vector< VnlVectorType >       m_ysim;
  VnlVectorType                      m_x;

  // Precomputed values for optimization

  VnlSparseMatrixType HtH;
  VnlVectorType HtY;
  double YtY;

  unsigned int m_Iterations;

  PixelType                   m_DefaultPixelValue; // default pixel value
                                               // if the point is
                                               // outside the image
  SpacingType                 m_OutputSpacing;     // output image spacing
  OriginPointType             m_OutputOrigin;      // output image origin
  DirectionType               m_OutputDirection;   // output image direction cosines
  IndexType                   m_OutputStartIndex;  // output image start index
  bool                        m_UseReferenceImage;

  unsigned int m_OptimizationMethod;

};


} // end namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSuperResolutionImageFilter.txx"
#endif

#endif
