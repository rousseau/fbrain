/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/03/2010 modified 04/06/2013
  Author(s): Estanislao Oubel (oubel@unistra.fr)
             Frederic Champ (champ(at)unistra.fr)

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

#ifndef __btkRBFInterpolateImageFunctionS2S_h
#define __btkRBFInterpolateImageFunctionS2S_h


// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

#include "itkTransformFileWriter.h"
#include "itkTransformFileReader.h"

//#include "nr3.h"
//#include "ludcmp.h"
//#include "interp_rbf.h"

#include "btkRBFInterpolation.h"
#include "btkGradientDirection.h"
#include "ANN.h"



namespace btk
{

using namespace itk;

/** \class LinearInterpolateImageFunction
 * \brief Linearly interpolate an image at specified positions.
 *
 * LinearInterpolateImageFunction linearly interpolates image intensity at
 * a non-integer pixel position. This class is templated
 * over the input image type and the coordinate representation type
 * (e.g. float or double).
 *
 * This function works for N-dimensional images.
 *
 * \warning This function work only for images with scalar pixel
 * types. For vector images use VectorLinearInterpolateImageFunction.
 *
 * \sa VectorLinearInterpolateImageFunction
 *
 * \ingroup ImageFunctions ImageInterpolators
 */
template <class TInputImage, class TCoordRep = double>
class RBFInterpolateImageFunctionS2S :
  public InterpolateImageFunction<TInputImage,TCoordRep>
{
public:
  /** Standard class typedefs. */
  typedef RBFInterpolateImageFunctionS2S                  Self;
  typedef InterpolateImageFunction<TInputImage,TCoordRep> Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  typedef ImageRegionConstIteratorWithIndex<TInputImage> IteratorType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(RBFInterpolateImageFunctionS2S, InterpolateImageFunction);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** OutputType typedef support. */
  typedef typename Superclass::OutputType OutputType;

  /** InputImageType typedef support. */
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename InputImageType::SizeType   ImageSizeType;
  typedef typename InputImageType::RegionType ImageRegionType;
  typedef typename InputImageType::IndexType  ImageIndexType;
  typedef typename InputImageType::PointType ImagePointType;
  typedef typename InputImageType::ConstPointer ImageConstPointer;
  typedef typename InputImageType::Pointer ImagePointer;

  typedef typename Superclass::PointType PointType;

  /** InputPixelType typedef support. */
  typedef typename Superclass::InputPixelType InputPixelType;

  /** RealType typedef support. */
  typedef typename Superclass::RealType RealType;

  /** Dimension underlying input image. */
  itkStaticConstMacro(ImageDimension, unsigned int,Superclass::ImageDimension);

  /** Index typedef support. */
  typedef typename Superclass::IndexType IndexType;

  /** ContinuousIndex typedef support. */
  typedef typename Superclass::ContinuousIndexType ContinuousIndexType;

  typedef AffineTransform< double, 3 >              TransformType;
  typedef TransformType::Pointer                    TransformPointer;
  typedef  std::vector<TransformPointer>            S2STransformType;
  typedef  std::vector<S2STransformType>            S2STransformArray;

  typedef Euler3DTransform< double >           RigidTransformType;

  /** Transform reader type. */
  typedef itk::TransformFileReader     TransformReaderType;
  typedef TransformReaderType::TransformListType * TransformListType;

  typedef vnl_matrix<double> VnlMatrixType;
  typedef itk::Matrix<double,3,3> MatrixType;
  typedef vnl_vector<double> VnlVectorType;

  /** Evaluate the function at a ContinuousIndex position
   *
   * Returns the linearly interpolated image intensity at a
   * specified point position. No bounds checking is done.
   * The point is assume to lie within the image buffer.
   *
   * ImageFunction::IsInsideBuffer() can be used to check bounds before
   * calling the method. */
  virtual OutputType EvaluateAtContinuousIndex(
    const ContinuousIndexType & index ) const;

  void SetGradientTable( const char* input );

  void SetGradientTable( std::vector<btk::GradientDirection> input );

  void RotateGradients( );

  void GetGradientDirection ( unsigned int index, double &theta, double &phi );

  virtual OutputType Evaluate( const PointType& point,
                               double theta, double phi,
                               double r_spa, double r_gra,
                               char init) const;

  virtual OutputType EvaluateAt( ImageIndexType index,
                               double theta, double phi,
                               double r_spa, double r_gra,
                               char init) const;

  void SetInputImage(const InputImageType *ptr);
  void SetInputImage(ImagePointer ptr);

  // The region used as argument for initialize is a region in the original
  // sequence. It's used to construct a tree with the points used for
  // interpolation
  void Initialize( ImageRegionType & region );
  void SetTransforms( const char * tpath );


protected:
  RBFInterpolateImageFunctionS2S();
  ~RBFInterpolateImageFunctionS2S();
  void PrintSelf(std::ostream& os, Indent indent) const;


private:
  RBFInterpolateImageFunctionS2S( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  /** Number of neighbors used in the interpolation */
  static const unsigned long  m_Neighbors;

  vnl_matrix< double > m_OriginalGradientTable;
  vnl_matrix< double > m_GradientTableCartesian;
//  vnl_matrix< double > m_GradientTableSpherical;

  ImageSizeType m_ImageSize;

//  unsigned int m_Sigma;
  unsigned int m_NumberOfGradients;

  unsigned int m_NumberOfSlices;
  unsigned int m_FirstSlice;
  unsigned int m_LastSlice;

  bool m_TransformsAreSet;

  ANNpointArray m_dataPtsSphere;
  ANNkd_tree*   m_kdTreeSphere;

  std::vector<ANNpointArray> m_dataPtsSpace;
  std::vector<ANNkd_tree*> m_kdTreeSpace;
  std::vector<ANNpointArray> m_dataVals;

  S2STransformArray                m_TransformArray;

//  TransformPointerArray m_InverseTransformArray;

//  mutable RBF_interp *m_RBFInterpolator;

};

} // end namespace btk

// Define instantiation macro for this template.
#define ITK_TEMPLATE_RBFInterpolateImageFunctionS2S(_, EXPORT, x, y) namespace itk { \
  _(2(class EXPORT RBFInterpolateImageFunctionS2S< ITK_TEMPLATE_2 x >)) \
  namespace Templates { typedef RBFInterpolateImageFunctionS2S< ITK_TEMPLATE_2 x > \
                                                  RBFInterpolateImageFunctionS2S##y; } \
  }

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/btkRBFInterpolateImageFunctionS2S+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "btkRBFInterpolateImageFunctionS2S.txx"
#endif

#endif
