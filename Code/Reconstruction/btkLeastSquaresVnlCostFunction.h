/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 28/04/2011
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

#ifndef BTKLEASTSQUARESVNLCOSTFUNCTION_H_
#define BTKLEASTSQUARESVNLCOSTFUNCTION_H_

#include "vnl/vnl_cost_function.h"
#include "vnl/vnl_matops.h"
#include "btkLinearInterpolateImageFunctionWithWeights.h"
#include "btkOrientedSpatialFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"

namespace btk
{
/** \class LeastSquaresVnlCostFunction
 * \brief Mean squared error function class to be used with vnl_conjugate_gradient
 *
 * LeastSquaresVnlCostFunction implements the mse cost function and its gradient
 * as required to be used with vnl_conjugate_gradient. H, Y, and HtY are stored
 * in this class to save memory.
 *
 * \ingroup Reconstruction
 */
template <class TImage>
class LeastSquaresVnlCostFunction : public vnl_cost_function
{
  public:

  typedef TImage   ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;

  typedef typename ImageType::RegionType  RegionType;
  typedef typename ImageType::SizeType    SizeType;
  typedef typename ImageType::IndexType   IndexType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::PointType   PointType;

  typedef ImageMaskSpatialObject< TImage::ImageDimension > MaskType;
  typedef typename MaskType::Pointer   MaskPointer;

  typedef Euler3DTransform<double> TransformType;
  typedef typename TransformType::Pointer TransformPointerType;

  typedef LinearInterpolateImageFunctionWithWeights<ImageType, double> InterpolatorType;
  typedef typename InterpolatorType::Pointer InterpolatorPointer;

  typedef ContinuousIndex<double, TImage::ImageDimension> ContinuousIndexType;

  /**Oriented spatial function typedef. */
  typedef OrientedSpatialFunction<double, 3, PointType> FunctionType;

  /**Const iterator typedef. */
  typedef ImageRegionConstIteratorWithIndex< ImageType >  ConstIteratorType;

  typedef vnl_vector<float> VnlVectorType;
  typedef vnl_sparse_matrix<float> VnlSparseMatrixType;

  struct img_size {
     unsigned int width;
     unsigned int height;
     unsigned int depth;
  } x_size;

  LeastSquaresVnlCostFunction(unsigned int dim);

  /**Cost function (MSE). */
  double f(const vnl_vector<double>& x);

  /**Gradient of the cost function, required for optimization with
    vnl_conjugate_gradient . */
  void gradf(const vnl_vector<double>& x, vnl_vector<double>& g);

  /**Construction of vector Y and matrix H of the observation model Y=H*X. */
  void Initialize();

  /**Sets the weight of the regularization term.*/
  void SetLambda(float value);

  /**Adds a low-resolution image.*/
  void AddImage( ImageType* image );

  /**Adds the image region of the corresponding low-resolution image.*/
  void AddRegion( RegionType region);

  /**Adds the image mask of the corresponding low-resolution image.*/
  void AddMask( MaskType *mask);

  /**Sets the reference image, i.e. an image providing the spatial positions
  where to compute the intensity values (each X value). This image is expected
  to have isotropic voxels, with a size equal to the in-plane voxel size of the
  low resolution images. However, any image can be provided to this end, the
  code has been written generically. */
  void SetReferenceImage( const ImageType * image );

  // TODO This function must be modified after modifying this class to use slice
  // by slice transforms.
  /**Sets the transforms obtained with the reconstruction method. These
  transformations correct the movements of the fetus during image acquisition.*/
  void SetTransform( int i, int j, TransformType* transform );

  /** Sets the type of PSF (Boxcar, Gaussian). */
  void SetPSF(unsigned int psf)
  {
    m_PSF = psf;
  }

  /** Gets the type of PSF (Boxcar, Gaussian). */
  itkGetMacro(PSF, unsigned int);


  private:

  void set(float * array, int size, float value);

  /** Mirror of the position pos. abs(pos) must not be > 2*(size-1) */
  int mirror(int pos, int size);

  /** Gets an image row as a vector. */
  void get_row(const vnl_vector<float>& image, img_size & size, int row,
      int frame, float * output);

  /** Sets an image row from a vector. */
  void set_row(vnl_vector<float>& image, img_size & size, int row, int frame,
      float * input);

  /** Gets an image column as a vector. */
  void get_col(const vnl_vector<float>& image, img_size & size, int col,
      int frame, float * output);

  /** Sets an image column from a vector. */
  void set_col(vnl_vector<float>& image, img_size & size, int col, int frame,
      float * input);

  /** Gets an image z-axis as a vector. */
  void get_spec(const vnl_vector<float>& image, img_size & size, int row,
      int col, float * output);

  /** Sets an image z-axis from a vector. */
  void set_spec(vnl_vector<float>& image, img_size & size, int row, int col,
      float * input);

  void convol1d(float * kernel, int ksize, float * src, int src_size, float * dest);

  // 3D convolution : over the rows
  void convol3dx(const vnl_vector<float>& image, vnl_vector<float>& image_conv,
      img_size& size, float * kernel, int ksize);

  // 3D convolution : over the columns
  void convol3dy(const vnl_vector<float>& image, vnl_vector<float>& image_conv,
      img_size & size, float * kernel, int ksize);

  // 3D convolution : over the spectra
  void convol3dz(const vnl_vector<float>& image, vnl_vector<float>& image_conv,
      img_size & size, float * kernel, int ksize);

  // Gets the first derivative of an image along the x axis.
  void get_derivative_x(const vnl_vector<double>& x, vnl_vector<float>& deriv_x);

  // Gets the first derivative of an image along the y axis.
  void get_derivative_y(const vnl_vector<double>& x, vnl_vector<float>& deriv_y);

  // Gets the first derivative of an image along the z axis.
  void get_derivative_z(const vnl_vector<double>& x, vnl_vector<float>& deriv_z);

  vnl_sparse_matrix<float> H;
  vnl_vector<float> HtY;
  vnl_vector<float> Y;

  float lambda;

  std::vector<ImagePointer>  m_Images;
  std::vector<RegionType>    m_Regions;
  std::vector<MaskPointer>   m_Masks;
  ImageConstPointer					 m_ReferenceImage;
  std::vector< std::vector<TransformPointerType> > m_Transforms;
  RegionType m_OutputImageRegion;
  unsigned int m_PSF;

};

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkLeastSquaresVnlCostFunction.txx"
#endif

#endif /* BTKLEASTSQUARESVNLCOSTFUNCTION_H_ */
