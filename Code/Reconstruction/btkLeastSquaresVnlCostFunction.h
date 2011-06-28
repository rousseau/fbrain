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
template <class TImage>
class LeastSquaresVnlCostFunction : public vnl_cost_function
{
  private:

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

  /**Blurring function typedef. */
  typedef OrientedSpatialFunction<double, 3, PointType> FunctionType;

  /**Const iterator typedef. */
  typedef ImageRegionConstIteratorWithIndex< ImageType >  ConstIteratorType;

  vnl_sparse_matrix<float> H;
  vnl_vector<float> HtY;
  vnl_vector<float> Y;

  float lambda;

  struct img_size {
     unsigned int width;
     unsigned int height;
     unsigned int depth;
  } x_size;

  std::vector<ImagePointer>  m_Images;
  std::vector<RegionType>    m_Regions;
  std::vector<MaskPointer>   m_Masks;
  ImageConstPointer					 m_ReferenceImage;
  std::vector< std::vector<TransformPointerType> > m_Transforms;
  RegionType m_OutputImageRegion;

  void set(float * array, int size, float value) {
    for (int i = 0; i < size; i++)
      array[i] = value;
  }

  // Mirror of the position pos. abs(pos) must not be > 2*(size-1)
  int mirror(int pos, int size) {
    int output = abs(pos);

    while (output >= size)
      output = abs(output - (output - size + 1) * 2);

    return output;
  }

  void get_row(const vnl_vector<float>& image, img_size & size, int row, int frame, float * output) {
    for (unsigned int i = 0; i < size.width; i++)
      output[i] = image[i + row * size.width + frame * size.width * size.height];
  }

  void set_row(vnl_vector<float>& image, img_size & size, int row, int frame, float * input) {
    for (unsigned int i = 0; i < size.width; i++)
      image[i + row * size.width + frame * size.width * size.height] = input[i];
  }

  void get_col(const vnl_vector<float>& image, img_size & size, int col, int frame, float * output) {
    for (unsigned int i = 0; i < size.height; i++)
      output[i] = image[col + i * size.width + frame * size.width * size.height];
  }

  void set_col(vnl_vector<float>& image, img_size & size, int col, int frame, float * input) {
    for (unsigned int i = 0; i < size.height; i++)
      image[col + i * size.width + frame * size.width * size.height] = input[i];
  }

  void get_spec(const vnl_vector<float>& image, img_size & size, int row, int col, float * output) {
    for (unsigned int i = 0; i < size.depth; i++)
      output[i] = image[col + row * size.width + i * size.width * size.height];
  }

  void set_spec(vnl_vector<float>& image, img_size & size, int row, int col, float * input) {
    for (unsigned int i = 0; i < size.depth; i++)
      image[col + row * size.width + i * size.width * size.height] = input[i];
  }


  void convol1d(float * kernel, int ksize, float * src, int src_size, float * dest)
  {
    int n2 = ksize / 2;
    int k;

    set(dest, src_size, 0);
    for (int i = 0; i < src_size; i++)
    {
      for (int j = 0; j < ksize; j++)
      {
        k = i - j + n2;
        k = mirror(k, src_size); // dcb|abcd    wxyz|yxw
        dest[i] += kernel[j] * src[k];
      }
    }
  }

  // 3D convolution : over the rows
  void convol3dx(const vnl_vector<float>& image, vnl_vector<float>& image_conv, img_size& size, float * kernel, int ksize) {
    float * drow = new float[size.width];
    float * srow = new float[size.width];

    for (unsigned int l = 0; l < size.depth; l++) {
      for (unsigned int py = 0; py < size.height; py++) {
        get_row(image, size, py, l, srow);
        convol1d(kernel, ksize, srow, size.width, drow);
        set_row(image_conv, size, py, l, drow);
      }
    }

    delete[] drow;
    delete[] srow;
  }

  // 3D convolution : over the columns
  void convol3dy(const vnl_vector<float>& image, vnl_vector<float>& image_conv, img_size & size, float * kernel, int ksize) {
    float * dcol = new float[size.height];
    float * scol = new float[size.height];

    for (unsigned int l = 0; l < size.depth; l++) {
      for (unsigned int px = 0; px < size.width; px++) {
        get_col(image, size, px, l, scol);
        convol1d(kernel, ksize, scol, size.height, dcol);
        set_col(image_conv, size, px, l, dcol);
      }
    }

    delete[] dcol;
    delete[] scol;
  }

  // 3D convolution : over the spectra
  void convol3dz(const vnl_vector<float>& image, vnl_vector<float>& image_conv, img_size & size, float * kernel, int ksize) {
    float * dspec = new float[size.depth];
    float * sspec = new float[size.depth];

    for (unsigned int py = 0; py < size.height; py++) {
      for (unsigned int px = 0; px < size.width; px++) {
        get_spec(image, size, py, px, sspec);
        convol1d(kernel, ksize, sspec, size.depth, dspec);
        set_spec(image_conv, size, py, px, dspec);
      }
    }

    delete[] dspec;
    delete[] sspec;
  }

  public: LeastSquaresVnlCostFunction(unsigned int dim): vnl_cost_function(dim) { lambda = 0.1; }

  // 3D convolution : over the spectra
  void get_derivative_x(const vnl_vector<double>& x, vnl_vector<float>& deriv_x) {

    vnl_vector<float>  x_float;
    x_float = vnl_matops::d2f(x);

    float* kernel = new float[2];
    kernel[0] = -1; kernel[1] = 1;

    float * drow = new float[x_size.width];
    float * srow = new float[x_size.width];

    for (unsigned int l = 0; l < x_size.depth; l++) {
      for (unsigned int py = 0; py < x_size.height; py++) {
        get_row(x_float, x_size, py, l, srow);
        convol1d(kernel, 2, srow, x_size.width, drow);
        set_row(deriv_x, x_size, py, l, drow);
      }
    }

    delete[] drow;
    delete[] srow;
  }

  // 3D convolution : over the spectra
  void get_derivative_y(const vnl_vector<double>& x, vnl_vector<float>& deriv_y) {

    vnl_vector<float>  x_float;
    x_float = vnl_matops::d2f(x);

    float* kernel = new float[2];
    kernel[0] = -1; kernel[1] = 1;

    float * dcol = new float[x_size.height];
    float * scol = new float[x_size.height];

    for (unsigned int l = 0; l < x_size.depth; l++) {
      for (unsigned int px = 0; px < x_size.width; px++) {
        get_col(x_float, x_size, px, l, scol);
        convol1d(kernel, 2, scol, x_size.height, dcol);
        set_col(deriv_y, x_size, px, l, dcol);
      }
    }

    delete[] dcol;
    delete[] scol;
  }

  // 3D convolution : over the spectra
  void get_derivative_z(const vnl_vector<double>& x, vnl_vector<float>& deriv_z) {

    vnl_vector<float>  x_float;
    x_float = vnl_matops::d2f(x);

    float* kernel = new float[2];
    kernel[0] = -1; kernel[1] = 1;

    float * dspec = new float[x_size.depth];
    float * sspec = new float[x_size.depth];

    for (unsigned int py = 0; py < x_size.height; py++) {
      for (unsigned int px = 0; px < x_size.width; px++) {
        get_spec(x_float, x_size, py, px, sspec);
        convol1d(kernel, 2, sspec, x_size.depth, dspec);
        set_spec(deriv_z, x_size, py, px, dspec);
      }
    }

    delete[] dspec;
    delete[] sspec;
  }

  typedef vnl_vector<float> VnlVectorType;
  typedef vnl_sparse_matrix<float> VnlSparseMatrixType;

  double f(const vnl_vector<double>& x)
  {
    vnl_vector<float>  x_float;
    x_float = vnl_matops::d2f(x);

    vnl_vector<float> Hx;
    H.mult(x_float,Hx);

    vnl_vector<float> HxMinusY;
    HxMinusY = Hx - Y;

    Hx.clear();

    double mse = HxMinusY.squared_magnitude() / HxMinusY.size();

    HxMinusY.clear();

    // Calculate the square of 1st derivatives along x, y, and z

    float* kernel = new float[3];
    kernel[0] = -1; kernel[1] = 2; kernel[2] = -1;

    vnl_vector<float> KxX;
    KxX.set_size( x_float.size() );
    convol3dx(x_float, KxX, x_size, kernel, 3);
    float dX2dx = dot_product(x_float, KxX);
    KxX.clear();

    vnl_vector<float> KyX;
    KyX.set_size( x_float.size() );
    convol3dy(x_float, KyX, x_size, kernel, 3);
    float dX2dy = dot_product(x_float, KyX);
    KyX.clear();

    vnl_vector<float> KzX;
    KzX.set_size( x_float.size() );
    convol3dz(x_float, KzX, x_size, kernel, 3);
    float dX2dz = dot_product(x_float, KzX);
    KzX.clear();

    delete[] kernel;

    // Calculate the cost function by combining both terms

    double reg = lambda*(dX2dx + dX2dy + dX2dz) / x_float.size();
    double value = mse + reg;

    std::cout << "error, mse, reg = " << value << " , " << mse << " , " << reg << std::endl;

    return value;

  }

  void gradf(const vnl_vector<double>& x, vnl_vector<double>& g)
  {
    vnl_vector<float>  x_float;

    x_float = vnl_matops::d2f(x);

    vnl_vector<float> Hx;
    H.mult(x_float,Hx);

    // Calculate Ht*Hx. Note that this is calculated as Hx*H since
    // Ht*Hx = (Hxt*H)t and for Vnl (Hxt*H)t = Hxt*H = Hx*H because
    // the vnl_vector doesn't have a 2nd dimension. This allows us
    // to save a lot of memory because we don't need to store Ht.

    vnl_vector<float> HtHx;
    H.pre_mult(Hx,HtHx);
    Hx.clear();

    double factor = 2.0 / Y.size();
    g = vnl_matops::f2d( (-HtY + HtHx)*factor );

    HtHx.clear();

    // regularization term

    factor = 2*lambda/x_float.size();

    float* kernel = new float[3];
    kernel[0] = -1; kernel[1] = 2; kernel[2] = -1;

    vnl_vector<float> KxX;
    KxX.set_size( x_float.size() );
    convol3dx(x_float, KxX, x_size, kernel, 3);

    g = g + vnl_matops::f2d( KxX * factor );
    KxX.clear();

    vnl_vector<float> KyX;
    KyX.set_size( x_float.size() );
    convol3dy(x_float, KyX, x_size, kernel, 3);

    g = g + vnl_matops::f2d( KyX * factor );
    KyX.clear();

    vnl_vector<float> KzX;
    KzX.set_size( x_float.size() );
    convol3dz(x_float, KzX, x_size, kernel, 3);

    g = g + vnl_matops::f2d( KzX * factor );
    KzX.clear();

    delete[] kernel;

  }

  void Initialize()
  {
    IndexType start_hr;
    SizeType  size_hr;

    m_OutputImageRegion = m_ReferenceImage -> GetLargestPossibleRegion();
    start_hr = m_OutputImageRegion.GetIndex();
    size_hr = m_OutputImageRegion.GetSize();

    x_size.width  = size_hr[0];
    x_size.height = size_hr[1];
    x_size.depth  = size_hr[2];

    IndexType end_hr;
    end_hr[0] = start_hr[0] + size_hr[0] - 1 ;
    end_hr[1] = start_hr[1] + size_hr[1] - 1 ;
    end_hr[2] = start_hr[2] + size_hr[2] - 1 ;

    // Interpolator for HR image

    InterpolatorPointer interpolator = InterpolatorType::New();
    interpolator -> SetInputImage( m_ReferenceImage );

    // Differential continuous indexes to perform the neighborhood iteration
    SpacingType spacing_lr = m_Images[0] -> GetSpacing();
    SpacingType spacing_hr = m_ReferenceImage -> GetSpacing();

    std::vector<ContinuousIndexType> deltaIndexes;
    int npoints =  spacing_lr[2] / (2.0 * spacing_hr[2]) ;
    ContinuousIndexType delta;
    delta[0] = 0.0; delta[1] = 0.0;

    for (int i = -npoints ; i <= npoints; i++ )
    {
      delta[2] = i * 0.5 / (double) npoints;
      deltaIndexes.push_back(delta);
    }

    // Set size of matrices
    unsigned int ncols = m_OutputImageRegion.GetNumberOfPixels();

    unsigned int nrows = 0;
    for(unsigned int im = 0; im < m_Images.size(); im++)
      nrows += m_Regions[im].GetNumberOfPixels();

    H.set_size(nrows, ncols);
    Y.set_size(nrows);
    Y.fill(0.0);

    unsigned int offset = 0;

    for(unsigned int im = 0; im < m_Images.size(); im++)
    {
      SpacingType inputSpacing = m_Images[im] -> GetSpacing();

      // PSF definition
      //TODO Give the possibility to choose the PDF
      typename FunctionType::Pointer function = FunctionType::New();
      function -> SetPSF( 1 );
      function -> SetDirection( m_Images[im] -> GetDirection() );
      function -> SetSpacing( m_Images[im] -> GetSpacing() );

      // Iteration over slices

      IndexType inputIndex = m_Regions[im].GetIndex();
      SizeType  inputSize  = m_Regions[im].GetSize();

      IndexType lrIndex;
      IndexType lrDiffIndex;
      unsigned int lrLinearIndex;

      IndexType hrIndex;
      IndexType hrDiffIndex;
      ContinuousIndexType hrContIndex;
      unsigned int hrLinearIndex;

      ContinuousIndexType nbIndex;

      PointType lrPoint;
      PointType nbPoint;
      PointType transformedPoint;

      for ( unsigned int i=inputIndex[2]; i < inputIndex[2] + inputSize[2]; i++ )
      {

        RegionType wholeSliceRegion;
        wholeSliceRegion = m_Regions[im];

        IndexType  wholeSliceRegionIndex = wholeSliceRegion.GetIndex();
        SizeType   wholeSliceRegionSize  = wholeSliceRegion.GetSize();

        wholeSliceRegionIndex[2]= i;
        wholeSliceRegionSize[2] = 1;

        wholeSliceRegion.SetIndex(wholeSliceRegionIndex);
        wholeSliceRegion.SetSize(wholeSliceRegionSize);

        ConstIteratorType fixedIt( m_Images[im], wholeSliceRegion);

        double lrValue;
        double hrValue;

        for(fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt)
        {
          lrIndex = fixedIt.GetIndex();
          m_Images[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );

          if ( m_Masks.size() > 0)
            if ( ! m_Masks[im] -> IsInside(lrPoint) )
              continue;


          if ( interpolator -> IsInsideBuffer( lrPoint ) )
          {

            lrDiffIndex[0] = lrIndex[0] - inputIndex[0];
            lrDiffIndex[1] = lrIndex[1] - inputIndex[1];
            lrDiffIndex[2] = lrIndex[2] - inputIndex[2];

            lrLinearIndex = lrDiffIndex[0] + lrDiffIndex[1]*inputSize[0] + lrDiffIndex[2]*inputSize[0]*inputSize[1];

            Y[lrLinearIndex + offset] = fixedIt.Get();

            m_Images[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );
            function -> SetCenter( lrPoint );

            for(unsigned int k=0; k<deltaIndexes.size(); k++)
            {
              nbIndex[0] = deltaIndexes[k][0] + lrIndex[0];
              nbIndex[1] = deltaIndexes[k][1] + lrIndex[1];
              nbIndex[2] = deltaIndexes[k][2] + lrIndex[2];

              m_Images[im] -> TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );
              lrValue = function -> Evaluate(nbPoint);

              if ( lrValue > 0)
              {
                transformedPoint = m_Transforms[im][i] -> TransformPoint( nbPoint);

                m_ReferenceImage -> TransformPhysicalPointToContinuousIndex( transformedPoint, hrContIndex );

                bool isInsideHR = true;

                // FIXME This checking should be done for all points first, and discard the point
                // if al least one point is out of the reference image

                if ( (hrContIndex[0] < start_hr[0]) || (hrContIndex[0] > end_hr[0]) ||
                     (hrContIndex[1] < start_hr[1]) || (hrContIndex[1] > end_hr[1]) ||
                     (hrContIndex[2] < start_hr[2]) || (hrContIndex[2] > end_hr[2]) )
                   isInsideHR = false;

                if ( isInsideHR )
                {
                  hrValue = interpolator -> Evaluate( transformedPoint );

                  for(unsigned int n=0; n<interpolator -> GetContributingNeighbors(); n++)
                  {
                    hrIndex = interpolator -> GetIndex(n);

                    hrDiffIndex[0] = hrIndex[0] - start_hr[0];
                    hrDiffIndex[1] = hrIndex[1] - start_hr[1];
                    hrDiffIndex[2] = hrIndex[2] - start_hr[2];

                    hrLinearIndex = hrDiffIndex[0] + hrDiffIndex[1]*size_hr[0] + hrDiffIndex[2]*size_hr[0]*size_hr[1];
                    H(lrLinearIndex + offset, hrLinearIndex) += interpolator -> GetOverlap(n)* lrValue;


                  }

                }

              }

            }

          }

        }

      }

      offset += m_Regions[im].GetNumberOfPixels();

    }

    // Normalize H
    for (unsigned int i = 0; i < H.rows(); i++)
    {
      double sum = H.sum_row(i);

      VnlSparseMatrixType::row & r = H.get_row(i);
      VnlSparseMatrixType::row::iterator col_iter = r.begin();

      for ( ;col_iter != r.end(); ++col_iter)
        (*col_iter).second = (*col_iter).second / sum;

    }

    // Precalcule Ht*Y. Note that this is calculated as Y*H since
    // Ht*Y = (Yt*H)t and for Vnl (Yt*H)t = (Yt*H) = Y*H because
    // the vnl_vector doesn't have a 2nd dimension. This allows us
    // to save a lot of memory because we don't need to store Ht.
    H.pre_mult(Y,HtY);

  }

  void SetLambda(float value)
  {
    lambda = value;
  }

  void AddImage( ImageType* image )
  {
    m_Images.push_back( image );

    // Add transforms for this image
    m_Transforms.resize( m_Transforms.size() + 1 );
    SizeType imageSize = image -> GetLargestPossibleRegion().GetSize();
    m_Transforms[m_Transforms.size()-1].resize( imageSize[2]);

    // Initialize transforms
    for (unsigned int i=0; i<imageSize[2]; i++)
      m_Transforms[m_Transforms.size()-1][i] = TransformType::New();

  }

  void AddRegion( RegionType region)
  {
    m_Regions.push_back( region );
  }

  void AddMask( MaskType *mask)
  {
    m_Masks.push_back( mask );
  }

  void SetReferenceImage( const ImageType * image )
  {
    m_ReferenceImage = image;
  }

  void SetTransform( int i, int j, TransformType* transform )
  {
    m_Transforms[i][j] = transform;
  }


};

} // namespace btk

#endif /* BTKLEASTSQUARESVNLCOSTFUNCTION_H_ */
