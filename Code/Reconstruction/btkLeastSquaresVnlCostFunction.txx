/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 12/09/2011
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

#ifndef __LeastSquaresVnlCostFunction_txx
#define __LeastSquaresVnlCostFunction_txx

namespace btk
{

template <class TImage>
LeastSquaresVnlCostFunction<TImage>::LeastSquaresVnlCostFunction(unsigned int dim):
vnl_cost_function(dim)
{
  lambda = 0.1;
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::set(float * array, int size, float value)
{
  for (int i = 0; i < size; i++)
  array[i] = value;
}

// Mirror of the position pos. abs(pos) must not be > 2*(size-1)
template <class TImage>
int
LeastSquaresVnlCostFunction<TImage>::mirror(int pos, int size)
{
  int output = abs(pos);

  while (output >= size)
    output = abs(output - (output - size + 1) * 2);

  return output;
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::get_row(const vnl_vector<float>& image,
    img_size & size, int row, int frame, float * output)
{
  for (unsigned int i = 0; i < size.width; i++)
  output[i] = image[i + row * size.width + frame * size.width * size.height];
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::set_row(vnl_vector<float>& image,
    img_size & size, int row, int frame, float * input)
{
  for (unsigned int i = 0; i < size.width; i++)
    image[i + row * size.width + frame * size.width * size.height] = input[i];
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::get_col(const vnl_vector<float>& image,
    img_size & size, int col, int frame, float * output)
{
  for (unsigned int i = 0; i < size.height; i++)
    output[i] = image[col + i * size.width + frame * size.width * size.height];
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::set_col(vnl_vector<float>& image,
    img_size & size, int col, int frame, float * input)
{
  for (unsigned int i = 0; i < size.height; i++)
    image[col + i * size.width + frame * size.width * size.height] = input[i];
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::get_spec(const vnl_vector<float>& image,
    img_size & size, int row, int col, float * output)
{
  for (unsigned int i = 0; i < size.depth; i++)
  output[i] = image[col + row * size.width + i * size.width * size.height];
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::set_spec(vnl_vector<float>& image,
    img_size & size, int row, int col, float * input)
{
  for (unsigned int i = 0; i < size.depth; i++)
    image[col + row * size.width + i * size.width * size.height] = input[i];
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::convol1d(float * kernel, int ksize,
    float * src, int src_size, float * dest)
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
template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::convol3dx(const vnl_vector<float>& image,
    vnl_vector<float>& image_conv, img_size& size, float * kernel, int ksize)
{
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
template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::convol3dy(const vnl_vector<float>& image,
    vnl_vector<float>& image_conv, img_size & size, float * kernel, int ksize)
{
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
template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::convol3dz(const vnl_vector<float>& image,
    vnl_vector<float>& image_conv, img_size & size, float * kernel, int ksize)
{
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

// 3D convolution : over the spectra
template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::get_derivative_x(const vnl_vector<double>& x,
    vnl_vector<float>& deriv_x)
{

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

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::get_derivative_y(const vnl_vector<double>& x,
    vnl_vector<float>& deriv_y)
{

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

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::get_derivative_z(const vnl_vector<double>& x,
    vnl_vector<float>& deriv_z)
{

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

template <class TImage>
double
LeastSquaresVnlCostFunction<TImage>::f(const vnl_vector<double>& x)
{
  // Calculate the error with respect to the low resolution images

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
  double reg = 0.0;

  float* kernel = new float[2];
  kernel[0] = -1; kernel[1] = 1;

  vnl_vector<float> DxX;
  DxX.set_size( x_float.size() );
  convol3dx(x_float, DxX, x_size, kernel, 2);
  for(unsigned int i=0; i<x_float.size(); i++)
    reg += DxX[i]*DxX[i] / x_float.size();
  DxX.clear();

  vnl_vector<float> DyX;
  DyX.set_size( x_float.size() );
  convol3dy(x_float, DyX, x_size, kernel, 2);
  for(unsigned int i=0; i<x_float.size(); i++)
    reg += DyX[i]*DyX[i] / x_float.size();
  DyX.clear();

  vnl_vector<float> DzX;
  DzX.set_size( x_float.size() );
  convol3dz(x_float, DzX, x_size, kernel, 2);
  for(unsigned int i=0; i<x_float.size(); i++)
    reg += DzX[i]*DzX[i] / x_float.size();
  DzX.clear();

  delete[] kernel;

  // Calculate the cost function by combining both terms

  double value = mse + lambda*reg;

  std::cout << "error, mse, reg = " << value << " , " << mse << " , "
      << lambda*reg << std::endl;

  return value;

}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::gradf(const vnl_vector<double>& x,
    vnl_vector<double>& g)
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

/*    factor = 2*lambda/x_float.size();

  float* kernel = new float[3];
  kernel[0] = -1; kernel[1] = 2; kernel[2] = -1;

  vnl_vector<float> DxX;
  DxX.set_size( x_float.size() );
  convol3dx(x_float, DxX, x_size, kernel, 3);
  g = g + vnl_matops::f2d( DxX * factor );
  DxX.clear();

  vnl_vector<float> DyX;
  DyX.set_size( x_float.size() );
  convol3dy(x_float, DyX, x_size, kernel, 3);
  g = g + vnl_matops::f2d( DyX * factor );
  DyX.clear();

  vnl_vector<float> DzX;
  DzX.set_size( x_float.size() );
  convol3dz(x_float, DzX, x_size, kernel, 3);
  g = g + vnl_matops::f2d( DzX * factor );
  DzX.clear();

  delete[] kernel; */

}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::Initialize()
{
  m_OutputImageRegion = m_ReferenceImage -> GetLargestPossibleRegion();
  IndexType start_hr  = m_OutputImageRegion.GetIndex();
  SizeType  size_hr   = m_OutputImageRegion.GetSize();

  //x_size : size of the SR image (used in other functions)
  x_size.width  = size_hr[0];
  x_size.height = size_hr[1];
  x_size.depth  = size_hr[2];

  IndexType end_hr;
  end_hr[0] = start_hr[0] + size_hr[0] - 1 ;
  end_hr[1] = start_hr[1] + size_hr[1] - 1 ;
  end_hr[2] = start_hr[2] + size_hr[2] - 1 ;

  // Differential continuous indexes to perform the neighborhood iteration
  SpacingType spacing_lr = m_Images[0] -> GetSpacing();
  SpacingType spacing_hr = m_ReferenceImage -> GetSpacing();

  //spacing_lr[2] is assumed to be the lowest resolution 
  //compute the index of the PSF in the LR image resolution
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

  unsigned int im;
  #pragma omp parallel for private(im) schedule(dynamic)

  for(im = 0; im < m_Images.size(); im++)
  {

    // Interpolator for HR image
    InterpolatorPointer interpolator = InterpolatorType::New();
    interpolator -> SetInputImage( m_ReferenceImage );

    SpacingType inputSpacing = m_Images[im] -> GetSpacing();

    // PSF definition
    //TODO Give the possibility to choose the PDF
    typename FunctionType::Pointer function = FunctionType::New();
    function -> SetPSF( 1 );
    function -> SetDirection( m_Images[im] -> GetDirection() );
    function -> SetSpacing( m_Images[im] -> GetSpacing() );

    //Define the ROI of the current LR image
    IndexType inputIndex = m_Regions[im].GetIndex();
    SizeType  inputSize  = m_Regions[im].GetSize();

    //Define all indexes needed for iteration over the slices
    IndexType lrIndex;              //index of a point in the LR image im
    IndexType lrDiffIndex;          //index of this point in the current ROI of the LR image im
    unsigned int lrLinearIndex;     //index lineaire de ce point dans le vecteur 

    IndexType hrIndex;
    IndexType hrDiffIndex;
    ContinuousIndexType hrContIndex;
    unsigned int hrLinearIndex;

    ContinuousIndexType nbIndex;

    PointType lrPoint;            //PSF center in world coordinates (PointType = worldcoordinate for ITK)
    PointType nbPoint;            //PSF point in world coordinate
    PointType transformedPoint;   //after applying current transform (registration)

    unsigned int offset = 0;
    for(unsigned int im2 = 0; im2 < im; im2++)
      offset += m_Regions[im2].GetNumberOfPixels();

    // Iteration over slices
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
        //Current index in the LR image
        lrIndex = fixedIt.GetIndex();
        
        //World coordinates of lrIndex using the image header
        m_Images[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );

        if ( m_Masks.size() > 0)
          if ( ! m_Masks[im] -> IsInside(lrPoint) )
            continue;

        //Compute the coordinates in the SR using the estimated registration
        transformedPoint = m_Transforms[im][i] -> TransformPoint( lrPoint );

        //check if this point is in the SR image (m_ReferenceImage)
        if ( ! interpolator -> IsInsideBuffer( transformedPoint ) )
          continue;

        //From the LR image coordinates to the LR ROI coordinates
        lrDiffIndex[0] = lrIndex[0] - inputIndex[0];
        lrDiffIndex[1] = lrIndex[1] - inputIndex[1];
        lrDiffIndex[2] = lrIndex[2] - inputIndex[2];
        
        //Compute the corresponding linear index of lrDiffIndex
        lrLinearIndex = lrDiffIndex[0] + lrDiffIndex[1]*inputSize[0] +
                        lrDiffIndex[2]*inputSize[0]*inputSize[1];

        //Get the intensity value in the LR image
        Y[lrLinearIndex + offset] = fixedIt.Get();

        //Set the center point of the PSF
        function -> SetCenter( lrPoint );

        //Loop over points of the PSF
        for(unsigned int k=0; k<deltaIndexes.size(); k++)
        {
          //Coordinates in the LR image
          nbIndex[0] = deltaIndexes[k][0] + lrIndex[0];
          nbIndex[1] = deltaIndexes[k][1] + lrIndex[1];
          nbIndex[2] = deltaIndexes[k][2] + lrIndex[2];

          //World coordinates using LR image header 
          m_Images[im] -> TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );
          
          //Compute the PSF value at this point
          lrValue = function -> Evaluate(nbPoint);

          if ( lrValue > 0)
          {
            //Compute the world coordinate of this point in the SR image
            transformedPoint = m_Transforms[im][i] -> TransformPoint( nbPoint );

            //Set this coordinate in continuous index in SR image space
            m_ReferenceImage -> TransformPhysicalPointToContinuousIndex(
                                transformedPoint, hrContIndex );

            bool isInsideHR = true;

            // FIXME This checking should be done for all points first, and
            // discard the point if al least one point is out of the reference
            // image

            if ( (hrContIndex[0] < start_hr[0]) || (hrContIndex[0] > end_hr[0]) ||
                 (hrContIndex[1] < start_hr[1]) || (hrContIndex[1] > end_hr[1]) ||
                 (hrContIndex[2] < start_hr[2]) || (hrContIndex[2] > end_hr[2]) )
               isInsideHR = false;

            if ( isInsideHR )
            {
              //Compute the corresponding value in the SR image -> useless
              //Allows to compute the set of contributing neighbors
              hrValue = interpolator -> Evaluate( transformedPoint );

              //Loop over points affected using the interpolation
              for(unsigned int n=0; n<interpolator -> GetContributingNeighbors();
                  n++)
              {
                //Index in the SR image
                hrIndex = interpolator -> GetIndex(n);

                //Index in the ROI of the SR index
                hrDiffIndex[0] = hrIndex[0] - start_hr[0];
                hrDiffIndex[1] = hrIndex[1] - start_hr[1];
                hrDiffIndex[2] = hrIndex[2] - start_hr[2];

                //Compute the corresponding linear index
                hrLinearIndex = hrDiffIndex[0] + hrDiffIndex[1]*size_hr[0] +
                                hrDiffIndex[2]*size_hr[0]*size_hr[1];
                                
                //Add the correct value in H !                
                H(lrLinearIndex + offset, hrLinearIndex) +=
                    interpolator -> GetOverlap(n)* lrValue;
              }

            }

          }

        }

      }

    }

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

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::SetLambda(float value)
{
  lambda = value;
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::AddImage( ImageType* image )
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

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::AddRegion( RegionType region)
{
  m_Regions.push_back( region );
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::AddMask( MaskType *mask)
{
  m_Masks.push_back( mask );
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::SetReferenceImage( const ImageType * image )
{
  m_ReferenceImage = image;
}

template <class TImage>
void
LeastSquaresVnlCostFunction<TImage>::SetTransform( int i, int j, TransformType* transform )
{
  m_Transforms[i][j] = transform;
}

} // namespace btk

#endif /* LeastSquaresVnlCostFunction_txx */
