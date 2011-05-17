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

namespace btk
{

class LeastSquaresVnlCostFunction : public vnl_cost_function
{
  private:

  vnl_sparse_matrix<float> H;
  vnl_sparse_matrix<float> HtH;
  vnl_vector<float> HtY;
  vnl_vector<float> Y;

  float lambda;

  struct img_size {
     unsigned int width;
     unsigned int height;
     unsigned int depth;
  } x_size;

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

  public: LeastSquaresVnlCostFunction(unsigned int dim): vnl_cost_function(dim) { lambda = 0.0; }

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
    std::cout << "In f(.)" << std::endl; std::cout.flush();

    vnl_vector<float>  x_float;
    x_float = vnl_matops::d2f(x);

    vnl_vector<float> Hx;
    H.mult(x_float,Hx);

    vnl_vector<float> HxMinusY;
    HxMinusY = Hx - Y;

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

    double mse = HxMinusY.squared_magnitude() / HxMinusY.size();
    double reg = lambda*(dX2dx + dX2dy + dX2dz) / x_float.size();

    std::cout << "mse, reg = " << mse << " , " << reg << std::endl;

    double value = mse + reg;

    std::cout << "error = " << value << std::endl;

    return value;

  }

  // TODO Hvec should be pass as a const argument, see how to improve this
  // (an error is obtained with get_row function)
  void SetParameters(const VnlSparseMatrixType & Hin, const VnlVectorType & Yin, const vnl_vector<double>& x, const vnl_vector<int>& x_size_in)
  {

    H = Hin;
    Y = Yin;

    x_size.width 	= x_size_in[0];
    x_size.height = x_size_in[1];
    x_size.depth  = x_size_in[2];

    vnl_sparse_matrix<float> Ht;

    unsigned int Hcols = H.cols();
    unsigned int Hrows = H.rows();

    Ht.set_size(Hcols,Hrows);

    // Calculate Ht
    std::cout << "Precomputing H, Ht, and Y" << std::endl; std::cout.flush();
    for( H.reset(); H.next(); )
      Ht( H.getcolumn(), H.getrow() ) = H.value();

    // precalcule Ht * H
    std::cout << "Precomputing Ht*H" << std::endl; std::cout.flush();
    HtH = Ht * H;

    // precalcule Ht * Y
    std::cout << "Precomputing Ht*Y" << std::endl; std::cout.flush();
    Ht.mult(Y,HtY);
    Ht.clear();

  }

  void gradf(const vnl_vector<double>& x, vnl_vector<double>& g)
  {
    vnl_vector<float>  x_float;
    vnl_vector<float>  g_float;

    x_float = vnl_matops::d2f(x);

    std::cout << "in gradf " << std::endl; std::cout.flush();

    VnlVectorType HtHx;

    HtH.mult(x_float,HtHx);
    g_float = (-HtY + HtHx)*2.0;

    // regularization term

    float* kernel = new float[3];
    kernel[0] = -1; kernel[1] = 2; kernel[2] = -1;

    vnl_vector<float> KxX;
    KxX.set_size( x_float.size() );
    convol3dx(x_float, KxX, x_size, kernel, 3);

    vnl_vector<float> KyX;
    KyX.set_size( x_float.size() );
    convol3dy(x_float, KyX, x_size, kernel, 3);

    vnl_vector<float> KzX;
    KzX.set_size( x_float.size() );
    convol3dz(x_float, KzX, x_size, kernel, 3);

    delete[] kernel;


    for (unsigned int i=0; i<g.size(); i++)
      g[i] = g_float[i] / Y.size();

    g_float = g_float + 2*lambda*(KxX + KyX + KzX)/ x_float.size();


    std::cout << "exiting of gradf " << std::endl; std::cout.flush();
  }
};

} // namespace btk

#endif /* BTKLEASTSQUARESVNLCOSTFUNCTION_H_ */
