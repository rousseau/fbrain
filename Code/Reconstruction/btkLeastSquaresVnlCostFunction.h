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

  public: LeastSquaresVnlCostFunction(): vnl_cost_function(1) {}

  typedef vnl_vector<float> VnlVectorType;
  typedef vnl_sparse_matrix<float> VnlSparseMatrixType;

  double f(const vnl_vector<float>& x)
  {
    vnl_vector<float> HxMinusY;
    H.mult(x,HxMinusY);

    return HxMinusY.squared_magnitude();
  }

  // TODO Hvec should be pass as a const argument, see how to improve this
  // (an error is obtained with get_row function)
  void SetParameters(std::vector< VnlSparseMatrixType >& Hvec, const std::vector< VnlVectorType >& yvec)
  {

    vnl_sparse_matrix<float> Ht;

    // Precompute some matrices
    unsigned int ncols = Hvec[0].cols();
    unsigned int nrows = 0;

    for(unsigned int im = 0; im < Hvec.size(); im++)
      nrows += Hvec[im].rows();

    H.set_size(nrows,ncols);
    Ht.set_size(ncols,nrows);
    Y.set_size(nrows);

    // Precalcule H, Ht, and Y
    std::cout << "Precomputing H, Ht, and Y" << std::endl; std::cout.flush();

    unsigned int offset = 0;

    // Create H, Ht, and Y

    for(unsigned int im = 0; im < Hvec.size(); im++)
    {
      for (unsigned int i = 0; i < Hvec[im].rows(); i++)
      {
        vnl_sparse_matrix<float>::row & r = Hvec[im].get_row(i);
        vnl_sparse_matrix<float>::row::iterator col_iter = r.begin();

        for ( ;col_iter != r.end(); ++col_iter)
        {
          H( i + offset, (*col_iter).first ) = (*col_iter).second;
          Ht( (*col_iter).first, i + offset ) = (*col_iter).second;
        }

        Y[i+offset] = yvec[im][i];
      }
      offset += Hvec[im].rows();
    }

    // precalcule Ht * H
    std::cout << "Precomputing Ht*H" << std::endl; std::cout.flush();
    HtH = Ht * H;

    // precalcule Ht * Y
    std::cout << "Precomputing Ht*Y" << std::endl; std::cout.flush();
    Ht.mult(Y,HtY);
    Ht.clear();

  }

  void gradf(const vnl_vector<float>& x, vnl_vector<float>& g)
  {
    std::cout << "in gradf " << std::endl; std::cout.flush();

    VnlVectorType HtHx;

    HtH.mult(x,HtHx);
    g = (HtY + HtHx)*2.0;

    std::cout << "exiting of gradf " << std::endl; std::cout.flush();
  }
};

} // namespace btk

#endif /* BTKLEASTSQUARESVNLCOSTFUNCTION_H_ */
