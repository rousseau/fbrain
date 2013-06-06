/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/05/2012
  Author(s): Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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
#ifndef __BTKRBFINTERPOLATION_H__
#define __BTKRBFINTERPOLATION_H__

#include "btkNumerical.h"
#include "btkLUDecomposition.h"
#include "VNLSparseLUSolverTraits.h"

namespace btk
{
//--------------------------------------------------------------------------------------------------
struct RBF_fn
{
    virtual double rbf(double r) = 0;
};
//--------------------------------------------------------------------------------------------------
struct RBF2_fn {
  virtual double rbf(const double *p1, const double *p2) = 0;
};
//--------------------------------------------------------------------------------------------------
struct RBF_interp
{
    int dim, n;
    const BtkMatrix<double> &pts;
    const BtkVector<double> &vals;
    BtkVector<double> w;
    RBF_fn &fn;
    bool norm;

    //*******************************************************************************************
    RBF_interp(BtkMatrix<double> &ptss, BtkVector<double> &valss, RBF_fn &func, Bool nrbf=false)
    : dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss), vals(valss),
    w(n), fn(func), norm(nrbf)
    {
        int i,j;
        double sum;
        BtkMatrix<double> rbf(n,n);
        BtkVector<double> rhs(n);

        for (i=0;i<n;i++)
        {
            sum = 0.;
            for (j=0;j<n;j++)
            {
                sum += (rbf[i][j] = fn.rbf(rad(&pts[i][0],&pts[j][0])));
            }

            if (norm)
                rhs[i] = sum*vals[i];

            else
                rhs[i] = vals[i];
        }

        LUDecomposition lu(rbf);
        lu.solve(rhs,w);
    }
    //*******************************************************************************************
    double interp(BtkVector<double> &pt)
    {
        double fval, sum=0., sumw=0.;

        if (pt.size() != dim)
            throw("RBF_interp bad pt size");

        for (Int i=0;i<n;i++)
        {
            fval = fn.rbf(rad(&pt[0],&pts[i][0]));
            sumw += w[i]*fval;
            sum += fval;
        }
        return norm ? sumw/sum : sumw;
    }
    //*******************************************************************************************
    double rad(const double *p1, const double *p2)
    {
        double sum = 0.;

        for (Int i=0;i<dim;i++)
            sum += SQR(p1[i]-p2[i]);

        return sqrt(sum);
    }
    //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
struct RBF2_interp
{
    Int dim, n;
    const BtkMatrix<double> &pts;
    const BtkVector<double> &vals;
    //vnl_vector< double > w;
    BtkVector<double> w;
    RBF2_fn &fn;
    Bool norm;

  RBF2_interp(BtkMatrix<double> &ptss, BtkVector<double> &valss, RBF2_fn &func, Bool nrbf=false)
  : dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss), vals(valss),
  w(n), fn(func), norm(nrbf)
  {
    int i,j;
    double sum;
//    vnl_sparse_matrix< double > rbf(n,n);
//    vnl_vector< double > rhs(n);
    BtkMatrix<double> rbf(n,n);
    BtkVector<double> rhs(n);

    for (i=0;i<n;i++)
    {
      sum = 0.;
      for (j=0;j<n;j++)
      {
        sum += (rbf[i][j] = fn.rbf(&pts[i][0],&pts[j][0]));
        //sum += (rbf(i,j) = fn.rbf(&pts[i][0],&pts[j][0]));



      }

      if (norm)
      {
          rhs[i] = sum*vals[i];
          //rhs(i) = sum*vals[i];
      }

      else
      {
          rhs[i] = vals[i];
          //rhs(i) = vals[i];
      }

    }

//    typedef VNLSparseLUSolverTraits<double> LUSolverType;
//    LUSolverType * LUSolver  = new LUSolverType;
//    LUSolver -> Solve(rbf,rhs,w);
//    delete LUSolver;

    LUDecomposition lu(rbf);
    lu.solve(rhs,w);
  }


  double interp(BtkVector<double> &pt)
  {
    double fval, sum=0., sumw=0.;

    if (pt.size() != dim)
        throw("RBF_interp bad pt size");

    for (int i=0;i<n;i++)
    {
      fval = fn.rbf(&pt[0],&pts[i][0]);
      sumw += w[i]*fval;
      // sumw += w(i)*fval;
      sum += fval;

    }
//     std::cout<<"sum"<<sum<<std::endl;
//     std::cout<<"sumw"<<sumw<<std::endl;
    return norm ? sumw/sum : sumw;
  }

    double rad(const double *p1, const double *p2)
    {
        double sum = 0.;

        for (int i=0;i<dim;i++)
            sum += SQR(p1[i]-p2[i]);

        return sqrt(sum);

    }
};
//--------------------------------------------------------------------------------------------------
struct RBF_multiquadric : RBF_fn
{
    double r02;
    //*******************************************************************************************
    RBF_multiquadric(double scale=1.) : r02(SQR(scale)) {}
    //*******************************************************************************************
    double rbf(double r)
    {
        return sqrt(SQR(r)+r02);
    }
    //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
struct RBF_thinplate : RBF_fn
{
    double r0;
    //*******************************************************************************************
    RBF_thinplate(double scale=1.) : r0(scale) {}
    //*******************************************************************************************
    double rbf(double r)
    {
        return r <= 0. ? 0. : SQR(r)*log(r/r0);
    }
    //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
struct RBF_gauss : RBF_fn
{
    double r0;
    //*******************************************************************************************
    RBF_gauss(double scale=1.) : r0(scale) {}
    //*******************************************************************************************
    double rbf(double r)
    {
        return exp(-0.5*SQR(r/r0));
    }
    //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
struct RBF2_gauss : RBF2_fn
{
  Doub r0_spa;
  Doub r0_ang;
  //*******************************************************************************************
  RBF2_gauss(Doub scale1=1, Doub scale2=1) : r0_spa(scale1),r0_ang(scale2) {}
  //*******************************************************************************************
  Doub rbf(const Doub *p1, const Doub *p2)
  {
    double r_spa = 0;
    double r_ang = 0;

    double q1[3];
    double q2[3];

    for (Int i=0;i<3;i++)
        r_spa += SQR(p1[i]-p2[i]);


    r_spa = sqrt( r_spa );

    q1[0] = std::sin(p1[3]) * std::cos(p1[4]);//m_GradientTableCartesian[1][0];
    q1[1] = std::sin(p1[3]) * std::sin(p1[4]);
    q1[2] = std::cos(p1[3]);

    q2[0] = std::sin(p2[3]) * std::cos(p2[4]);//m_GradientTableCartesian[1][0];
    q2[1] = std::sin(p2[3]) * std::sin(p2[4]);
    q2[2] = std::cos(p2[3]);

    for (int i=0;i<3;i++)
        r_ang += q1[i]*q2[i];

    // value can be higher than 1 at machine precision
    if ( abs(r_ang) > 1  )
    {
//      std::cout << " entro !!! : " << r_ang << " " << abs(r_ang) << " " << (abs(r_ang) - 1) << std::endl;
      r_ang = 1;
    }
    r_ang = std::acos(abs( r_ang ));
    //std::cout<<r_spa<<" - "<<r_ang<<std::endl;
    double value = exp(-0.5*SQR(r_spa/r0_spa))*exp(-0.5*SQR(r_ang/r0_ang));
    return value;
  }
  //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
struct RBF_inversemultiquadric : RBF_fn
{
    double r02;
    //*******************************************************************************************
    RBF_inversemultiquadric(double scale=1.) : r02(SQR(scale)) {}
    //*******************************************************************************************
    double rbf(double r)
    {
        return 1./sqrt(SQR(r)+r02);
    }
    //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
struct Shep_interp
{
    Int dim, n;
    const BtkMatrix<double> &pts;
    const BtkVector<double> &vals;
    double pneg;
    //*******************************************************************************************
    Shep_interp(BtkMatrix<double> &ptss, BtkVector<double> &valss, double p=2.)
    : dim(ptss.ncols()), n(ptss.nrows()) , pts(ptss),
    vals(valss), pneg(-p) {}
    //*******************************************************************************************
    double interp(BtkVector<double> &pt)
    {
        double r, w, sum=0., sumw=0.;

        if (pt.size() != dim)
            throw("RBF_interp bad pt size");

        for (Int i=0;i<n;i++)
        {
            if ((r=rad(&pt[0],&pts[i][0])) == 0.)
                return vals[i];

            sum += (w = pow(r,pneg));
            sumw += w*vals[i];
        }
        return sumw/sum;
    }
    //*******************************************************************************************
    double rad(const double *p1, const double *p2)
    {
        double sum = 0.;

        for (Int i=0;i<dim;i++)
            sum += SQR(p1[i]-p2[i]);

        return sqrt(sum);
    }
    //*******************************************************************************************
};
//--------------------------------------------------------------------------------------------------
}



#endif
