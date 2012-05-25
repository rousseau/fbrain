#include "btkLUDecomposition.h"


namespace btk
{
//----------------------------------------------------------------------------------------
LUDecomposition::LUDecomposition(BtkMatrix<double> &a) : n(a.nrows()), lu(a), aref(a), indx(n)
{
    const Doub TINY=1.0e-40;
    int i,imax,j,k;
    double big,temp;
    BtkVector<double> vv(n);
    d=1.0;

    for (i=0;i<n;i++)
    {
        big=0.0;

        for (j=0;j<n;j++)
            if ((temp=abs(lu[i][j])) > big)
                big=temp;

        if (big == 0.0)
            throw("Singular matrix in LUdcmp");

        vv[i]=1.0/big;
    }

    for (k=0;k<n;k++)
    {
        big=0.0;
        for (i=k;i<n;i++)
        {
            temp=vv[i]*abs(lu[i][k]);

            if (temp > big)
            {
                big=temp;
                imax=i;
            }
        }

        if (k != imax)
        {
            for (j=0;j<n;j++)
            {
                temp=lu[imax][j];
                lu[imax][j]=lu[k][j];
                lu[k][j]=temp;
            }

            d = -d;
            vv[imax]=vv[k];
        }

        indx[k]=imax;

        if (lu[k][k] == 0.0)
            lu[k][k]=TINY;

        for (i=k+1;i<n;i++)
        {
            temp=lu[i][k] /= lu[k][k];

            for (j=k+1;j<n;j++)
                lu[i][j] -= temp*lu[k][j];
        }
    }
}
//----------------------------------------------------------------------------------------
void LUDecomposition::solve(BtkVector<double> &b, BtkVector<double> &x)
{
    Int i,ii=0,ip,j;
    Doub sum;
    if (b.size() != n || x.size() != n)
        throw("LUdcmp::solve bad sizes");

    for (i=0;i<n;i++)
        x[i] = b[i];

    for (i=0;i<n;i++)
    {
        ip=indx[i];
        sum=x[ip];
        x[ip]=x[i];

        if (ii != 0)
            for (j=ii-1;j<i;j++)
                sum -= lu[i][j]*x[j];

        else if (sum != 0.0)
            ii=i+1;

        x[i]=sum;
    }
    for (i=n-1;i>=0;i--)
    {
        sum=x[i];

        for (j=i+1;j<n;j++)
            sum -= lu[i][j]*x[j];

        x[i]=sum/lu[i][i];
    }
}
//----------------------------------------------------------------------------------------
void LUDecomposition::solve(BtkMatrix<double> &b, BtkMatrix<double> &x)
{
    int i,j,m=b.ncols();
    if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
        throw("LUdcmp::solve bad sizes");
    BtkVector<double> xx(n);

    for (j=0;j<m;j++)
    {
        for (i=0;i<n;i++)
            xx[i] = b[i][j];

        solve(xx,xx);

        for (i=0;i<n;i++)
            x[i][j] = xx[i];
    }
}
//----------------------------------------------------------------------------------------
void LUDecomposition::inverse(BtkMatrix<double> &ainv)
{
    int i,j;
    ainv.resize(n,n);

    for (i=0;i<n;i++)
    {
        for (j=0;j<n;j++)
            ainv[i][j] = 0.;

        ainv[i][i] = 1.;
    }

    solve(ainv,ainv);
}
//----------------------------------------------------------------------------------------
Doub LUDecomposition::det()
{
    Doub dd = d;

    for (Int i=0;i<n;i++)
        dd *= lu[i][i];

    return dd;
}
//----------------------------------------------------------------------------------------
void LUDecomposition::mprove(BtkVector<double> &b, BtkVector<double> &x)
{
    int i,j;
    BtkVector<double> r(n);

    for (i=0;i<n;i++)
    {
        long double sdp = -b[i];

        for (j=0;j<n;j++)
            sdp += (long double)aref[i][j] * (long double)x[j];

        r[i]=sdp;
    }
    solve(r,r);

    for (i=0;i<n;i++)
        x[i] -= r[i];
}
//----------------------------------------------------------------------------------------
}


