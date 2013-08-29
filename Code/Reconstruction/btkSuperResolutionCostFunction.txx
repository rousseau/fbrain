/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 29/05/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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

#include "btkSuperResolutionCostFunction.hxx"

namespace btk
{
template< class TImage >
SuperResolutionCostFunction< TImage >::SuperResolutionCostFunction()
{
}
//-------------------------------------------------------------------------------------------------
template < class TImage >
double SuperResolutionCostFunction< TImage >::
GetValue(const vnl_vector<double> &_x)
{
    // Calculate the error with respect to the low resolution images

    vnl_vector<PrecisionType>  x_float;
    x_float = vnl_matops::d2f(_x);

    vnl_vector<PrecisionType> Hx;
    this->m_H.mult(x_float,Hx);

    m_X_Size.width = m_SRSize[0];
    m_X_Size.height = m_SRSize[1];
    m_X_Size.depth = m_SRSize[2];

//    vnl_vector<PrecisionType> Hx;
//    this->m_H.mult(_x,Hx);

    //vnl_vector<float> HxMinusY;
    vnl_vector<PrecisionType> HxMinusY;
    HxMinusY = Hx - this->m_Y;
    //HxMinusY = Y - Hx;

    Hx.clear();

    double mse = HxMinusY.squared_magnitude()/ HxMinusY.size();
    //double mse = HxMinusY.magnitude()/ HxMinusY.size();
    //double mse = HxMinusY.two_norm()/ HxMinusY.size();
    //double mse = HxMinusY.one_norm()/ HxMinusY.size();
    //double mse = HxMinusY.rms()/ HxMinusY.size();

    HxMinusY.clear();

    // Calculate the square of 1st derivatives along x, y, and z
    double reg = 0.0;
    double regCH = 0.0;

    //float* kernel = new float[2];
    PrecisionType* kernel = new PrecisionType[2];
    kernel[0] = -1; kernel[1] = 1;


    vnl_vector<float> DxX;
    DxX.set_size( x_float.size() );

/*    vnl_vector<PrecisionType> DxX;
    DxX.set_size( _x.size() )*/;

    Convol3dx(x_float, DxX, m_X_Size, kernel, 2);
    //Convol3dx(_x, DxX, m_X_Size, kernel, 2);

       //for(unsigned int i=0; i<x_float.size(); i++)
    for(unsigned int i=0; i<_x.size(); i++)
    {
        //reg += DxX[i]*DxX[i] / x_float.size();
        reg += DxX[i]*DxX[i] / _x.size();
        regCH += 2* sqrt(1 + (DxX[i]*DxX[i]/ _x.size()))-2;
    }
    DxX.clear();

    vnl_vector<float> DyX;
    DyX.set_size( x_float.size() );

//    vnl_vector<PrecisionType> DyX;
//    DyX.set_size( _x.size() );

    Convol3dy(x_float, DyX, m_X_Size, kernel, 2);
    //Convol3dy(_x, DyX, m_X_Size, kernel, 2);
    //for(unsigned int i=0; i<x_float.size(); i++)
    for(unsigned int i=0; i<_x.size(); i++)
    {
        //reg += DyX[i]*DyX[i] / x_float.size();
        reg += DyX[i]*DyX[i] / _x.size();
        regCH += 2 * sqrt(1 + (DyX[i]*DyX[i]/ _x.size())) -2;
    }
    DyX.clear();

    vnl_vector<float> DzX;
    DzX.set_size( x_float.size() );

//    vnl_vector<PrecisionType> DzX;
//    DzX.set_size( _x.size() );

    Convol3dz(x_float, DzX, m_X_Size, kernel, 2);
    //Convol3dz(_x, DzX, m_X_Size, kernel, 2);
    //for(unsigned int i=0; i<x_float.size(); i++)
    for(unsigned int i=0; i<_x.size(); i++)
    {
        //reg += DzX[i]*DzX[i] / x_float.size();
        reg += DzX[i]*DzX[i] / _x.size();
        regCH += 2 * sqrt(1 + (DzX[i]*DzX[i]/ _x.size())) -2;
    }
    DzX.clear();

    delete[] kernel;

    // Calculate the cost function by combining both terms
    //std::cout<<"Reg : "<<reg<<", CH reg : "<<regCH<<std::endl;

    double value = mse + m_Lambda*regCH;

    //std::cout << "error, mse, reg = " << value << " , " << mse << " , " <<reg<<" , "<< m_Lambda*reg << std::endl;

    return value;
}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void SuperResolutionCostFunction< TImage >::
GetGradient(const vnl_vector< double > & _x, vnl_vector< double >& _g)
{
    vnl_vector<float>  x_float;

    x_float = vnl_matops::d2f(_x);

    vnl_vector<float> Hx;
    m_H.mult(x_float,Hx);

//    vnl_vector< PrecisionType > Hx;
//    m_H.mult(_x,Hx);


    // Calculate Ht*Hx. Note that this is calculated as Hx*H since
    // Ht*Hx = (Hxt*H)t and for Vnl (Hxt*H)t = Hxt*H = Hx*H because
    // the vnl_vector doesn't have a 2nd dimension. This allows us
    // to save a lot of memory because we don't need to store Ht.

    //vnl_vector<float> HtHx;
    vnl_vector<PrecisionType> HtHx;
    m_H.pre_mult(Hx,HtHx);
    Hx.clear();

    double factor = 2.0 / m_Y.size();
    _g = vnl_matops::f2d( (-m_HtY + HtHx)*factor );
    //_g = (-m_HtY + HtHx)*factor ;



    HtHx.clear();
}

//-------------------------------------------------------------------------------------------------
template < class TImage >
int SuperResolutionCostFunction< TImage >::Mirror(int _pos, int _size)
{
    int output = std::abs(_pos);

    while(output >= _size)
    {
        output = std::abs(output - (output - _size + 1) * 2);
    }

    return output;
}

//-------------------------------------------------------------------------------------------------
template < class TImage >
void SuperResolutionCostFunction< TImage >::
Convol1d(PrecisionType *_kernel, int _ksize, PrecisionType *_src, int _src_size, PrecisionType *_dest)
{
    int n2 = _ksize / 2;
    int k;

    Set(_dest, _src_size, 0);
    for (int i = 0; i < _src_size; i++)
    {
        for (int j = 0; j < _ksize; j++)
        {
            k = i - j + n2;
            k = this->Mirror(k, _src_size);
            _dest[i] += _kernel[j] * _src[k];
        }
    }
}

//-------------------------------------------------------------------------------------------------
template < class TImage >
void SuperResolutionCostFunction< TImage >::
Convol3dx(const vnl_vector<PrecisionType> &_image, vnl_vector<PrecisionType> &_image_conv, IMG_SIZE &_size, PrecisionType *_kernel, int _ksize)
{
    PrecisionType * drow = new PrecisionType[_size.width];
    PrecisionType * srow = new PrecisionType[_size.width];

    for (unsigned int l = 0; l < _size.depth; l++)
    {
        for (unsigned int py = 0; py < _size.height; py++)
        {
            GetRow(_image, _size, py, l, srow);
            Convol1d(_kernel, _ksize, srow, _size.width, drow);
            SetRow(_image_conv, _size, py, l, drow);
        }
    }

    delete[] drow;
    delete[] srow;
}

//-------------------------------------------------------------------------------------------------
template < class TImage >
void SuperResolutionCostFunction< TImage >::
Convol3dy(const vnl_vector<PrecisionType> &_image, vnl_vector<PrecisionType> &_image_conv, IMG_SIZE &_size, PrecisionType *_kernel, int _ksize)
{
    PrecisionType * dcol = new PrecisionType[_size.height];
    PrecisionType * scol = new PrecisionType[_size.height];

    for (unsigned int l = 0; l < _size.depth; l++)
    {
        for (unsigned int px = 0; px < _size.width; px++)
        {
            GetCol(_image, _size, px, l, scol);
            Convol1d(_kernel, _ksize, scol, _size.height, dcol);
            SetCol(_image_conv, _size, px, l, dcol);
        }
    }

    delete[] dcol;
    delete[] scol;
}

//-------------------------------------------------------------------------------------------------
template < class TImage >
void SuperResolutionCostFunction< TImage >::
Convol3dz(const vnl_vector<PrecisionType> &_image, vnl_vector<PrecisionType> &_image_conv, IMG_SIZE &_size, PrecisionType *_kernel, int _ksize)
{
    PrecisionType * dspec = new PrecisionType[_size.depth];
    PrecisionType * sspec = new PrecisionType[_size.depth];

    for (unsigned int py = 0; py < _size.height; py++)
    {
        for (unsigned int px = 0; px < _size.width; px++)
        {
            GetSpec(_image, _size, py, px, sspec);
            Convol1d(_kernel, _ksize, sspec, _size.depth, dspec);
            SetSpec(_image_conv, _size, py, px, dspec);
        }
    }

    delete[] dspec;
    delete[] sspec;
}
//-------------------------------------------------------------------------------------------------
}
