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
#ifndef __BTKNUMERICAL_H__
#define __BTKNUMERICAL_H__

#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>


namespace btk
{
using namespace std;

//----------------------------------------------------------------------------------------------
template<class T>
inline T SQR(const T _a)
{
    return _a*_a;
}
//----------------------------------------------------------------------------------------------
template<class T>
inline const T &MAX(const T &_a, const T &_b)
{
    return _b > _a ? (_b) : (_a);
}
//----------------------------------------------------------------------------------------------
inline float MAX(const double &_a, const float &_b)
{
    return _b > _a ? (_b) : float(_a);
}
//----------------------------------------------------------------------------------------------
inline float MAX(const float &_a, const double &_b)
{
    return _b > _a ? float(_b) : (_a);
}
//----------------------------------------------------------------------------------------------
template<class T>
inline const T &MIN(const T &_a, const T &_b)
{
    return _b < _a ? (_b) : (_a);
}
//----------------------------------------------------------------------------------------------
inline float MIN(const double &_a, const float &_b)
{
    return _b < _a ? (_b) : float(_a);
}
//----------------------------------------------------------------------------------------------
inline float MIN(const float &_a, const double &_b)
{
    return _b < _a ? float(_b) : (_a);
}
//----------------------------------------------------------------------------------------------
template<class T>
inline T SIGN(const T &_a, const T &_b)
{
    return _b >= 0 ? (_a >= 0 ? _a : -_a) : (_a >= 0 ? -_a : _a);
}
//----------------------------------------------------------------------------------------------
inline float SIGN(const float &_a, const double &_b)
{
    return _b >= 0 ? (_a >= 0 ? _a : -_a) : (_a >= 0 ? -_a : _a);
}
//----------------------------------------------------------------------------------------------
inline float SIGN(const double &_a, const float &_b)
{
    return (float)(_b >= 0 ? (_a >= 0 ? _a : -_a) : (_a >= 0 ? -_a : _a));
}
//----------------------------------------------------------------------------------------------
template<class T>
inline void SWAP(T &_a, T &_b)
{
    T tmp = _a;
    _a = _b;
    _b = tmp;
}
//----------------------------------------------------------------------------------------------
//TODO: Create a file and a class for btkError
//FIXME : Compilations errors when include btkNumerical and do a throw !!

//#ifndef _USEBTKERRORCLASS_
//#define throw(message) \
//{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}
//#else
//struct BtkError
//{
//    char *message;
//    char *file;
//    int line;
//    BtkError(char *m, char *f, int l) : message(m), file(f), line(l) {}
//};
//#define throw(message) throw(BtkError(message,__FILE__,__LINE__));
//void BtkCatch(BtkError err)
//{
//    printf("ERROR: %s\n     in file %s at line %d\n",
//        err.message, err.file, err.line);
//    exit(1);
//}
//#endif
//----------------------------------------------------------------------------------------------
#ifdef _USESTDVECTOR_
#define BtkVector vector
#else
//----------------------------------------------------------------------------------------------
template <class T>
class BtkVector
{
public:
    BtkVector();
    explicit BtkVector(int n);		// Zero-based array
    BtkVector(int n, const T &a);	//initialize to constant value
    BtkVector(int n, const T *a);	// Initialize to array
    BtkVector(const BtkVector &rhs);	// Copy constructor
    BtkVector & operator=(const BtkVector &rhs);	//assignment
    typedef T value_type; // make T available externally
    inline T & operator[](const int i);	//i'th element
    inline const T & operator[](const int i) const;
    inline int size() const;
    void resize(int newn); // resize (contents not preserved)
    void assign(int newn, const T &a); // resize and assign a constant value
    ~BtkVector();

private:
    int m_size;	// size of array. upper index is m_size-1
    T *m_pData; //Pointer to data
};
//----------------------------------------------------------------------------------------------
// BtkVector definitions
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T>::BtkVector() : m_size(0), m_pData(NULL) {}
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T>::BtkVector(int n) : m_size(n), m_pData(n>0 ? new T[n] : NULL) {}
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T>::BtkVector(int n, const T& a) : m_size(n), m_pData(n>0 ? new T[n] : NULL)
{
    for(int i=0; i<n; i++)
        m_pData[i] = a;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T>::BtkVector(int n, const T *a) : m_size(n), m_pData(n>0 ? new T[n] : NULL)
{
    for(int i=0; i<n; i++)
        m_pData[i] = *a++;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T>::BtkVector(const BtkVector<T> &rhs) : m_size(rhs.m_size), m_pData(m_size>0 ? new T[m_size] : NULL)
{
    for(int i=0; i<m_size; i++)
        m_pData[i] = rhs[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T> & BtkVector<T>::operator=(const BtkVector<T> &rhs)
{
    if (this != &rhs)
    {
        if (m_size != rhs.m_size)
        {
            if (m_pData != NULL)
                delete [] (m_pData);
            m_size=rhs.m_size;
            m_pData= m_size>0 ? new T[m_size] : NULL;
        }
        for (int i=0; i<m_size; i++)
            m_pData[i]=rhs[i];
    }
    return *this;
}
//----------------------------------------------------------------------------------------------
template <class T>
inline T & BtkVector<T>::operator[](const int i)	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=m_size)
{
    throw("BtkVector subscript out of bounds");
}
#endif
    return m_pData[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
inline const T & BtkVector<T>::operator[](const int i) const	//subscripting
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=m_size)
{
    throw("BtkVector subscript out of bounds");
}
#endif
    return m_pData[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
inline int BtkVector<T>::size() const
{
    return m_size;
}
//----------------------------------------------------------------------------------------------
template <class T>
void BtkVector<T>::resize(int newn)
{
    if (newn != m_size) {
        if (m_pData != NULL)
            delete[] (m_pData);
        m_size = newn;
        m_pData = m_size > 0 ? new T[m_size] : NULL;
    }
}
//----------------------------------------------------------------------------------------------
template <class T>
void BtkVector<T>::assign(int newn, const T& a)
{
    if (newn != m_size) {
        if (m_pData != NULL)
            delete[] (m_pData);
        m_size = newn;
        m_pData = m_size > 0 ? new T[m_size] : NULL;
    }
    for (int i=0;i<m_size;i++)
        m_pData[i] = a;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkVector<T>::~BtkVector()
{
    if (m_pData != NULL)
        delete[] (m_pData);
}
//----------------------------------------------------------------------------------------------


#endif //ifdef _USESTDVECTOR_
//----------------------------------------------------------------------------------------------
template <class T>
class BtkMatrix
{

public:
    BtkMatrix();
    BtkMatrix(int n, int m);			// Zero-based array
    BtkMatrix(int n, int m, const T &a);	//Initialize to constant
    BtkMatrix(int n, int m, const T *a);	// Initialize to array
    BtkMatrix(const BtkMatrix &rhs);		// Copy constructor
    BtkMatrix & operator=(const BtkMatrix &rhs);	//assignment
    typedef T value_type; // make T available externally
    inline T* operator[](const int i);	//subscripting: pointer to row i
    inline const T* operator[](const int i) const;
    inline int nrows() const;
    inline int ncols() const;
    void resize(int newn, int newm); // resize (contents not preserved)
    void assign(int newn, int newm, const T &a); // resize and assign a constant value
    ~BtkMatrix();

private:
    int m_NSize;
    int m_MSize;
    T **m_pData;
};
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T>::BtkMatrix() : m_NSize(0), m_MSize(0), m_pData(NULL) {}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T>::BtkMatrix(int n, int m) : m_NSize(n), m_MSize(m), m_pData(n>0 ? new T*[n] : NULL)
{
    int i,nel=m*n;

    if (m_pData)
        m_pData[0] = nel>0 ? new T[nel] : NULL;

    for (i=1;i<n;i++)
        m_pData[i] = m_pData[i-1] + m;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T>::BtkMatrix(int n, int m, const T &a) : m_NSize(n), m_MSize(m), m_pData(n>0 ? new T*[n] : NULL)
{
    int i,j,nel=m*n;

    if (m_pData)
        m_pData[0] = nel>0 ? new T[nel] : NULL;

    for (i=1; i< n; i++)
        m_pData[i] = m_pData[i-1] + m;

    for (i=0; i< n; i++)
        for (j=0; j<m; j++)
            m_pData[i][j] = a;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T>::BtkMatrix(int n, int m, const T *a) : m_NSize(n), m_MSize(m), m_pData(n>0 ? new T*[n] : NULL)
{
    int i,j,nel=m*n;

    if (m_pData)
        m_pData[0] = nel>0 ? new T[nel] : NULL;

    for (i=1; i< n; i++)
        m_pData[i] = m_pData[i-1] + m;

    for (i=0; i< n; i++)
        for (j=0; j<m; j++)
            m_pData[i][j] = *a++;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T>::BtkMatrix(const BtkMatrix &rhs) : m_NSize(rhs.m_NSize), m_MSize(rhs.m_MSize), m_pData(m_NSize>0 ? new T*[m_NSize] : NULL)
{
    int i,j,nel=m_MSize*m_NSize;

    if (m_pData)
        m_pData[0] = nel>0 ? new T[nel] : NULL;

    for (i=1; i< m_NSize; i++)
        m_pData[i] = m_pData[i-1] + m_MSize;

    for (i=0; i< m_NSize; i++)
        for (j=0; j<m_MSize; j++)
            m_pData[i][j] = rhs[i][j];
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T> & BtkMatrix<T>::operator=(const BtkMatrix<T> &rhs)
{
    if (this != &rhs)
    {
        int i,j,nel;
        if (m_NSize != rhs.m_NSize || m_MSize != rhs.m_MSize)
        {
            if (m_pData != NULL)
            {
                delete[] (m_pData[0]);
                delete[] (m_pData);
            }

            m_NSize=rhs.m_NSize;
            m_MSize=rhs.m_MSize;
            m_pData = m_NSize>0 ? new T*[m_NSize] : NULL;
            nel = m_MSize*m_NSize;

            if (m_pData)
                m_pData[0] = nel>0 ? new T[nel] : NULL;

            for (i=1; i< m_NSize; i++)
                m_pData[i] = m_pData[i-1] + m_MSize;
        }

        for (i=0; i< m_NSize; i++)
            for (j=0; j<m_MSize; j++)
                m_pData[i][j] = rhs[i][j];
    }
    return *this;
}
//----------------------------------------------------------------------------------------------
template <class T>
inline T* BtkMatrix<T>::operator[](const int i)	//subscripting: pointer to row i
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=m_NSize)
{
    throw("BtkMatrix subscript out of bounds");
}
#endif
    return m_pData[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
inline const T* BtkMatrix<T>::operator[](const int i) const
{
#ifdef _CHECKBOUNDS_
if (i<0 || i>=m_NSize)
{
    throw("BtkMatrix subscript out of bounds");
}
#endif
    return m_pData[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
inline int BtkMatrix<T>::nrows() const
{
    return m_NSize;
}
//----------------------------------------------------------------------------------------------
template <class T>
inline int BtkMatrix<T>::ncols() const
{
    return m_MSize;
}
//----------------------------------------------------------------------------------------------
template <class T>
void BtkMatrix<T>::resize(int newn, int newm)
{
    int i,nel;
    if (newn != m_NSize || newm != m_MSize)
    {
        if (m_pData != NULL)
        {
            delete[] (m_pData[0]);
            delete[] (m_pData);
        }

        m_NSize = newn;
        m_MSize = newm;
        m_pData = m_NSize>0 ? new T*[m_NSize] : NULL;
        nel = m_MSize*m_NSize;

        if (m_pData)
            m_pData[0] = nel>0 ? new T[nel] : NULL;

        for (i=1; i< m_NSize; i++)
            m_pData[i] = m_pData[i-1] + m_MSize;
    }
}
//----------------------------------------------------------------------------------------------
template <class T>
void BtkMatrix<T>::assign(int newn, int newm, const T& a)
{
    int i,j,nel;
    if (newn != m_NSize || newm != m_MSize)
    {
        if (m_pData != NULL)
        {
            delete[] (m_pData[0]);
            delete[] (m_pData);
        }
        m_NSize = newn;
        m_MSize = newm;
        m_pData = m_NSize>0 ? new T*[m_NSize] : NULL;
        nel = m_MSize*m_NSize;

        if (m_pData)
            m_pData[0] = nel>0 ? new T[nel] : NULL;

        for (i=1; i< m_NSize; i++)
            m_pData[i] = m_pData[i-1] + m_MSize;
    }

    for (i=0; i< m_NSize; i++)
        for (j=0; j<m_MSize; j++)
            m_pData[i][j] = a;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMatrix<T>::~BtkMatrix()
{
    if (m_pData != NULL)
    {
        delete[] (m_pData[0]);
        delete[] (m_pData);
    }
}
//----------------------------------------------------------------------------------------------
template <class T>
class BtkMat3D
{

public:
    BtkMat3D();
    BtkMat3D(int n, int m, int k);
    inline T** operator[](const int i);	//subscripting: pointer to row i
    inline const T* const * operator[](const int i) const;
    inline int dim1() const;
    inline int dim2() const;
    inline int dim3() const;
    ~BtkMat3D();

private:
    int m_NSize;
    int m_MSize;
    int m_KSize;
    T ***m_pData;
};
//----------------------------------------------------------------------------------------------
template <class T>
BtkMat3D<T>::BtkMat3D(): m_NSize(0), m_MSize(0), m_KSize(0), m_pData(NULL) {}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMat3D<T>::BtkMat3D(int n, int m, int k) : m_NSize(n), m_MSize(m), m_KSize(k), m_pData(new T**[n])
{
    int i,j;
    m_pData[0] = new T*[n*m];
    m_pData[0][0] = new T[n*m*k];

    for(j=1; j<m; j++)
        m_pData[0][j] = m_pData[0][j-1] + k;

    for(i=1; i<n; i++)
    {
        m_pData[i] = m_pData[i-1] + m;
        m_pData[i][0] = m_pData[i-1][0] + m*k;

        for(j=1; j<m; j++)
            m_pData[i][j] = m_pData[i][j-1] + k;
    }
}
//----------------------------------------------------------------------------------------------
template <class T>
inline T** BtkMat3D<T>::operator[](const int i) //subscripting: pointer to row i
{
    return m_pData[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
inline const T* const * BtkMat3D<T>::operator[](const int i) const
{
    return m_pData[i];
}
//----------------------------------------------------------------------------------------------
template <class T>
inline int BtkMat3D<T>::dim1() const
{
    return m_NSize;
}
//----------------------------------------------------------------------------------------------
template <class T>
inline int BtkMat3D<T>::dim2() const
{
    return m_MSize;
}
//----------------------------------------------------------------------------------------------
template <class T>
inline int BtkMat3D<T>::dim3() const
{
    return m_KSize;
}
//----------------------------------------------------------------------------------------------
template <class T>
BtkMat3D<T>::~BtkMat3D()
{
    if (m_pData != NULL)
    {
        delete[] (m_pData[0][0]);
        delete[] (m_pData[0]);
        delete[] (m_pData);
    }
}
//----------------------------------------------------------------------------------------------

// Typedefs like NR, but not very helpfull
//TODO: Remove useless ones, and rename others !

typedef int Int; // 32 bit integer
typedef unsigned int Uint;

#ifdef _MSC_VER
typedef __int64 Llong; // 64 bit integer
typedef unsigned __int64 Ullong;
#else
typedef long long int Llong; // 64 bit integer
typedef unsigned long long int Ullong;
#endif

typedef char Char; // 8 bit integer
typedef unsigned char Uchar;

typedef double Doub; // default floating type
typedef long double Ldoub;

typedef complex<double> Complex; // default complex type

typedef bool Bool;

// NaN: uncomment one of the following 3 methods of defining a global NaN
// you can test by verifying that (NaN != NaN) is true

static const Doub NaN = numeric_limits<Doub>::quiet_NaN();

//Uint proto_nan[2]={0xffffffff, 0x7fffffff};
//double NaN = *( double* )proto_nan;

//Doub NaN = sqrt(-1.);

// vector types

typedef const BtkVector<Int> VecInt_I;
typedef BtkVector<Int> VecInt, VecInt_O, VecInt_IO;

typedef const BtkVector<Uint> VecUint_I;
typedef BtkVector<Uint> VecUint, VecUint_O, VecUint_IO;

typedef const BtkVector<Llong> VecLlong_I;
typedef BtkVector<Llong> VecLlong, VecLlong_O, VecLlong_IO;

typedef const BtkVector<Ullong> VecUllong_I;
typedef BtkVector<Ullong> VecUllong, VecUllong_O, VecUllong_IO;

typedef const BtkVector<Char> VecChar_I;
typedef BtkVector<Char> VecChar, VecChar_O, VecChar_IO;

typedef const BtkVector<Char*> VecCharp_I;
typedef BtkVector<Char*> VecCharp, VecCharp_O, VecCharp_IO;

typedef const BtkVector<Uchar> VecUchar_I;
typedef BtkVector<Uchar> VecUchar, VecUchar_O, VecUchar_IO;

typedef const BtkVector<Doub> VecDoub_I;
typedef BtkVector<Doub> VecDoub, VecDoub_O, VecDoub_IO;

typedef const BtkVector<Doub*> VecDoubp_I;
typedef BtkVector<Doub*> VecDoubp, VecDoubp_O, VecDoubp_IO;

typedef const BtkVector<Complex> VecComplex_I;
typedef BtkVector<Complex> VecComplex, VecComplex_O, VecComplex_IO;

typedef const BtkVector<Bool> VecBool_I;
typedef BtkVector<Bool> VecBool, VecBool_O, VecBool_IO;

// matrix types

typedef const BtkMatrix<Int> MatInt_I;
typedef BtkMatrix<Int> MatInt, MatInt_O, MatInt_IO;

typedef const BtkMatrix<Uint> MatUint_I;
typedef BtkMatrix<Uint> MatUint, MatUint_O, MatUint_IO;

typedef const BtkMatrix<Llong> MatLlong_I;
typedef BtkMatrix<Llong> MatLlong, MatLlong_O, MatLlong_IO;

typedef const BtkMatrix<Ullong> MatUllong_I;
typedef BtkMatrix<Ullong> MatUllong, MatUllong_O, MatUllong_IO;

typedef const BtkMatrix<Char> MatChar_I;
typedef BtkMatrix<Char> MatChar, MatChar_O, MatChar_IO;

typedef const BtkMatrix<Uchar> MatUchar_I;
typedef BtkMatrix<Uchar> MatUchar, MatUchar_O, MatUchar_IO;

typedef const BtkMatrix<Doub> MatDoub_I;
typedef BtkMatrix<Doub> MatDoub, MatDoub_O, MatDoub_IO;

typedef const BtkMatrix<Bool> MatBool_I;
typedef BtkMatrix<Bool> MatBool, MatBool_O, MatBool_IO;

// 3D matrix types

typedef const BtkMat3D<Doub> Mat3DDoub_I;
typedef BtkMat3D<Doub> Mat3DDoub, Mat3DDoub_O, Mat3DDoub_IO;

// Floating Point Exceptions for Microsoft compilers
//----------------------------------------------------------------------------------------------
#ifdef _TURNONFPES_
#ifdef _MSC_VER
struct turn_on_floating_exceptions
{
    turn_on_floating_exceptions()
    {
        int cw = _controlfp( 0, 0 );
        cw &=~(EM_INVALID | EM_OVERFLOW | EM_ZERODIVIDE );
        _controlfp( cw, MCW_EM );
    }
};
turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif /* _MSC_VER */
#endif /* _TURNONFPES */
//----------------------------------------------------------------------------------------------


}
//----------------------------------------------------------------------------------------------


#endif
