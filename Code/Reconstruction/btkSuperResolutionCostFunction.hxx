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

#ifndef BTKSUPERRESOLUTIONCOSTFUNCTION_H
#define BTKSUPERRESOLUTIONCOSTFUNCTION_H

#include "itkObject.h"
#include "vnl/vnl_matops.h"
#include "vnl/vnl_sparse_matrix.h"

#include "btkMacro.h"

namespace btk
{
template < class TImage >
class SuperResolutionCostFunction: public itk::Object
{
    public:
        typedef btk::SuperResolutionCostFunction< TImage >        Self;
        typedef itk::Object                             Superclass;

        typedef itk::SmartPointer< Self >               Pointer;
        typedef itk::SmartPointer< const Self >         ConstPointer;

        typedef float                                   PrecisionType;
        //typedef double                                  PrecisionType;

        typedef TImage   ImageType;

        double GetValue(const vnl_vector< double >& _x);

        void GetGradient(const vnl_vector< double >& _x, vnl_vector< double >& _g);


        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(btk::SuperResolutionCostFunction, itk::Object);


        btkSetMacro(H,vnl_sparse_matrix< PrecisionType >&);

        btkSetMacro(HtY,vnl_vector< PrecisionType >&);

        btkSetMacro(Y,vnl_vector< PrecisionType >&);

        btkSetMacro(Lambda, PrecisionType);

        btkSetMacro(NumberOfParameters,unsigned int);

        btkSetMacro(SRSize, typename ImageType::SizeType);


        struct IMG_SIZE
        {
           unsigned int width;
           unsigned int height;
           unsigned int depth;
        };


    protected:
        SuperResolutionCostFunction();
        virtual ~SuperResolutionCostFunction(){}

    private:

        void
        Set(PrecisionType * _array, int _size, PrecisionType _value)
        {
            for (unsigned int i = 0; i < _size; i++)
            {
                _array[i] = _value;
            }

        }


        /** Gets an image row as a vector. */
        void GetRow(const vnl_vector<PrecisionType>& _image, IMG_SIZE & _size, int _row, int _frame, PrecisionType * _output)
        {
            for (unsigned int i = 0; i < _size.width; i++)
            {
                _output[i] = _image[i + _row * _size.width + _frame * _size.width * _size.height];
            }
        }

        /** Sets an image row from a vector. */
        void SetRow(vnl_vector<PrecisionType>& _image, IMG_SIZE & _size, int _row, int _frame, PrecisionType * _input)
        {
            for (unsigned int i = 0; i < _size.width; i++)
            {
                _image[i + _row * _size.width + _frame * _size.width * _size.height] = _input[i];
            }
        }

        /** Gets an image column as a vector. */
        void GetCol(const vnl_vector<PrecisionType>& _image, IMG_SIZE & _size, int _col, int _frame, PrecisionType * _output)
        {
            for (unsigned int i = 0; i < _size.height; i++)
            {
                _output[i] = _image[_col + i * _size.width + _frame * _size.width * _size.height];
            }
        }

        /** Sets an image column from a vector. */
        void SetCol(vnl_vector<PrecisionType>& _image, IMG_SIZE & _size, int _col, int _frame, PrecisionType * _input)
        {
            for (unsigned int i = 0; i < _size.height; i++)
            {
                _image[_col + i * _size.width + _frame * _size.width * _size.height] = _input[i];
            }
        }

        /** Gets an image z-axis as a vector. */
        void GetSpec(const vnl_vector<PrecisionType>& _image, IMG_SIZE & _size, int _row, int _col, PrecisionType * _output)
        {
            for (unsigned int i = 0; i < _size.depth; i++)
            {
                _output[i] = _image[_col + _row * _size.width + i * _size.width * _size.height];
            }
        }

        /** Sets an image z-axis from a vector. */
        void SetSpec(vnl_vector<PrecisionType>& _image, IMG_SIZE & _size, int _row, int _col, PrecisionType * _input)
        {
            for (unsigned int i = 0; i < _size.depth; i++)
            {
                _image[_col + _row * _size.width + i * _size.width * _size.height] = _input[i];
            }
        }


        /** Mirror of the position pos. abs(pos) must not be > 2*(size-1) */
        int Mirror(int _pos, int _size);

        void Convol1d(PrecisionType * _kernel, int _ksize, PrecisionType * _src, int _src_size, PrecisionType * _dest);

        /** 3D convolution : over the rows */
        void Convol3dx(const vnl_vector<PrecisionType>& _image, vnl_vector<PrecisionType>& _image_conv,
            IMG_SIZE& _size, PrecisionType * _kernel, int _ksize);

        /** 3D convolution : over the columns */
        void Convol3dy(const vnl_vector<PrecisionType>& _image, vnl_vector<PrecisionType>& _image_conv,
            IMG_SIZE & _size, PrecisionType * _kernel, int _ksize);

        /** 3D convolution : over the spectra */
        void Convol3dz(const vnl_vector<PrecisionType>& _image, vnl_vector<PrecisionType>& _image_conv,
            IMG_SIZE & _size, PrecisionType * _kernel, int _ksize);


    private:

        vnl_sparse_matrix< PrecisionType > m_H;

        vnl_vector< PrecisionType > m_HtY;

        vnl_vector< PrecisionType > m_Y;

        float m_Lambda;

        unsigned int m_NumberOfParameters;

        IMG_SIZE        m_X_Size;

        typename ImageType::SizeType    m_SRSize;

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSuperResolutionCostFunction.txx"
#endif


#endif // BTKSUPERRESOLUTIONCOSTFUNCTION_H
