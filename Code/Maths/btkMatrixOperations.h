/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 28/10/2011
  Author(s): Julien Pontabry (pontabry@unistra.fr)
  
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

#ifndef BTK_MATRIX_OPERATIONS_H
#define BTK_MATRIX_OPERATIONS_H

// ITK includes
#include "itkVariableSizeMatrix.h"

namespace btk
{

/**
 * @brief Defines mathematical operations on matrices.
 * @author Julien Pontabry
 * @ingroup Maths
 */
class MatrixOperations
{
    public:
        typedef MatrixOperations                 Self;
        typedef itk::VariableSizeMatrix< double > Matrix;

        /**
         * @brief Compute the spectral norm (or the euclidian norm) of a matrix
         * @param matrix Matrix to compute the norm
         * @return Euclidian norm of the matrix
         */
        static double Norm(Self::Matrix &matrix);

        /**
         * @brief Compute the matrix exponential
         * @param matrix Matrix to compute the exponential
         * @param maxIterations Number of iterations of Padé approximation with scaling
         * @return Matrix which is the exponential of parameter matrix
         */
        static Self::Matrix Exponential(Self::Matrix &matrix, unsigned int maxIterations=6);

        /**
         * @brief Compute the matrix logarithm
         * @param matrix Matrix to compute the logarithm
         * @param epsilon Convergence threshold
         * @param maxNbOfIterations Maximum number of iterations (if numerical problems leads to poor convergence, it stops the algorithm before falling into infinite loops)
         * @return Matrix which is the logarithm of parameter matrix
         */
        static Self::Matrix Logarithm(Self::Matrix &matrix, double epsilon=0.0001, unsigned int maxNbOfIterations=10);

        /**
         * @brief Compute the matrix square root
         * @param matrix Matrix to compute the square root
         * @param epsilon Convergence threshold
         * @param maxNbOfIterations Maximum number of iterations (if numerical problems leads to poor convergence, it stops the algorithm before falling into infinite loops)
         * @return Matrix which is the square root of parameter matrix
         */
        static Self::Matrix Sqrt(Self::Matrix &matrix, double epsilon=0.0001, unsigned int maxNbOfIterations=10);
};

} // namespace btk

#endif // BTK_MATRIX_OPERATIONS_H
