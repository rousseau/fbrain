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

#include "btkMatrixOperations.h"


// STL includes
#include "cmath"
#include "complex"

// VNL (ITK) includes
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/vnl_power.h"


namespace btk
{

double MatrixOperations::Norm(Self::Matrix &matrix)
{
    // Get the conjugate transpose matrix
    Self::Matrix matrixT, matrixTmatrix;
    matrixT = matrix.GetTranspose();
    matrixTmatrix = matrixT * matrix;

    // Compute eigen decomposition
    vnl_real_eigensystem eig(matrixTmatrix.GetVnlMatrix());

    // Search maximal eigen value
    float lambdaMax = 0;

    for(unsigned int i = 0; i < matrix.Rows(); i++)
    {
        float lambda = eig.D(i).real();

        if(lambda > lambdaMax)
            lambdaMax = lambda;
    }

    return std::sqrt(lambdaMax);
}

//-----------------------------------------------------------------------------------------------------------

MatrixOperations::Matrix MatrixOperations::Exponential(Self::Matrix &matrix, unsigned int maxIterations)
{
    float j = std::max( 0.0, 1.0 + std::floor(std::log(MatrixOperations::Norm(matrix))/std::log(2)) );
    Self::Matrix A = matrix * std::pow(2,-j);

    Self::Matrix D(matrix.Rows(),matrix.Cols()), N(matrix.Rows(),matrix.Cols()), X(matrix.Rows(),matrix.Cols());
    D.SetIdentity(); N.SetIdentity(); X.SetIdentity();
    float c = 1;

    for(unsigned int k = 1; k < maxIterations; k++)
    {
        c = c * ((float)maxIterations - (float)k + 1.0) / ((float)k * (2.0*(float)maxIterations - (float)k + 1.0));
        X = A * X;
        N = N + X * c;
        D = D + X * c * ((k % 2 == 0) ? 1.0 : -1.0);
    }

    X = D.GetInverse() * N.GetVnlMatrix();

    Self::Matrix expM;
    expM = vnl_power(X.GetVnlMatrix(),(std::pow(2,j)));

    return expM;
}

//-----------------------------------------------------------------------------------------------------------

MatrixOperations::Matrix MatrixOperations::Logarithm(Self::Matrix &matrix, double epsilon, unsigned int maxNbOfIterations)
{
    unsigned int   k = 0;
    Self::Matrix I(matrix.Rows(),matrix.Cols()); I.SetIdentity();
    Self::Matrix   A = matrix;
    Self::Matrix AmI = A - I;

    double     normMatrixMinusI = MatrixOperations::Norm(AmI);
    unsigned int nbOfIterations = 0;

    while(normMatrixMinusI > 0.5 && nbOfIterations <= maxNbOfIterations)
    {
        A = MatrixOperations::Sqrt(A,epsilon);
        k++;
        AmI = A - I;
        normMatrixMinusI = MatrixOperations::Norm(AmI);
        nbOfIterations++;
    }

    A = I - A;
    Self::Matrix Z = A;
    Self::Matrix X = A;
    unsigned int i = 1;

    float    normZ = MatrixOperations::Norm(Z);
    nbOfIterations = 0;

    while(normZ > epsilon  && nbOfIterations <= maxNbOfIterations)
    {
        Z = Z * A;
        i++;
        X = X + Z/(float)i;
        normZ = MatrixOperations::Norm(Z);
        nbOfIterations++;
    }

    X = X * std::pow(2.0,(float)k);

    return X*(-1);
}

//-----------------------------------------------------------------------------------------------------------

MatrixOperations::Matrix MatrixOperations::Sqrt(Self::Matrix &matrix, double epsilon, unsigned int maxNbOfIterations)
{
    Self::Matrix X = matrix;
    Self::Matrix Y(matrix.Rows(),matrix.Cols()); Y.SetIdentity();

    unsigned int nbOfIterations = 0;
    Self::Matrix X2mMatrix = X * X - matrix;
    Self::Matrix iX(X.Rows(),X.Cols()), iY(Y.Rows(),Y.Cols());

    while(MatrixOperations::Norm(X2mMatrix) > epsilon  && nbOfIterations <= maxNbOfIterations)
    {
        iX = X.GetInverse();
        iY = Y.GetInverse();

        X = (X + iY) * 0.5;
        Y = (Y + iX) * 0.5;

        X2mMatrix = X * X - matrix;
        nbOfIterations++;
    }

    return X;
}

} // namespace btk
