/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

15 june 2010
< pontabry at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#include "btkDTFPImportanceDensity.h"


// STL includes
#include "cassert"
#include "cmath"

// ITK includes
#include "itkDiffusionTensor3D.h"

// default 0.25
#define Tau 0.25


namespace btk
{

DTFPImportanceDensity::DTFPImportanceDensity(DTFPSignal *model)
{
    assert(model);


    m_model          = model;
    m_2PI            = 2.0 * M_PI;
    m_log4PI         = std::log(4.0 * M_PI);
}

Direction DTFPImportanceDensity::computeMeanDirection(Point xk, Direction ukm1)
{
    itk::DiffusionTensor3D<Real> tensor = m_model->DiffusionTensorAt(xk);

    itk::FixedArray<Real,3> eigenValues;
    itk::Matrix<Real,3,3> eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues,eigenVectors);
//    Pr(eigenValues[0]); Pr(eigenValues[1]); Pr(eigenValues[2]);
//    Pr(eigenVectors(0,0)); Pr(eigenVectors(0,1)); Pr(eigenVectors(0,2));
//    Pr(eigenVectors(1,0)); Pr(eigenVectors(1,1)); Pr(eigenVectors(1,2));
//    Pr(eigenVectors(2,0)); Pr(eigenVectors(2,1)); Pr(eigenVectors(2,2));

    Real cl = (eigenValues[2] - eigenValues[1]) / std::sqrt(eigenValues[2]*eigenValues[2] + eigenValues[1]*eigenValues[1] + eigenValues[0]*eigenValues[0]);

    if(cl > Tau)
    {
        Vector vkm1     = ukm1.toVector();
        Real dotProduct = vkm1.x()*eigenVectors(2,0) + vkm1.y()*eigenVectors(2,1) + vkm1.z()*eigenVectors(2,2);
        Real norm       = std::sqrt( (vkm1.x()*vkm1.x() + vkm1.y()*vkm1.y() + vkm1.z()*vkm1.z()) * (eigenVectors(2,0)*eigenVectors(2,0) + eigenVectors(2,1)*eigenVectors(2,1) + eigenVectors(2,2)*eigenVectors(2,2)) );
        Real angle      = std::acos(dotProduct/norm);

        if(angle <= M_PI/2.0)
            return Vector(eigenVectors(2,0), eigenVectors(2,1), eigenVectors(2,2)).toDirection();
        else
            return Vector(-eigenVectors(2,0), -eigenVectors(2,1), -eigenVectors(2,2)).toDirection();
    }
    else // cl <= 0.25
    {
        return ukm1;
    }
}

Direction DTFPImportanceDensity::simulate(Direction mean, Real kappa)
{
    // Sample random scalar
    Real y = (Real)std::rand() / (Real)RAND_MAX;
    Real w = 1.0/kappa * std::log(std::exp(-kappa) + kappa * 2.0/kappa * std::sinh(kappa) * y);

    // Sample random angle (to get random unit vector)
    Real angle = ((Real)std::rand() / (Real)RAND_MAX) * m_2PI;

    // Concatenate to obtain unit vector with vmf distribution with mean = (0,0,1)
    Real cst = std::sqrt(1-w*w);
    Vector X(cst*std::cos(angle), cst*std::sin(angle), w);

    Real theta, phi;
    Vector vmean = mean.toVector();

    // Rotate if necessary
    if(!(vmean.x() == 0 && vmean.y() == 0 && vmean.z() == 1)) // mean not (0,0,1)
    {
        // Rotations from (0,0,1) to mean
        Real cosTheta = std::cos(mean.theta());
        Real sinTheta = std::sin(mean.theta());
        Real cosPhi   = std::cos(mean.phi());
        Real sinPhi   = std::sin(mean.phi());

        // Ry(theta)
        // | cos(theta)  0  sin(theta)  |   |x|
        // |     0       1       0      | . |y|
        // |-sin(theta)  0   cos(theta) |   |z|
        Vector tmp(
                cosTheta * X.x() + sinTheta * X.z(),
                X.y(),
                -sinTheta * X.x() + cosTheta * X.z()
        );

        // Rz(phi)
        // | cos(phi)  -sin(phi)  0 |   |x|
        // | sin(phi)   cos(phi)  0 | . |y|
        // |    0          0      1 |   |z|
        Vector tmp2(
                cosPhi * tmp.x() - sinPhi * tmp.y(),
                sinPhi * tmp.x() + cosPhi * tmp.y(),
                tmp.z()
        );

        theta = tmp2.toSphericalCoordinates().theta();
        phi   = tmp2.toSphericalCoordinates().phi();
    }
    else // vkm1 is (0,0,1)
    {
        theta = X.toSphericalCoordinates().theta();
        phi   = X.toSphericalCoordinates().phi();
    }

    return Direction(theta,phi);
}

Real DTFPImportanceDensity::compute(Direction uk, Direction mean, Real kappa)
{
    Vector vk   = uk.toVector();
    Vector vmean = mean.toVector();

    return std::log(kappa) - m_log4PI - std::log(std::sinh(kappa)) + kappa * (vk.x()*vmean.x() + vk.y()*vmean.y() + vk.z()*vmean.z());
}

Real DTFPImportanceDensity::computeConcentration(Direction mu, Point xk)
{
    itk::DiffusionTensor3D<Real> tensor = m_model->DiffusionTensorAt(xk);

    itk::FixedArray<Real,3> eigenValues;
    tensor.ComputeEigenValues(eigenValues);
//    Pr(eigenValues[0]); Pr(eigenValues[1]); Pr(eigenValues[2]);

    Real cl = (eigenValues[2] - eigenValues[1]) / std::sqrt(eigenValues[2]*eigenValues[2] + eigenValues[1]*eigenValues[1] + eigenValues[0]*eigenValues[0]);

    if(cl > Tau)
    {
        Real alpha = 15;
        Real gamma = 0.42;
        Real FA    = tensor.GetFractionalAnisotropy();
        Real expon = FA/gamma;

        return alpha + std::exp(expon*expon);
    }
    else // cl <= Tau
    {
        return 30;
    }
}

} // namespace btk

