/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 13/08/2012
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

#include "btkOrientationDiffusionFunctionModel.h"


// Local includes
#include "btkSphericalHarmonics.h"
#include "btkLegendrePolynomial.h"


namespace btk
{

OrientationDiffusionFunctionModel::OrientationDiffusionFunctionModel() : m_UseSharpModel(false)
{
    // ----
}

//----------------------------------------------------------------------------------------

OrientationDiffusionFunctionModel::~OrientationDiffusionFunctionModel()
{
    // ----
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::Update()
{
    Self::Superclass::Update();

    // Create interpolate function on model image
    m_ModelImageFunction = InterpolateModelFunction::New();
    m_ModelImageFunction->SetInputImage(m_InputModelImage);

    // Set up variables
    m_NumberOfSHCoefficients  = m_InputModelImage->GetVectorLength();
    m_SphericalHarmonicsOrder = -1.5 + std::sqrt( 2.25 - 2.*(1. - m_NumberOfSHCoefficients) );

    // Compute matrices
    Self::ComputeLegendreMatrix();
    Self::ComputeSphericalHarmonicsMatrix();

    if(m_UseSharpModel)
    {
        Self::ComputeModelSharpMatrix();
    }
}

//----------------------------------------------------------------------------------------

inline OrientationDiffusionFunctionModel::ContinuousIndex OrientationDiffusionFunctionModel::TransformPhysicalPointToContinuousIndex(PhysicalPoint point)
{
    ContinuousIndex cindex;
    m_InputModelImage->TransformPhysicalPointToContinuousIndex(point, cindex);

    return cindex;
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::ModelAt(ModelImage::PixelType shCoefficients, btk::GradientDirection direction)
{
    // TODO : this code use the analytical form to compute exact signal in direction.
    // We may use interpolation to improve computation speed...


    /*/ Matrix form
    Self::Matrix Y(1, m_NumberOfSHCoefficients);

    btk::SphericalDirection u = direction.GetSphericalDirection();
    unsigned int            j = 0;

    for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
    {
        for(int m = -(int)l; m <= (int)l; m++)
        {
            Y(0,j++) = btk::SphericalHarmonics::ComputeBasis(u,l,m);
        } // for m
    } // for l

    Matrix C(m_NumberOfSHCoefficients, 1);

    for(unsigned int i = 0; i < m_NumberOfSHCoefficients; i++)
    {
        C(i,0) = shCoefficients[i];
    }

    float response = 0.f;

    if(m_UseSharpModel)
    {
        response = (Y * (m_LegendreMatrix * (m_ModelSharpMatrix * C)))(0,0);
    }
    else // m_UseSharpModel = false
    {
        response = (Y * (m_LegendreMatrix * C))(0,0);
    }

    return response;//*/


    // Direct form
    btk::SphericalDirection u = direction.GetSphericalDirection();
    float            response = 0.f;
    unsigned int            i = 0;

    if(m_UseSharpModel)
    {
        for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
        {
            for(int m = -(int)l; m <= (int)l; m++)
            {
                response += btk::SphericalHarmonics::ComputeBasis(u,l,m) * m_LegendreMatrix(i,i) * m_ModelSharpMatrix(i,i) * shCoefficients[i++];
            } // for m
        } // for l
    }
    else // m_UseSharpModel = false
    {
        for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
        {
            for(int m = -(int)l; m <= (int)l; m++)
            {
                response += btk::SphericalHarmonics::ComputeBasis(u,l,m) * m_LegendreMatrix(i,i) * shCoefficients[i++];
            } // for m
        } // for l
    }

    return response;//*/
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::ModelAt(ContinuousIndex cindex, GradientDirection direction)
{
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    return Self::ModelAt(shCoefficients, direction);
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::ModelAt(PhysicalPoint point, GradientDirection direction)
{
    return Self::ModelAt(Self::TransformPhysicalPointToContinuousIndex(point), direction);
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::ModelAt(ContinuousIndex cindex)
{
    /*/ Reuse form
    std::vector< float > response;

    for(unsigned int i = 0; i < m_Directions.size(); i++)
    {
        response.push_back(Self::ModelAt(cindex, m_Directions[i]));
    }

    return response;//*/

    /*/ Matrix form
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    Self::Matrix C(m_NumberOfSHCoefficients, 1);

    for(unsigned int i = 0; i < m_NumberOfSHCoefficients; i++)
    {
        C(i,0) = shCoefficients[i];
    }

    Self::Matrix Psi;

    if(m_UseSharpModel)
    {
        Psi = m_SphericalHarmonicsBasisMatrix * (m_LegendreMatrix * (m_ModelSharpMatrix * C));
    }
    else // m_UseSharpModel = false
    {
        Psi = m_SphericalHarmonicsBasisMatrix * (m_LegendreMatrix * C);
    }

    std::vector< float > response;

    for(unsigned int i = 0; i < Psi.Rows(); i++)
    {
        response.push_back(Psi(i,0));
    }

    return response;//*/

    // Direct form
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    std::vector< float > response;
    unsigned int numberOfDirections = m_Directions.size();

    if(m_UseSharpModel)
    {
        for(unsigned int i = 0; i < numberOfDirections; i++)
        {
            float value = 0.f;

            for(unsigned int j = 0; j < m_NumberOfSHCoefficients; j++)
            {
                value += m_SphericalHarmonicsBasisMatrix(i,j) * m_ModelSharpMatrix(j,j) * m_LegendreMatrix(j,j) * shCoefficients[j];
            }

            response.push_back(value);
        }
    }
    else // m_UseSharpModel = false
    {
        for(unsigned int i = 0; i < numberOfDirections; i++)
        {
            float value = 0.f;

            for(unsigned int j = 0; j < m_NumberOfSHCoefficients; j++)
            {
                value += m_SphericalHarmonicsBasisMatrix(i,j) * m_LegendreMatrix(j,j) * shCoefficients[j];
            }

            response.push_back(value);
        }
    }

    return response;//*/
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::ModelAt(PhysicalPoint point)
{
    return Self::ModelAt(Self::TransformPhysicalPointToContinuousIndex(point));
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::SignalAt(ModelImage::PixelType shCoefficients, btk::GradientDirection direction)
{
    // TODO : this code use the analytical form to compute exact signal in direction.
    // We may use interpolation to improve computation speed...


    /*/ Matrix form
    Self::Matrix Y(1,m_NumberOfSHCoefficients);

    btk::SphericalDirection u = direction.GetSphericalDirection();
    unsigned int            j = 0;

    for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
    {

        for(int m = -(int)l; m <= (int)l; m++)
        {
            Y(0,j++) = btk::SphericalHarmonics::ComputeBasis(u,l,m);
        } // for m
    } // for l

    Self::Matrix C(m_NumberOfSHCoefficients, 1);

    for(unsigned int i = 0; i < m_NumberOfSHCoefficients; i++)
    {
        C(i,0) = shCoefficients[i];
    }

    return (Y * C)(0,0);//*/


    // Direct form
    btk::SphericalDirection u = direction.GetSphericalDirection();
    float            response = 0.f;
    unsigned int            i = 0;

    for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
    {
        for(int m = -(int)l; m <= (int)l; m++)
        {
            response += btk::SphericalHarmonics::ComputeBasis(u,l,m) * shCoefficients[i++];
        } // for m
    } // for l

    return response;//*/
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::SignalAt(ContinuousIndex cindex, GradientDirection direction)
{
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    return Self::SignalAt(shCoefficients, direction);
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::SignalAt(PhysicalPoint point, GradientDirection direction)
{
    return Self::SignalAt(Self::TransformPhysicalPointToContinuousIndex(point), direction);
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::SignalAt(ContinuousIndex cindex)
{
    /*/ Reuse form
    std::vector< float > response;

    for(unsigned int i = 0; i < m_Directions.size(); i++)
    {
        response.push_back(Self::SignalAt(cindex, m_Directions[i]));
    }

    return response;//*/

    /*/ Matrix form
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    Self::Matrix C(m_NumberOfSHCoefficients, 1);

    for(unsigned int i = 0; i < m_NumberOfSHCoefficients; i++)
    {
        C(i,0) = shCoefficients[i];
    }

    Self::Matrix S = m_SphericalHarmonicsBasisMatrix * C;

    std::vector< float > response;

    for(unsigned int i = 0; i < S.Rows(); i++)
    {
        response.push_back(S(i,0));
    }

    return response;//*/

    // Direct form
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    std::vector< float > response;
    unsigned int numberOfDirections = m_Directions.size();

    for(unsigned int i = 0; i < numberOfDirections; i++)
    {
        float value = 0.f;

        for(unsigned int j = 0; j < m_NumberOfSHCoefficients; j++)
        {
            value += m_SphericalHarmonicsBasisMatrix(i,j) * shCoefficients[j];
        }

        response.push_back(value);
    }

    return response;//*/
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::SignalAt(PhysicalPoint point)
{
    return Self::SignalAt(Self::TransformPhysicalPointToContinuousIndex(point));
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(ContinuousIndex cindex)
{
    // TODO : test
    std::vector< float > Psi = Self::ModelAt(cindex);

    float min = Psi[0], max = Psi[0];
    for(unsigned int i = 1; i < Psi.size(); i++)
    {
        if(min > Psi[i])
            min = Psi[i];

        if(max < Psi[i])
            max = Psi[i];
    }

    // Compute the regular step on the unit sphere
    unsigned int elevationResolution = static_cast< unsigned int >(std::ceil( std::sqrt(static_cast< float >(m_SphericalResolution)/2.f) ));
    float                       step = M_PI / static_cast< float >(elevationResolution);

    std::vector< btk::GradientDirection > meanDirections;

    unsigned int thetaRes = step+1;
    unsigned int phiRes   = step;


    for(unsigned int i = 1; i < thetaRes-1; i++)
    {
        for(unsigned int j = 0; j < phiRes; j++)
        {
            float odf1, odf2, odf3, odf4, odf5, odf6, odf7, odf8;
            float odf = Psi[phiRes*i + j];

            if(j == 0)
            {
                odf1 = Psi[phiRes*(i-1)+(phiRes-1)];
                odf2 = Psi[phiRes*(i-1)+j];
                odf3 = Psi[phiRes*(i-1)+(j+1)];
                odf4 = Psi[phiRes*i+(phiRes-1)];
                odf5 = Psi[phiRes*i+(j+1)];
                odf6 = Psi[phiRes*(i+1)+(phiRes-1)];
                odf7 = Psi[phiRes*(i+1)+j];
                odf8 = Psi[phiRes*(i+1)+(j+1)];
            }
            else if(j == phiRes-1)
            {
                odf1 = Psi[phiRes*(i-1)+(j-1)];
                odf2 = Psi[phiRes*(i-1)+j];
                odf3 = Psi[phiRes*(i-1)+(0)];
                odf4 = Psi[phiRes*i+(j-1)];
                odf5 = Psi[phiRes*i+(0)];
                odf6 = Psi[phiRes*(i+1)+(j-1)];
                odf7 = Psi[phiRes*(i+1)+j];
                odf8 = Psi[phiRes*(i+1)+(0)];
            }
            else // j != 0 && j != phiRes-1
            {
                odf1 = Psi[phiRes*(i-1)+(j-1)];
                odf2 = Psi[phiRes*(i-1)+j];
                odf3 = Psi[phiRes*(i-1)+(j+1)];
                odf4 = Psi[phiRes*i+(j-1)];
                odf5 = Psi[phiRes*i+(j+1)];
                odf6 = Psi[phiRes*(i+1)+(j-1)];
                odf7 = Psi[phiRes*(i+1)+j];
                odf8 = Psi[phiRes*(i+1)+(j+1)];
            }

            if(odf > odf1 && odf > odf2 && odf > odf3 && odf > odf4 && odf > odf5 && odf > odf6 && odf > odf7 && odf > odf8)
            {
                if((Psi[phiRes*i+j]-min)/(max-min) > 0.9)
                    meanDirections.push_back(btk::GradientDirection(i*step,j*step));
            }
        } // for each phi
    } // for each theta

    return meanDirections;
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(ContinuousIndex cindex, GradientDirection vector, float angle)
{
    std::vector< btk::GradientDirection > meanDirections = Self::MeanDirectionsAt(cindex);

    std::vector< btk::GradientDirection > restrictedMeanDirections;

    for(unsigned int i = 0; i < meanDirections.size(); i++)
    {
        float dotProduct = meanDirections[i]*vector;
        float      alpha = std::acos( dotProduct / (meanDirections[i].GetNorm()*vector.GetNorm()) );

        // Check if the current direction is in the solid angle
        if(alpha <= angle)
        {
            restrictedMeanDirections.push_back(meanDirections[i]);
        }
    }

    return restrictedMeanDirections;
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(PhysicalPoint point)
{
    return Self::MeanDirectionsAt(Self::TransformPhysicalPointToContinuousIndex(point));
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(PhysicalPoint point, GradientDirection vector, float angle)
{
    return Self::MeanDirectionsAt(Self::TransformPhysicalPointToContinuousIndex(point), vector, angle);
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::UseSharpModelOn()
{
    if(m_SphericalHarmonicsOrder <= 8)
    {
        m_UseSharpModel = true;
    }
    else
    {
        m_UseSharpModel = false;
    }
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::UseSharpModelOff()
{
    m_UseSharpModel = false;
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::ComputeLegendreMatrix()
{
    // Resize the matrix (rows: number of SH coefficients, columns: number of SH coefficients).
    m_LegendreMatrix = Self::Matrix(m_NumberOfSHCoefficients, m_NumberOfSHCoefficients);

    // This matrix is diagonal, so we fill the other coefficients with 0.
    m_LegendreMatrix.Fill(0);

    // Compute the diagonal coefficients.
    unsigned int i = 0;

    for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
    {
        for(int m = -(int)l; m <= (int)l; m++)
        {
            m_LegendreMatrix(i,i) = 2.0 * M_PI * btk::LegendrePolynomial::ComputeInZero(l);
            i++;
        } // for m
    } // for l
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::ComputeSphericalHarmonicsMatrix()
{
    // Resize the matrix (rows: number of gradient directions, columns: number of SH coefficients).
    m_SphericalHarmonicsBasisMatrix = Self::Matrix(m_Directions.size(), m_NumberOfSHCoefficients);

    // Compute the basis
    for(unsigned int u = 0; u < m_Directions.size(); u++)
    {
        unsigned int j = 0;

        for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
        {
            for(int m = -(int)l; m <= (int)l; m++)
            {
                m_SphericalHarmonicsBasisMatrix(u,j++) = btk::SphericalHarmonics::ComputeBasis(m_Directions[u].GetSphericalDirection(), l, m);
            } // for each m
        } // for each even order
    } // for each gradient direction
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::ComputeModelSharpMatrix()
{
    m_ModelSharpMatrix = Self::Matrix(m_NumberOfSHCoefficients,m_NumberOfSHCoefficients);
    m_ModelSharpMatrix.Fill(0);

    double lambda1 = 1.0;
    double lambda2 = 0.1;

    double alpha = 1.0 - lambda2/lambda1;

    double alpha_2 = alpha * alpha;
    double alpha_3 = alpha_2 * alpha;
    double alpha_4 = alpha_3 * alpha;

    double sqrtAlpha                  = std::sqrt(alpha);
    double asinSqrtAlpha              = std::asin(sqrtAlpha);
    double sqrtOneMinusAlpha          = std::sqrt(1.0 - alpha);
    double sqrtOneMinusAlphaSqrtAlpha = sqrtOneMinusAlpha * sqrtAlpha;

    double coefficient = 4.0 * m_BValue * std::sqrt(lambda2 * lambda1);

    if(m_SphericalHarmonicsOrder >= 0)
    {
        double coefficientOfOrder = (1.0 / std::pow(alpha,0.5)) * 2.0 * asinSqrtAlpha;
        m_ModelSharpMatrix(0,0) = coefficient / coefficientOfOrder;
    }

    if(m_SphericalHarmonicsOrder >= 2)
    {
        double coefficientOfOrder = (-1.0 / ( 2.0 * std::pow(alpha,1.5) )) *
                             (
                                 (-3.0 + 2.0 * alpha) * asinSqrtAlpha + 3.0 * sqrtOneMinusAlphaSqrtAlpha
                             );

        for(unsigned int i = 1; i < 6; i++)
        {
            m_ModelSharpMatrix(i,i) = coefficient / coefficientOfOrder;
        }
    }

    if(m_SphericalHarmonicsOrder >= 4)
    {
        double coefficientOfOrder = (1.0 / ( 32.0 * std::pow(alpha,2.5) )) *
                             (
                                 (105.0 - 102.0 * alpha + 24.0 * alpha_2) * asinSqrtAlpha + (-105.0 + 50.0 * alpha) * sqrtOneMinusAlphaSqrtAlpha
                             );

        for(unsigned int i = 6; i < 15; i++)
        {
            m_ModelSharpMatrix(i,i) = coefficient / coefficientOfOrder;
        }
    }

    if(m_SphericalHarmonicsOrder >= 6)
    {
        double coefficientOfOrder = (-1.0 / ( 128.0 * std::pow(alpha,3.5) )) *
                             (
                                 (-1155.0 + 1890.0 * alpha - 840.0 * alpha_2 + 80.0 * alpha_3) * asinSqrtAlpha + (1155.0 - 1120.0 * alpha + 196.0 * alpha_2) * sqrtOneMinusAlphaSqrtAlpha
                             );

        for(unsigned int i = 15; i < 28; i++)
        {
            m_ModelSharpMatrix(i,i) = coefficient / coefficientOfOrder;
        }
    }

    if(m_SphericalHarmonicsOrder >= 8)
    {
        double coefficientOfOrder = (1.0 / ( 8192.0 * std::pow(alpha,4.5) )) *
                             (
                                 (225225.0 - 480480.0 * alpha + 332640.0 * alpha_2 - 80640.0 * alpha_3 + 4480.0 * alpha_4) * asinSqrtAlpha + (-225225.0 + 330330.0 * alpha - 132440.0 * alpha_2 + 12176.0 * alpha_3) * sqrtOneMinusAlphaSqrtAlpha
                             );

        for(unsigned int i = 28; i < 45; i++)
        {
            m_ModelSharpMatrix(i,i) = coefficient / coefficientOfOrder;
        }
    }
}

} // namespace btk
