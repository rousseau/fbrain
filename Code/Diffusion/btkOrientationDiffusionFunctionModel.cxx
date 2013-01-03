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


#define ACCESS(i,j) (phiRes*(i) + (j) + 1)


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
    this->ComputeLegendreMatrix();
    this->ComputeSphericalHarmonicsMatrix();

    if(m_UseSharpModel)
    {
        this->ComputeModelSharpMatrix();
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

    return response;
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::ModelAt(ContinuousIndex cindex, GradientDirection direction)
{
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    return this->ModelAt(shCoefficients, direction);
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::ModelAt(PhysicalPoint point, GradientDirection direction)
{
    return this->ModelAt(this->TransformPhysicalPointToContinuousIndex(point), direction);
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::ModelAt(ContinuousIndex cindex)
{
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

    return response;
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::ModelAt(PhysicalPoint point)
{
    return this->ModelAt(this->TransformPhysicalPointToContinuousIndex(point));
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::SignalAt(ModelImage::PixelType shCoefficients, btk::GradientDirection direction)
{
    // TODO : this code use the analytical form to compute exact signal in direction.
    // We may use interpolation to improve computation speed...

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

    return response;
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::SignalAt(ContinuousIndex cindex, GradientDirection direction)
{
    ModelImage::PixelType shCoefficients = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    return this->SignalAt(shCoefficients, direction);
}

//----------------------------------------------------------------------------------------

float OrientationDiffusionFunctionModel::SignalAt(PhysicalPoint point, GradientDirection direction)
{
    return this->SignalAt(this->TransformPhysicalPointToContinuousIndex(point), direction);
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::SignalAt(ContinuousIndex cindex)
{
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

    return response;
}

//----------------------------------------------------------------------------------------

std::vector< float > OrientationDiffusionFunctionModel::SignalAt(PhysicalPoint point)
{
    return this->SignalAt(this->TransformPhysicalPointToContinuousIndex(point));
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(ContinuousIndex cindex)
{
    // TODO: improve speed with analytical space reduction (or interpolation of pre-computed data)

    std::vector< float > Psi = this->ModelAt(cindex);

    float min = Psi[0], max = Psi[0];
    for(unsigned int i = 1; i < Psi.size(); i++)
    {
        if(min > Psi[i])
            min = Psi[i];

        if(max < Psi[i])
            max = Psi[i];
    }

    // Compute the regular step on the unit sphere
    unsigned int thetaRes = static_cast< unsigned int >(std::ceil( std::sqrt(static_cast< float >(m_SphericalResolution)/2.f) ));
    unsigned int   phiRes = 2*thetaRes;

    std::vector< btk::GradientDirection > meanDirections;

    for(unsigned int i = 0; i < thetaRes; i++)
    {
        for(unsigned int j = 0; j < phiRes; j++)
        {
            float odf1, odf2, odf3, odf4, odf5, odf6, odf7, odf8;
            float odf = Psi[ACCESS(i,j)];

            if(j == 0)
            {
                if(i == 0)
                {
                    odf1 = odf2 = odf3 = Psi[0];
                }
                else // i != 0
                {
                    odf1 = Psi[ACCESS(i-1,phiRes-1)];
                    odf2 = Psi[ACCESS(i-1,0)];
                    odf3 = Psi[ACCESS(i-1,1)];
                }

                odf4 = Psi[ACCESS(i,phiRes-1)];
                odf5 = Psi[ACCESS(i,1)];

                if(i == thetaRes-1)
                {
                    odf6 = odf7 = odf8 = Psi[Psi.size()-1];
                }
                else // i != thetaRes-1
                {
                    odf6 = Psi[ACCESS(i+1,phiRes-1)];
                    odf7 = Psi[ACCESS(i+1,0)];
                    odf8 = Psi[ACCESS(i+1,1)];
                }
            }
            else if(j == phiRes-1)
            {
                if(i == 0)
                {
                    odf1 = odf2 = odf3 = Psi[0];
                }
                else // i != 0
                {
                    odf1 = Psi[ACCESS(i-1,j-1)];
                    odf2 = Psi[ACCESS(i-1,j)];
                    odf3 = Psi[ACCESS(i-1,0)];
                }

                odf4 = Psi[ACCESS(i,j-1)];
                odf5 = Psi[ACCESS(i,0)];

                if(i == thetaRes-1)
                {
                    odf6 = odf7 = odf8 = Psi[Psi.size()-1];
                }
                else // i != thetaRes-1
                {
                    odf6 = Psi[ACCESS(i+1,j-1)];
                    odf7 = Psi[ACCESS(i+1,j)];
                    odf8 = Psi[ACCESS(i+1,j+1)];
                }
            }
            else // j != 0 && j != phiRes-1
            {
                if(i == 0)
                {
                    odf1 = odf2 = odf3 = Psi[0];
                }
                else // i != 0
                {
                    odf1 = Psi[ACCESS(i-1,j-1)];
                    odf2 = Psi[ACCESS(i-1,j)];
                    odf3 = Psi[ACCESS(i-1,j+1)];
                }

                odf4 = Psi[ACCESS(i,j-1)];
                odf5 = Psi[ACCESS(i,j+1)];

                if(i == thetaRes-1)
                {
                    odf6 = odf7 = odf8 = Psi[Psi.size()-1];
                }
                else // i != thetaRes-1
                {
                    odf6 = Psi[ACCESS(i+1,j-1)];
                    odf7 = Psi[ACCESS(i+1,j)];
                    odf8 = Psi[ACCESS(i+1,j+1)];
                }
            }

            if(odf > odf1 && odf > odf2 && odf > odf3 && odf > odf4 && odf > odf5 && odf > odf6 && odf > odf7 && odf > odf8)
            {
                if((odf-min)/(max-min) > 0.9)
                {
                    meanDirections.push_back(m_Directions[ACCESS(i,j)]);
                }
            }
        } // for each phi
    } // for each theta

    // North pole
    bool maxima = true;
    unsigned int k = 1;

    float odf = Psi[0];
    while(maxima && k <= phiRes)
    {
        if(odf < Psi[k++])
        {
            maxima = false;
        }
    }

    if(maxima && (odf-min)/(max-min) > 0.9)
    {
        meanDirections.push_back(m_Directions[0]);
    }


    // South pole
    maxima = true;
    k = Psi.size()-1;

    odf = Psi[Psi.size()-1];
    while(maxima && k >= Psi.size()-1-phiRes)
    {
        if(odf < Psi[k--])
        {
            maxima = false;
        }
    }

    if(maxima && (odf-min)/(max-min) > 0.9)
    {
        meanDirections.push_back(m_Directions[m_Directions.size()-1]);
    }

    return meanDirections;
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(ContinuousIndex cindex, GradientDirection vector, float angle)
{
    std::vector< btk::GradientDirection > meanDirections = this->MeanDirectionsAt(cindex);

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
    return this->MeanDirectionsAt(this->TransformPhysicalPointToContinuousIndex(point));
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > OrientationDiffusionFunctionModel::MeanDirectionsAt(PhysicalPoint point, GradientDirection vector, float angle)
{
    return this->MeanDirectionsAt(this->TransformPhysicalPointToContinuousIndex(point), vector, angle);
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void OrientationDiffusionFunctionModel::UseSharpModelOn()
{
    m_UseSharpModel = true;
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
