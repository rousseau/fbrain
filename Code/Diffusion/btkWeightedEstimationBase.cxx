/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/06/2013
  Author(s): Frederic Champ (champ(at)unistra.fr)
  
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

#include "btkWeightedEstimationBase.h"


// Local includes
#include "btkSphericalHarmonics.h"

// ITK includes
#include "math.h"
#include "float.h"
#include "algorithm"


namespace btk
{

WeightedEstimationBase::WeightedEstimationBase() : Superclass()
{
    m_SphericalHarmonicsOrder = 4;
    m_RegularizationParameter = 0.006;
    m_Sigma = 1.0;
    m_SphericalResolution = 100;
}

//----------------------------------------------------------------------------------------

WeightedEstimationBase::~WeightedEstimationBase()
{
    // ----
}

//----------------------------------------------------------------------------------------

void WeightedEstimationBase::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void WeightedEstimationBase::ComputeWeightedMatrix()
{

    m_WeightsMatrix = Self::Matrix(m_numberOfNeighbors, m_numberOfNeighbors);
    m_WeightsMatrix.Fill(0.0);

    for(unsigned int i = 0; i<m_numberOfNeighbors; i++)
    {

        float SpatialDistance=0.0;
        double    weight;
        ////////////////////////////////////////////////////////////////////////////
        //
        // Normalized spatial distance
        //
        if(m_DMax !=  m_DMin)
        {
            SpatialDistance = ( m_DMax - m_NeighborsDistances[i]) / (m_DMax-m_DMin);
        }

        ////////////////////////////////////////////////////////////////////////////
        //
        // kernel computation
        //
        SpatialDistance/=m_Sigma; // Fixed at 1 ...

        //Gaussian kernel
        weight = std::exp(-0.5*  (SpatialDistance*SpatialDistance));

        m_WeightsMatrix[i][i]= weight;
    }

}


//----------------------------------------------------------------------------------------
void WeightedEstimationBase::ComputeSphericalHarmonicsBasisMatrix()
{

    // Resize the matrix (rows: number of gradient directions, columns: number of SH coefficients).
    m_SphericalHarmonicsBasisMatrix = Self::Matrix(m_numberOfNeighbors, m_NumberOfSHCoefficients);

    // Compute the basis
    for(unsigned int u = 0; u < m_numberOfNeighbors; u++)
    {
        unsigned int j = 0;

        for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
        {
            for(int m = -(int)l; m <= (int)l; m++)
            {
                m_SphericalHarmonicsBasisMatrix(u,j++) = btk::SphericalHarmonics::ComputeBasis(m_GradientTable[u].GetSphericalDirection(), l, m);
            } // for each m
        } // for each even order
    } // for each gradient direction
}

//----------------------------------------------------------------------------------------

void WeightedEstimationBase::ComputeRegularizationMatrix()
{
    // Resize the matrix (rows: number of SH coefficients, columns: number of SH coefficients).
    m_RegularizationMatrix = Self::Matrix(m_NumberOfSHCoefficients, m_NumberOfSHCoefficients);

    // This matrix is diagonal, so we fill the other coefficients with 0.
    m_RegularizationMatrix.Fill(0);

    // Compute the diagonal coefficients.
    unsigned int i = 0;

    for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
    {
        float   lp1 = l+1;
        float value = l*l * lp1*lp1;

        for(int m = -(int)l; m <= (int)l; m++)
        {
            m_RegularizationMatrix(i,i) = value;
            i++;
        } // for each m
    } // for each even order
}

//----------------------------------------------------------------------------------------

void WeightedEstimationBase::ComputeTransitionMatrix()
{

    // Resize  the matrix (rows: number of SH coefficients, columns: number of gradient directions).
    m_TransitionMatrix = Matrix(m_NumberOfSHCoefficients, m_numberOfNeighbors);

    // Compute the transpose of the spherical harmonics basis matrix.
    Self::Matrix SphericalHarmonicsBasisTranspose;
    SphericalHarmonicsBasisTranspose = m_SphericalHarmonicsBasisMatrix.GetTranspose();

    // Compute the inverse.
    Self::Matrix Inverse;
    Inverse = ((SphericalHarmonicsBasisTranspose* m_WeightsMatrix * m_SphericalHarmonicsBasisMatrix) + (m_RegularizationMatrix * m_RegularizationParameter)).GetInverse();

    // Compute the transition matrix.
    m_TransitionMatrix = Inverse * SphericalHarmonicsBasisTranspose * m_WeightsMatrix;


}


//----------------------------------------------------------------------------------------
void WeightedEstimationBase::Initialize()
{

    m_numberOfNeighbors = m_NeighborsDistances.size();

    // Compute the number of coefficients
    m_NumberOfSHCoefficients = 0.5 * (m_SphericalHarmonicsOrder+1) * (m_SphericalHarmonicsOrder+2);

    // Compute needed matrices
    this->ComputeSphericalHarmonicsBasisMatrix();
    this->ComputeRegularizationMatrix();

    //btkCoutMacro("SignaSignalVectorlVector & ComputeRegularizationMatrix = done.");
    m_SignalMatrix = Matrix(m_NeighborsDistances.size(), 1);

    for(unsigned int i = 0; i < m_SignalMatrix.Rows(); i++)
    {
        m_SignalMatrix(i,0) = m_SignalValues[i]; // / static_cast< float >(signal[0]);
    }



}

//----------------------------------------------------------------------------------------

float WeightedEstimationBase::SignalAt(btk::GradientDirection direction)
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
            response += btk::SphericalHarmonics::ComputeBasis(u,l,m) * m_Coefficients[i];
            i++;
        } // for m
    } // for l

    return response;
}

//----------------------------------------------------------------------------------------
void WeightedEstimationBase::Update()
{

    ////////////////////////////////////////////////////////////////////////////
    //
    // Compute the regular step on the unit sphere for model estimation
    //
    unsigned int elevationResolution = static_cast< unsigned int >(std::ceil( std::sqrt(static_cast< float >(m_SphericalResolution)/2.f) ));
    float                       step = M_PI / static_cast< float >(elevationResolution);
    const double M_2MPI = 2.0 * M_PI;

    ////////////////////////////////////////////////////////////////////////////
    //
    // Max-Min search of spatial distance
    //
    m_DMax = m_NeighborsDistances[0];
    m_DMin = m_NeighborsDistances[0];

    for(unsigned int j = 0; j < m_NeighborsDistances.size(); j++)
    {

        if(m_NeighborsDistances[j]>m_DMax)
        {
            m_DMax=m_NeighborsDistances[j];
        }
        if(m_NeighborsDistances[j]<m_DMin)
        {
            m_DMin=m_NeighborsDistances[j];
        }

    }


    ////////////////////////////////////////////////////////////////////////////
    //
    // Compute Matrix
    //
    this->ComputeWeightedMatrix();
    this->ComputeTransitionMatrix();

    Self::Matrix coefficientsMatrix(m_NumberOfSHCoefficients, 1);
    coefficientsMatrix = m_TransitionMatrix * m_SignalMatrix;

    ////////////////////////////////////////////////////////////////////////////
    //
    // Set coefficients in pixel
    //
    m_Coefficients.resize(m_NumberOfSHCoefficients);

    for(unsigned int i = 0; i < coefficientsMatrix.Rows(); i++)
    {
        m_Coefficients[i] = coefficientsMatrix(i,0);
    }

}

} // namespace btk
