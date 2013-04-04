/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 07/03/2013
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

#include "btkImportanceDensity.h"


// STL includes
#include "cfloat"

// Local includes
#include "btkVonMisesFisherProbabilityDensity.h"


namespace btk
{

ImportanceDensity::ImportanceDensity() : m_Model(NULL), m_Signal(NULL), m_AngleThreshold(M_PI/6.0)
{
    // ----
}

//----------------------------------------------------------------------------------------

ImportanceDensity::ImportanceDensity(DiffusionModel::Pointer model, DiffusionSignal::Pointer signal, double angleThreshold) : m_Model(model), m_Signal(signal), m_AngleThreshold(angleThreshold/2.0)
{
    // ----
}

//----------------------------------------------------------------------------------------

GradientDirection ImportanceDensity::GetMeanDirection(PhysicalPoint pk, GradientDirection vkm1)
{
    std::vector< GradientDirection > maxima = m_Model->MeanDirectionsAt(pk, vkm1, m_AngleThreshold);

    GradientDirection maximum;

    if(maxima.size() > 1)
    {
        // Compute
        double sum = 0.0;
        std::vector< double > values;

        for(unsigned int i = 0; i < maxima.size(); i++)
        {
            values.push_back(m_Model->ModelAt(pk, vkm1));
            sum += values.back();
        }

        // Simulate x~U(0,1)
        double x = (static_cast< double >(std::rand()) / static_cast< double >(RAND_MAX)) * sum;

        // compare to intervals and choose the mean direction
        bool      found = false;
        unsigned int  i = 0;
        double interval = 0;

        while(!found && i < values.size())
        {
            interval += values[i];

            if(x < interval)
                found = true;
            else
                i++;
        }

        maximum = maxima[i];
    }
    else if(maxima.size() == 1)
    {
        maximum = maxima.front();
    }
//    else // maxima.size() == 0
//    {
//        maximum = GradientDirection();
//    }

    return maximum;
}

//----------------------------------------------------------------------------------------

GradientDirection ImportanceDensity::Simulate(GradientDirection meanDirection, double concentration)
{
    // Simulate the vMF distribution with mean (0,0,1)
    double        y = static_cast< double >(std::rand()) / static_cast< double >(RAND_MAX);
    double        w = (1.0 / concentration) * std::log( std::exp(-concentration) + concentration * ((2.0/concentration) * (std::exp(concentration) - std::exp(-concentration))/2.0) * y);
    double    theta = ( static_cast< double >(std::rand()) / static_cast< double >(RAND_MAX)) * m_2pi;
    double constant = std::sqrt(1-w*w);

    GradientDirection x(constant*std::cos(theta), constant*std::sin(theta), w);

    // Rotate if necessary
    if(meanDirection[0] != 0 || meanDirection[1] != 0 || meanDirection[2] != 1) // modal direction is not (0,0,1)
    {
        SphericalDirection mu = meanDirection.GetSphericalDirection();

        // Rotations from (0,0,1) to current modal direction
        double cosTheta = std::cos(mu[0]);
        double sinTheta = std::sin(mu[0]);
        double cosPhi   = std::cos(mu[1]);
        double sinPhi   = std::sin(mu[1]);

        // Ry(theta)
        // | cos(theta)  0  sin(theta)  |   |x|
        // |     0       1       0      | . |y|
        // |-sin(theta)  0   cos(theta) |   |z|
        double Ryx = cosTheta * x[0] + sinTheta * x[2];
        double Ryy = x[1];
        double Ryz = -sinTheta * x[0] + cosTheta * x[2];

        // Rz(phi)
        // | cos(phi)  -sin(phi)  0 |   |x|
        // | sin(phi)   cos(phi)  0 | . |y|
        // |    0          0      1 |   |z|
        x[0] = cosPhi * Ryx - sinPhi * Ryy;
        x[1] = sinPhi * Ryx + cosPhi * Ryy;
        x[2] = Ryz;

        x.UpdateSphericalCoordinates();
    }

    return x;
}

//----------------------------------------------------------------------------------------

double ImportanceDensity::EstimateConcentrationParameter(GradientDirection mu, PhysicalPoint pk)
{
    std::vector< float >        model = m_Model->SignalAt(pk);
    DiffusionSignal::PixelType signal = m_Signal->SignalAt(pk);

    std::vector< GradientDirection > directions = m_Signal->GetGradientTable();

    double sumSquare = 0.0;
    double min = DBL_MAX, max = DBL_MIN;
    unsigned int numberOfUsedDirections = 0;

    for(unsigned int i = 0; i < directions.size(); i++)
    {
        GradientDirection g = directions[i];

        double scalaProduct = g * mu;
        double        angle = std::acos( scalaProduct / std::sqrt(g.GetSquaredNorm() * mu.GetSquaredNorm()));

        if(angle <= m_AngleThreshold)
        {
            double diff = model[i] - signal[i];
            sumSquare  += diff*diff;

            if(min > model[i])
                min = model[i];

            if(min > signal[i])
                min = signal[i];

            if(max < model[i])
                max = model[i];

            if(max < signal[i])
                max = signal[i];

            numberOfUsedDirections++;
        }
    }

    double NRMSE = std::sqrt( sumSquare / numberOfUsedDirections ) / (max - min);
    double kappa = (NRMSE > 0.3) ? 10 : 435.1455649013809 * std::exp(-25.74351897609373 * NRMSE + 44.9437993453943 * NRMSE*NRMSE);

    return (kappa <= 200) ? kappa : 200;
}

//----------------------------------------------------------------------------------------

double ImportanceDensity::Evaluate(GradientDirection vk, GradientDirection vkm1, double concentration)
{
    double normalizationCoefficient = std::log(concentration) - ( std::log(m_2pi) + std::log(std::exp(concentration) - std::exp(-concentration)));

    return normalizationCoefficient + (vkm1 * vk * concentration);
}

} // namespace btk
