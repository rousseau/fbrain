/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 
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

#include "btkGaussianPSF.h"

namespace btk
{

GaussianPSF::GaussianPSF()
{
    m_Direction.set_size(3,3);

    m_Center.set_size(3);
    m_Center.fill(0.0);

    m_Spacing.set_size(3);
    m_Spacing.fill(1);

    // Gaussian setup (ITK object)
//    m_Gaussian = GaussianFunctionType::New();
//    m_Gaussian -> SetNormalized(alse);//false

    ArrayType mean;
    mean[0] = 0; mean[1] = 0; mean[2] = 0;
    //m_Gaussian -> SetMean( mean );

    ArrayType sigma;

    //Compute sigma of the Gaussian PSF
//    m_Sigma[0] = sqrt(m_Spacing[0]*m_Spacing[0]/(8*log(2)));
//    m_Sigma[1] = sqrt(m_Spacing[1]*m_Spacing[1]/(8*log(2)));
//    m_Sigma[2] = sqrt(m_Spacing[2]*m_Spacing[2]/(8*log(2)));

    m_Sigma[0] = m_Spacing[0] / 2.3548;
    m_Sigma[1] = m_Spacing[1] / 2.3548;
    m_Sigma[2] = m_Spacing[2] / 2.3548;

//    sigma[0] = 0.5*m_Spacing[0];
//    sigma[1] = 0.5*m_Spacing[1];
//    sigma[2] = 0.5*m_Spacing[2];

    //m_Gaussian -> SetSigma( m_Sigma );
}
//-------------------------------------------------------------------------------------------------
GaussianPSF::OutputType
GaussianPSF::Evaluate(const InputType & position) const
{

    vnl_vector<double> diff = position.GetVnlVector() - m_Center;
    PointType diffPoint;
    double x , y, z;

    //Dot product between image direction and point vector (in PSF space)
    double icoor = dot_product(diff,m_idir);
    double jcoor = dot_product(diff,m_jdir);
    double kcoor = dot_product(diff,m_kdir);

    x = diffPoint[0] = icoor;
    y = diffPoint[1] = jcoor;
    z = diffPoint[2] = kcoor;

    //TODO: Do an oversampling over Point !



    double value = 0.0;

    value = (x*x)/(2*m_Sigma[0]*m_Sigma[0]) + (y*y)/(2*m_Sigma[1]*m_Sigma[1]) + (z*z)/(2*m_Sigma[2]*m_Sigma[2]);

    value = exp(-value);

   // std::cout<<value<<std::endl;
    //value = m_Gaussian->Evaluate(diffPoint);

    return(OutputType)value;
}
//-------------------------------------------------------------------------------------------------
void GaussianPSF::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Direction: " << m_Direction << std::endl;
    os << indent << "Center: " << m_Center << std::endl;
    os << indent << "Spacing: " << m_Spacing << std::endl;
}

}// end namespace
