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

#include "btkBoxCarPSF.h"

namespace btk
{

BoxCarPSF::BoxCarPSF()
{

}

//-------------------------------------------------------------------------------------------------

BoxCarPSF::OutputType
BoxCarPSF::Evaluate(const InputType &position) const
{
    vnl_vector<double> diff = position.GetVnlVector() - m_Center;
    PointType diffPoint;

    //Dot product between image direction and point vector (in PSF space)
    double icoor = dot_product(diff,m_idir);
    double jcoor = dot_product(diff,m_jdir);
    double kcoor = dot_product(diff,m_kdir);

    diffPoint[0] = icoor;
    diffPoint[1] = jcoor;
    diffPoint[2] = kcoor;

    double value = 0.0;

    if(( fabs(icoor) <= 0.5 * m_Spacing[0] ) &&
       ( fabs(jcoor) <= 0.5 * m_Spacing[1] ) &&
       ( fabs(kcoor) <= 0.5 * m_Spacing[2]))
    {
        value = 1.0;
    }

    return (OutputType)value;

}
//-------------------------------------------------------------------------------------------------
void BoxCarPSF::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Direction: " << m_Direction << std::endl;
    os << indent << "Center: " << m_Center << std::endl;
    os << indent << "Spacing: " << m_Spacing << std::endl;
}
//-------------------------------------------------------------------------------------------------

}
