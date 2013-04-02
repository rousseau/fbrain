/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 05/07/2012
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

#include "btkDiffusionSequence.h"


// VNL includes
#include "vnl/vnl_inverse.h"


namespace btk
{

DiffusionSequence::DiffusionSequence()
{
    // ----
}

//----------------------------------------------------------------------------------------

DiffusionSequence::~DiffusionSequence()
{
    // ----
}

//----------------------------------------------------------------------------------------

void DiffusionSequence::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "Number of gradients: " << m_GradientTable.size() << std::endl;
    os << indent << "Gradient table: " << std::endl;
    for(unsigned int i = 0; i < m_GradientTable.size(); i++)
    {
        os << indent << indent << m_GradientTable[i] << std::endl;
    }

    os << indent << "B-Value: " << m_BValues[1] << std::endl;
}

//----------------------------------------------------------------------------------------

void DiffusionSequence::ConvertGradientTableToPhysicalCoordinates()
{
    for(unsigned int i = 0; i < m_GradientTable.size(); i++)
    {
        vnl_matrix< double > direction = this->GetDirection().GetVnlMatrix().extract(3,3,0,0);
        vnl_vector< double >  gradient = m_GradientTable[i].GetVnlVector();

        gradient = direction * gradient;
        gradient.normalize();
        m_GradientTable[i] = GradientDirection(gradient[0], gradient[1], gradient[2]);
    }
}

//----------------------------------------------------------------------------------------

void DiffusionSequence::ConvertGradientTableToImageCoordinates()
{
    for(unsigned int i = 0; i < m_GradientTable.size(); i++)
    {
        vnl_matrix< double > directionInverse = vnl_inverse(this->GetDirection().GetVnlMatrix().extract(3,3,0,0));
        vnl_vector< double >         gradient = m_GradientTable[i].GetVnlVector();

        gradient = directionInverse * gradient;
        gradient.normalize();
        m_GradientTable[i] = GradientDirection(gradient[0], gradient[1], gradient[2]);
    }
}

} // namespace btk
