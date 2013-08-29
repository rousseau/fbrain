/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
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

#include "btkDiffusionSlice.h"
#include "btkImageHelper.h"

namespace btk
{

DiffusionSlice::DiffusionSlice()
{
    m_BValue = 0;
    m_OutlierStatus = false;

}

//----------------------------------------------------------------------------------------

DiffusionSlice::~DiffusionSlice()
{

    // ----
}


//----------------------------------------------------------------------------------------

void DiffusionSlice::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    //Superclass::PrintSelf(os,indent);
    os<<indent<< "Origin: "<<this->m_Origin << std::endl;
    os<<indent<< "Spacing: "<<this->m_Spacing << std::endl;
    os<<indent<< "Direction: "<<this->m_Direction << std::endl;
    os<<indent<< "Region: "<<this->GetLargestPossibleRegion()<< std::endl;

    os << indent << "Gradient direction: "<< m_GradientDirection<<std::endl;
    os << indent << "B-Value: " << m_BValue << std::endl;
    os << indent << "Outlier status: " << m_OutlierStatus << std::endl;


}

void DiffusionSlice::SetImage(Superclass::ConstPointer image)
{

    this->Initialize();
    this->Graft(image);

    // Correct physical point use as origin;
    DiffusionSlice::RegionType region= this->GetLargestPossibleRegion();
    DiffusionSlice::IndexType originIndex = region.GetIndex();

    DiffusionSlice::PointType originPoint;
    originPoint[0] = this->m_Spacing[0]*originIndex[0];
    originPoint[1] = this->m_Spacing[1]*originIndex[1];
    originPoint[2] = this->m_Spacing[2]*originIndex[2];


    originPoint =this->m_Direction*originPoint; // rotate in the physical space
    this->m_Origin+=originPoint.GetVectorFromOrigin(); // set offset

    region.SetIndex(2,0); // change index of image
    this->SetRegions(region);

}

void DiffusionSlice::Transform(TransformPointer transform)
{
    TransformPointer inverse = TransformType::New(); // itk transform norm: output-to-input transform
    transform->GetInverse(inverse);

    TransformType::OffsetType offset = inverse -> GetOffset(); // Ge Offset (translation)
    TransformType::MatrixType matrix = inverse -> GetMatrix(); // Get Rotation Matrix

    PointType offsetPoint;
    offsetPoint[0]=offset[0];
    offsetPoint[1]=offset[1];
    offsetPoint[2]=offset[2];

    // Transform origin
    this->m_Origin =   matrix *m_Origin + offset;

    // Basis change and transform direction
    this->m_Direction = m_Direction*(m_InverseDirection * matrix * m_Direction);
    this->m_InverseDirection = m_Direction.GetInverse();
    this->ComputeIndexToPhysicalPointMatrices();

    // Rotate gradient direction
    this->m_GradientDirection = matrix*this->m_GradientDirection;
    this->m_GradientDirection.Normalize();

}


} //end namespace btk
