/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 27/05/2013
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

#include "btkSincPSF.h"

#include "btkMathFunctions.h"

namespace btk
{

SincPSF::SincPSF():Superclass::PSF()
{
//    m_Direction.set_size(3,3);

//    m_Center.set_size(3);
//    m_Center.fill(0.0);

//    m_Spacing.set_size(3);
//    m_Spacing.fill(1);

    m_PsfImage = ImageType::New();
}
//-------------------------------------------------------------------------------------------------
void SincPSF::Initialize()
{
    m_idir = m_Direction.get_column(0);
    m_jdir = m_Direction.get_column(1);
    m_kdir = m_Direction.get_column(2);
}
//-------------------------------------------------------------------------------------------------
SincPSF::OutputType
SincPSF::Evaluate(const InputType &position) const
{
    //NOTE: This fonction is deprecated
    vnl_vector<double> diff = position.GetVnlVector() - m_Center;
    PointType diffPoint;
    double x,y,z;

    //Dot product between image direction and point vector (in PSF space)
    double icoor = dot_product(diff,m_idir);
    double jcoor = dot_product(diff,m_jdir);
    double kcoor = dot_product(diff,m_kdir);

    x = diffPoint[0] = icoor;
    y = diffPoint[1] = jcoor;
    z = diffPoint[2] = kcoor;

    double value = 0.0;

    value = sin((x*x) + (y*y) + (z*z)) / ((x*x) + (y*y) + (z*z));

    return (OutputType)value;

}
//-------------------------------------------------------------------------------------------------
void SincPSF::ConstructImage()
{
    ImageType::RegionType region;

    //ImageType::SizeType size;

    ImageType::IndexType index;

    //size = _size;

    index[0] = 0;
    index[1] = 0;
    index[2] = 0;

    region.SetSize(m_Size);
    region.SetIndex(index);

    ImageType::SpacingType spc;
    spc[0] = m_Spacing[0];
    spc[1] = m_Spacing[1];
    spc[2] = m_Spacing[2];

    m_PsfImage->SetRegions(region);
    m_PsfImage->SetSpacing(spc);
    m_PsfImage->Allocate();
    m_PsfImage->FillBuffer(0.0);
    double sum = 0;

    ImageType::IndexType hrIndex;

    itkIteratorWithIndex itPSF(m_PsfImage,m_PsfImage->GetLargestPossibleRegion());
    itkContinuousIndex hrIndexCenter;
    hrIndexCenter[0] = (m_Size[0]-1)/2.0;
    hrIndexCenter[1] = (m_Size[1]-1)/2.0;
    hrIndexCenter[2] = (m_Size[2]-1)/2.0;
    for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
    {
        hrIndex = itPSF.GetIndex();
        float x = hrIndex[0]- hrIndexCenter[0];
        float y = hrIndex[1]- hrIndexCenter[1];
        float z = hrIndex[2]- hrIndexCenter[2];
        // old version (like the definition for a sinc in multidimension)
        //float value = sin((x*x) + (y*y) + (z*z)) / ((x*x) + (y*y) + (z*z));
        //new version
        float value = MathFunctions::Sinc(x) * MathFunctions::Sinc(y)
                      * MathFunctions::Sinc(z);

        //std::cout<<"value : "<<value<<std::endl;
        itPSF.Set(value);

        sum += itPSF.Get();
    }


    if(sum>0)
    {
        for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
        {
            if(itPSF.Get() > 0)
            {
                itPSF.Set( itPSF.Get() / sum );
            }
        }
    }




}
//-------------------------------------------------------------------------------------------------
void SincPSF::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os,indent);

    os << indent << "Direction: " << m_Direction << std::endl;
    os << indent << "Center: " << m_Center << std::endl;
    os << indent << "Spacing: " << m_Spacing << std::endl;
}
}
