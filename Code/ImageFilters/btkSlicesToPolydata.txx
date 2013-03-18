/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 08/03/2013
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

#ifndef BTKSLICESTOPOLYDATA_TXX
#define BTKSLICESTOPOLYDATA_TXX

#include "btkSlicesToPolydata.h"

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPointData.h"
#include "vtkPolyLine.h"
#include "vtkCommand.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataWriter.h"

namespace btk
{
template<typename TImage>
SlicesToPolyData< TImage >::SlicesToPolyData()
{
    m_Output = NULL;
    m_Transform = NULL;
    m_Input = NULL;
}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
SlicesToPolyData< TImage >::~SlicesToPolyData()
{
    m_Output = NULL;
    m_Transform = NULL;
    m_Input = NULL;
}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void
SlicesToPolyData< TImage >::Initialize()
{
    m_Output = vtkSmartPointer< vtkPolyData >::New();
    if(m_Input.IsNull())
    {
        btkException("Missing input image !");
    }
    if(m_Transform == NULL)
    {
       btkException("Missing transform !");
    }
}
//-------------------------------------------------------------------------------------------------
template<typename TImage>
void
SlicesToPolyData< TImage >::Update()
{
    this->Initialize();

    typename ImageType::SizeType size = m_Input->GetLargestPossibleRegion().GetSize();
    vtkSmartPointer< vtkPoints > vtkP = vtkSmartPointer< vtkPoints > ::New();
    vtkSmartPointer< vtkCellArray > Lines = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkPolyLine > line = vtkSmartPointer< vtkPolyLine >::New();

    for(unsigned int z = 0; z < size[2]; z++)
    {
        typename ImageType::IndexType corner1, corner2, corner3, corner4;

        corner1[0] = 0; corner1[1] = 0; corner1[2] = z;
        corner2[0] = 0; corner2[1] = size[1]-1; corner2[2] = z;
        corner3[0] = size[0]-1; corner3[1] = size[1] -1 ; corner3[2] = z;
        corner4[0] = size[0]-1; corner4[1] = 0; corner4[2] = z;

        typename ImageType::PointType point1, point2, point3, point4;
        typename ImageType::PointType Tpoint1, Tpoint2, Tpoint3, Tpoint4;

        m_Input->TransformIndexToPhysicalPoint(corner1,point1);
        m_Input->TransformIndexToPhysicalPoint(corner2,point2);
        m_Input->TransformIndexToPhysicalPoint(corner3,point3);
        m_Input->TransformIndexToPhysicalPoint(corner4,point4);

        Tpoint1 = m_Transform->TransformPoint(point1);
        Tpoint2 = m_Transform->TransformPoint(point2);
        Tpoint3 = m_Transform->TransformPoint(point3);
        Tpoint4 = m_Transform->TransformPoint(point4);



        line->GetPointIds()->SetNumberOfIds(5);
        line->GetPointIds()->SetId(0, vtkP->InsertNextPoint(Tpoint1[0],
                                   Tpoint1[1],
                Tpoint1[2]));

        line->GetPointIds()->SetId(1, vtkP->InsertNextPoint(Tpoint2[0],
                                   Tpoint2[1],
                Tpoint2[2]));

        line->GetPointIds()->SetId(2, vtkP->InsertNextPoint(Tpoint3[0],
                                   Tpoint3[1],
                Tpoint3[2]));

        line->GetPointIds()->SetId(3, vtkP->InsertNextPoint(Tpoint4[0],
                                   Tpoint4[1],
                Tpoint4[2]));
        line->GetPointIds()->SetId(4, vtkP->InsertNextPoint(Tpoint1[0],
                                   Tpoint1[1],
                Tpoint1[2]));

        Lines->InsertNextCell(line);



    }

    m_Output->SetPoints(vtkP);
    m_Output->SetLines(Lines);


}
}
#endif
