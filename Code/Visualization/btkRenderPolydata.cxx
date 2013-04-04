/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date:08/03/2013
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

#include "btkRenderPolydata.h"


#include "vtkWindowToImageFilter.h"

#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"




namespace btk
{

RenderPolyData::RenderPolyData():m_NumberOfPolyData(0)
{
    this->m_RenderWindow = vtkSmartPointer< vtkRenderWindow >::New();
    this->m_Renderer = vtkSmartPointer< vtkRenderer >::New();
    this->m_Iren = vtkSmartPointer< vtkRenderWindowInteractor >::New();
    this->m_Interactor = vtkSmartPointer< vtkInteractorStyleTrackballCamera >::New();
    m_RenderWindow->SetSize(400,400);
}

//-------------------------------------------------------------------------------------------------
RenderPolyData::~RenderPolyData()
{

}
//-------------------------------------------------------------------------------------------------
void RenderPolyData::SetNumberOfPolyData(unsigned int _n)
{
    m_Inputs.resize(_n);
    m_Colors.resize(_n);
    this->m_NumberOfPolyData = _n;
}
//-------------------------------------------------------------------------------------------------
void RenderPolyData::
SetNthPolyData(unsigned int _n, vtkSmartPointer< vtkPolyData > _PolyData )
{
    if(_n > m_NumberOfPolyData)
    {
        //btkException("The "<<_n<<"th PolyData does not exist, only "
        //      <<this->GetNumberOfPolyData()<<" available ! check SetNumberOfPolyDatas() !");
    }

    m_Inputs[_n] = _PolyData;
}
//-------------------------------------------------------------------------------------------------
void RenderPolyData::
SetNthPolyDataColor(unsigned int _n, std::vector<double> _color )
{
    if(_n > m_NumberOfPolyData)
    {
        //btkException("The "<<_n<<"th PolyData does not exist, only "
        //<<this->GetNumberOfPolyData()<<" available ! check SetNumberOfPolyDatas() !");
    }

    m_Colors[_n] = _color;
    //std::cout<<"color :"<<_n<<" : "<<_color<<std::endl;
}
//-------------------------------------------------------------------------------------------------
void RenderPolyData::Render()
{
    std::cout<<"Render"<<std::endl;
    std::cout<<m_NumberOfPolyData<<std::endl;
    for(unsigned int i=0; i< m_NumberOfPolyData; i++)
    {
        vtkSmartPointer<vtkPolyDataMapper> mapper_data =
                vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
        mapper_data->SetInput(m_Inputs[i]);
#else
        mapper_data->SetInputData(m_Inputs[i]);
#endif

        vtkSmartPointer<vtkActor> actor_data =
                vtkSmartPointer<vtkActor>::New();
        actor_data->SetMapper(mapper_data);

        double color[3] = {m_Colors[i][0],m_Colors[i][1],m_Colors[i][2]};

        actor_data->GetProperty()->SetColor(color);
        m_Renderer->AddActor(actor_data);

    }

    m_RenderWindow->AddRenderer(m_Renderer);
    m_Iren->SetRenderWindow(m_RenderWindow);
    m_Iren->SetInteractorStyle(m_Interactor);

    m_RenderWindow->Render();
    m_Iren->Start();

}
}//end namespace
