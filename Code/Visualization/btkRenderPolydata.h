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

#ifndef BTKRENDERPOLYDATA_H
#define BTKRENDERPOLYDATA_H

/* ITK */
#include "itkMacro.h"
#include "itkObjectFactory.h"
#include "itkObject.h"
#include "itkSmartPointer.h"

/* BTK */
#include "btkMacro.h"

/* VTK */
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"

/* Others */
#include "vector"

namespace btk
{
/**
 * @class RenderPolyData
 * @brief Class for render one or several vtk::PolyDatas
 * First we Set the number of PolyData to render with SetNumberOfPolyData.
 * Next we Set the Polydata with SetNthPolyData (inner a for loop for example) and we add a color with SetNthPolyDataColor
 * At the end why call Render().
 *
 * Note that if we recall Render() it will be update the render window (if a polydata change for example).
 *
 * @author Marc Schweitzer
 * @ingroup Visualization
 */
class RenderPolyData : public itk::Object
{
    public :
        /** typedefs */
        typedef RenderPolyData   Self;
        typedef itk::Object Superclass;
        typedef itk::SmartPointer< Self > Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Render Method */
        virtual void Render();
        /** Set the nth PolyData (_n need to be <= m_NumberOfPolyData) */
        virtual void SetNthPolyData(unsigned int _n, vtkSmartPointer< vtkPolyData > _PolyData );
        /** Set the color of the nth PolyData (_n need to be <= m_NumberOfPolyData) */
        virtual void SetNthPolyDataColor(unsigned int _n, std::vector<double>);
        /** Set the width of the nth PolyData (_n need to be <= m_NumberOfPolyData) */
        virtual void SetNthPolyDataLineWidth(unsigned int _n, float _w);
        /** Set the size of the nth PolyData (_n need to be <= m_NumberOfPolyData) */
        virtual void SetNthPolyDataPointSize(unsigned int _n, float _s);
        /** Set/Get Number of PolyData to render */
        virtual void SetNumberOfPolyData(unsigned int _n);
        btkGetMacro(NumberOfPolyData,unsigned int);

    protected:
        /** Constructor */
        RenderPolyData();
        /** Destructor */
        virtual ~RenderPolyData();

    private :

        std::vector< vtkSmartPointer< vtkPolyData > > m_Inputs; /** vector of inputs polydatas */
        unsigned int m_NumberOfPolyData; /** Number of PolyData to render */
        std::vector< std::vector<double> > m_Colors; /** Corresponding color of polydatas */
        std::vector< float >    m_Widths; /** Widths of the PolyDatas (exemples if a line) */
        std::vector< float >    m_PointSizes; /** size of points */

        vtkSmartPointer< vtkRenderWindow > m_RenderWindow; /** render window */
        vtkSmartPointer< vtkInteractorStyleTrackballCamera > m_Interactor; /** Interactor */
        vtkSmartPointer< vtkRenderWindowInteractor > m_Iren; /** render window interactor */
        vtkSmartPointer< vtkRenderer > m_Renderer; /** renderer */

};

} //end namespace
#endif // BTKRENDERPOLYDATA_H
