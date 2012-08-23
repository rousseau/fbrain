/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 22/08/2012
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

#include "btkStreamlineTractographyAlgorithm.h"


// VTK includes
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyLine.h"

// Local includes
#include "btkGradientDirection.h"


namespace btk
{

StreamlineTractographyAlgorithm::StreamlineTractographyAlgorithm() : Superclass()
{
    // ----
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PropagateSeed(Self::PhysicalPoint point)
{
    //
    // Estimate Fiber
    //

    // Usefull constants
    float stepSize_2 = m_StepSize / 2.f;
    float stepSize_6 = m_StepSize / 6.f;

    bool stop = false;
    std::vector< Self::PhysicalPoint > points;
    points.push_back(point);

    btk::GradientDirection nextEigenVector = m_DiffusionModel->MeanDirectionsAt(point)[0];

    do
    {
        // This use the RK4 to compute the next point (Runge-Kutta, order 4)
        Self::PhysicalPoint lastPoint = points.back();

        btk::GradientDirection k1 = nextEigenVector;
        btk::GradientDirection k2 = m_DiffusionModel->MeanDirectionsAt(lastPoint + (k1*stepSize_2))[0];
        if(k2*k1 < 0) k2 *= -1;
        btk::GradientDirection k3 = m_DiffusionModel->MeanDirectionsAt(lastPoint + (k2*stepSize_2))[0];
        if(k3*k1 < 0) k3 *= -1;
        btk::GradientDirection k4 = m_DiffusionModel->MeanDirectionsAt(lastPoint + (k3*m_StepSize))[0];
        if(k4*k1 < 0) k4 *= -1;

        Self::PhysicalPoint nextPoint = lastPoint + (k1 + (k2*2.f) + (k3*2.f) + k4) * stepSize_6;//*/

        /*/ This use the RK1 (Euler) to compute the next point (Runge-Kutta, order 1, Euler method)
        btk::GradientDirection     k1 = nextEigenVector;
        Self::PhysicalPoint nextPoint = points.back() + (k1*m_StepSize);//*/

        // Check if the physical point is in the mask
        Self::MaskImage::IndexType maskIndex;
        m_Mask->TransformPhysicalPointToIndex(nextPoint, maskIndex);

        if(m_Mask->GetPixel(maskIndex) == 0)
        {
            stop = true;
        }
        else // m_Mask->GetPixel(maskIndex) != 0
        {
            nextEigenVector = m_DiffusionModel->MeanDirectionsAt(nextPoint)[0];

            // Check the consistency of the direction of the main eigenvector
            if(nextEigenVector*k1 < 0)
            {
                nextEigenVector *= -1;
            }

            // Add the new point
            points.push_back(nextPoint);
        }
    } while(!stop);


    //
    // Build graphical fiber
    //

    if(points.size() > 0)
    {
        m_CurrentFiber = vtkSmartPointer< vtkPolyData >::New();

        // Graphical representation structures
        vtkSmartPointer< vtkPoints >  vpoints = vtkSmartPointer< vtkPoints >::New();
        vtkSmartPointer< vtkCellArray > lines = vtkSmartPointer< vtkCellArray >::New();
        vtkSmartPointer< vtkPolyLine >   line = vtkSmartPointer< vtkPolyLine >::New();

        line->GetPointIds()->SetNumberOfIds(points.size());

        for(unsigned int i = 0; i < points.size(); i++)
        {
            vpoints->InsertNextPoint(-points[i][0], -points[i][1], points[i][2]);
            line->GetPointIds()->SetId(i,i);
        }

        lines->InsertNextCell(line);
        m_CurrentFiber->SetPoints(vpoints);
        m_CurrentFiber->SetLines(lines);
    }
}

} // namespace btk
