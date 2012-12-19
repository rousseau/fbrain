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
#include "vtkPolyDataWriter.h"

// Local includes
#include "btkGradientDirection.h"


namespace btk
{

StreamlineTractographyAlgorithm::StreamlineTractographyAlgorithm() : m_StepSize(0.5), m_UseRungeKuttaOrder4(false), m_ThresholdAngle(M_PI/3.0f), Superclass()
{
    // ----
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

vtkSmartPointer< vtkPolyData > StreamlineTractographyAlgorithm::PropagateSeed(Self::PhysicalPoint point)
{
    // Diffusion directions provided by the model at point
    std::vector< btk::GradientDirection > nextDirections = m_DiffusionModel->MeanDirectionsAt(point);

    // Graphical representation structures
    vtkSmartPointer< vtkPoints >  vpoints = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray > lines = vtkSmartPointer< vtkCellArray >::New();

    vtkSmartPointer< vtkPolyData > currentFiber = vtkSmartPointer< vtkPolyData >::New();

    // Processing
    for(std::vector< btk::GradientDirection >::iterator it = nextDirections.begin(); it != nextDirections.end(); it++)
    {
        //
        // Estimate Fiber
        //

        std::vector< Self::PhysicalPoint > points;
        points.push_back(point);

        btk::GradientDirection nextDirection = *it;

        if(m_UseRungeKuttaOrder4)
        {
            Self::PropagateSeedRK4(points, nextDirection);
        }
        else // m_UseRungeKuttaOrder4 = false
        {
            Self::PropagateSeedRK1(points, nextDirection);
        }

        //
        // Build graphical fiber
        //

        if(points.size() > 1)
        {
            vtkSmartPointer< vtkPolyLine > line = vtkSmartPointer< vtkPolyLine >::New();

            line->GetPointIds()->SetNumberOfIds(points.size());

            unsigned int nbOfPreviousPoints = vpoints->GetNumberOfPoints();

            for(unsigned int i = 0; i < points.size(); i++)
            {
                vpoints->InsertNextPoint(-points[i][0], -points[i][1], points[i][2]);
                line->GetPointIds()->SetId(i,i+nbOfPreviousPoints);
            }

            lines->InsertNextCell(line);
        }
    } // for each direction

    // Build the whole trajectories
    currentFiber->SetPoints(vpoints);
    currentFiber->SetLines(lines);


    return currentFiber;
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PropagateSeedRK4(std::vector< Self::PhysicalPoint > &points, GradientDirection nextDirection)
{
    // Usefull constants
    float stepSize_2 = m_StepSize / 2.f;
    float stepSize_6 = m_StepSize / 6.f;

    bool stop = false;

    do
    {
        // This use the RK4 to compute the next point (Runge-Kutta, order 4)
        Self::PhysicalPoint lastPoint = points.back();

        btk::GradientDirection k1 = nextDirection;
        std::vector< btk::GradientDirection > directions = m_DiffusionModel->MeanDirectionsAt(lastPoint + (k1*stepSize_2), k1, m_ThresholdAngle);
        btk::GradientDirection k2 = Self::SelectClosestDirection(directions, k1);

        if(!k2.IsNull())
        {
                           directions = m_DiffusionModel->MeanDirectionsAt(lastPoint + (k2*stepSize_2), k2, m_ThresholdAngle);
            btk::GradientDirection k3 = Self::SelectClosestDirection(directions, k2);

            if(!k3.IsNull())
            {
                               directions = m_DiffusionModel->MeanDirectionsAt(lastPoint + (k3*m_StepSize), k3, m_ThresholdAngle);
                btk::GradientDirection k4 = Self::SelectClosestDirection(directions, k3);

                if(!k4.IsNull())
                {
                    Self::PhysicalPoint nextPoint = lastPoint + (k1 + (k2*2.f) + (k3*2.f) + k4) * stepSize_6;

                    // Check if the physical point is in the mask
                    Self::MaskImage::IndexType maskIndex;
                    m_Mask->TransformPhysicalPointToIndex(nextPoint, maskIndex);

                    if(m_Mask->GetPixel(maskIndex) == 0)
                    {
                        stop = true;
                    }
                    else // m_Mask->GetPixel(maskIndex) != 0
                    {
                        // Add the new point
                        points.push_back(nextPoint);

                        // Search next direction
                        std::vector< btk::GradientDirection > meanDirections = m_DiffusionModel->MeanDirectionsAt(nextPoint, k1, m_ThresholdAngle);
                        nextDirection = Self::SelectClosestDirection(meanDirections, k1);

                        if(nextDirection.IsNull())
                        {
                            stop = true;
                        }
                    }
                }
                else
                {
                    stop = true;
                }
            }
            else
            {
                stop = true;
            }
        }
        else
        {
            stop = true;
        }
    } while(!stop);
}

//----------------------------------------------------------------------------------------

void StreamlineTractographyAlgorithm::PropagateSeedRK1(std::vector< Self::PhysicalPoint > &points, GradientDirection nextDirection)
{
    bool stop = false;

    do
    {
        // This use the RK1 (Euler) to compute the next point (Runge-Kutta, order 1, Euler method)
        btk::GradientDirection     k1 = nextDirection;
        Self::PhysicalPoint nextPoint = points.back() + (k1*m_StepSize);

        // Check if the physical point is in the mask
        Self::MaskImage::IndexType maskIndex;
        m_Mask->TransformPhysicalPointToIndex(nextPoint, maskIndex);

        if(m_Mask->GetPixel(maskIndex) == 0)
        {
            stop = true;
        }
        else // m_Mask->GetPixel(maskIndex) != 0
        {
            // Add the new point
            points.push_back(nextPoint);

            // Search next direction
            std::vector< btk::GradientDirection > meanDirections = m_DiffusionModel->MeanDirectionsAt(nextPoint, k1, m_ThresholdAngle);
            nextDirection = Self::SelectClosestDirection(meanDirections, k1);

            if(nextDirection.IsNull())
            {
                stop = true;
            }
        }
    } while(!stop);
}

//----------------------------------------------------------------------------------------

btk::GradientDirection StreamlineTractographyAlgorithm::SelectClosestDirection(std::vector< btk::GradientDirection > &meanDirections, btk::GradientDirection &previousVector)
{
    unsigned int meanDirectionsSize = meanDirections.size();
    btk::GradientDirection nextDirection(0,0,0);

    if(meanDirectionsSize == 1)
    {
        nextDirection = meanDirections[0];
    }
    else if(meanDirectionsSize > 1)
    {
        // The next direction is choosen to be the closer to the previous one.
        float   minDotProduct = std::abs(meanDirections[0]*previousVector);
        unsigned int minIndex = 0;

        for(unsigned int i = 1; i < meanDirections.size(); i++)
        {
            float dotProduct = std::abs(meanDirections[i]*previousVector);

            if(dotProduct < minDotProduct)
            {
                minDotProduct = dotProduct;
                minIndex      = i;
            }
        } // for each mean direction

        nextDirection = meanDirections[minIndex];
    }

    return nextDirection;
}

} // namespace btk
