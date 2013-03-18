/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 04/01/2013
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


#include "btkParticleFilteringTractographyAlgorithm.h"


namespace btk
{

ParticleFilteringTractographyAlgorithm::ParticleFilteringTractographyAlgorithm()
{
    // ----
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

vtkSmartPointer< vtkPolyData > ParticleFilteringTractographyAlgorithm::PropagateSeed(Self::PhysicalPoint point)
{
    // Diffusion directions provided by the model at point
    std::vector< btk::GradientDirection > nextDirections = m_DiffusionModel->MeanDirectionsAt(point);

    // Graphical representation structures
    vtkSmartPointer< vtkPoints >        vpoints = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray >       lines = vtkSmartPointer< vtkCellArray >::New();
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

        this->PropagateSeed(points, nextDirection);

        //
        // Build graphical fiber
        //

        if(points.size() > 1)
        {
            vtkSmartPointer< vtkPolyLine > line = vtkSmartPointer< vtkPolyLine >::New();

            line->GetPointIds()->SetNumberOfIds(points.size());

            unsigned int nbOfPreviousPoints = vpoints->GetNumberOfPoints();

            vpoints->InsertNextPoint(-points[0][0], -points[0][1], points[0][2]);
            line->GetPointIds()->SetId(0,nbOfPreviousPoints);

            for(unsigned int i = 1; i < points.size(); i++)
            {
                // Insert cells (points and lines)
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

void ParticleFilteringTractographyAlgorithm::PropagateSeed(std::vector<Self::PhysicalPoint> &points, GradientDirection nextDirection)
{
    // TODO reflexions sur la meilleure manière d'aborder la programmation du problème
}

} // namespace btk
