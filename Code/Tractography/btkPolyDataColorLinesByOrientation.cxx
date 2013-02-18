/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 03/01/2013
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

#include "btkPolyDataColorLinesByOrientation.h"


// VTK includes
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPolyData.h"
#include "vtkDataObject.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPolyLine.h"
#include "vtkCommand.h"

// ITK includes
#include "itkVector.h"

// Local includes
#include "btkMacro.h"


namespace btk
{

vtkCxxRevisionMacro(PolyDataColorLinesByOrientation, "$Revision: 1.0 $");
vtkStandardNewMacro(PolyDataColorLinesByOrientation);

//----------------------------------------------------------------------------------------

PolyDataColorLinesByOrientation::PolyDataColorLinesByOrientation() : m_ColorOrientation(COLOR_MEAN_ORIENTATION)
{
    // ----
}

//----------------------------------------------------------------------------------------

PolyDataColorLinesByOrientation::~PolyDataColorLinesByOrientation()
{
    // ----
}

//----------------------------------------------------------------------------------------

int PolyDataColorLinesByOrientation::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    // Get the information objects
    vtkInformation  *inputInfo = inputVector[0]->GetInformationObject(0);
    vtkInformation *outputInfo = outputVector->GetInformationObject(0);

    // Get the input and output
    vtkPolyData *input = vtkPolyData::SafeDownCast(inputInfo->Get(vtkDataObject::DATA_OBJECT()));
    vtkPolyData *output = vtkPolyData::SafeDownCast(outputInfo->Get(vtkDataObject::DATA_OBJECT()));


    // Define structures
    vtkSmartPointer< vtkCellArray > inputLines = input->GetLines();

    vtkSmartPointer< vtkPoints >     outputPoints = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray >   outputLines = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkFloatArray > outputColors = vtkSmartPointer< vtkFloatArray >::New();

    switch(m_ColorOrientation)
    {
        case COLOR_MEAN_ORIENTATION:
            this->ColorByMeanOrientation(input, output, inputLines, outputPoints, outputLines, outputColors);
            break;

        case COLOR_LOCAL_ORIENTATION:
            this->ColorByLocalOrientation(input, output, inputLines, outputPoints, outputLines, outputColors);
            break;

        default:
            btkException("Unknown orientation type for fiber tracts coloration !");
    }

    // Update progress bar
    this->SetProgress(1.0);
    this->InvokeEvent(vtkCommand::ProgressEvent);


    return 1;
}

//----------------------------------------------------------------------------------------

inline void PolyDataColorLinesByOrientation::ColorByMeanOrientation(vtkSmartPointer< vtkPolyData > input, vtkSmartPointer< vtkPolyData > output, vtkSmartPointer< vtkCellArray > inputLines, vtkSmartPointer< vtkPoints > outputPoints, vtkSmartPointer< vtkCellArray > outputLines, vtkSmartPointer< vtkFloatArray > outputColors)
{
    // Initialize the colors
    outputColors->SetNumberOfComponents(1);
    outputColors->SetName("MeanOrientation");


    // Initialize the progress bar
    this->SetProgress(0.0);
    this->InvokeEvent(vtkCommand::ProgressEvent);
    unsigned int numberOfTracts = inputLines->GetNumberOfCells();
    double         progressStep = 1.0 / numberOfTracts;


    vtkIdType numberOfPoints, *pointIds;

    // Compute colors for each fiber
    while(inputLines->GetNextCell(numberOfPoints, pointIds) != 0)
    {
        itk::Vector< float,3 > displacement;
        displacement.Fill(0);

        // Get first point
        double firstPoint[3];
        input->GetPoint(pointIds[0], firstPoint);

        // Initialize current line and set first point
        vtkSmartPointer< vtkPolyLine > line = vtkSmartPointer< vtkPolyLine >::New();
        line->GetPointIds()->SetNumberOfIds(numberOfPoints);
        line->GetPointIds()->SetId(0, outputPoints->InsertNextPoint(firstPoint));

        // Compute mean direction based on each points
        for(unsigned int i = 1; i < numberOfPoints; i++)
        {
            // Get points' coordinates
            double previousPoint[3], currentPoint[3];
            input->GetPoint(pointIds[i-1], previousPoint);
            input->GetPoint(pointIds[i], currentPoint);

            // Add current point to output
            line->GetPointIds()->SetId(i, outputPoints->InsertNextPoint(currentPoint));

            // Add current displacement vector
            displacement[0] += std::abs(-currentPoint[0]+previousPoint[0]); displacement[1] += std::abs(-currentPoint[1]+previousPoint[1]); displacement[2] += std::abs(currentPoint[2]-previousPoint[2]);
        } // for each point

        // Add line to polydata
        outputLines->InsertNextCell(line);

        // Compute the mean orientation and the associated color
        displacement[0] /= numberOfPoints; displacement[1] /= numberOfPoints; displacement[2] /= numberOfPoints;

        // Compute rgb color index
        float color[1] = { 0 };
        displacement.Normalize();
        displacement[0] *= 255; displacement[1] *= 255; displacement[2] *= 255;
        color[0] = btk::RGBtoIndex(displacement[0], displacement[1], displacement[2]);

        // Add one scalar per point
        for(unsigned int i = 0; i < numberOfPoints; i++)
        {
            outputColors->InsertNextTupleValue(color);
        }

        // Update progress
        this->SetProgress(this->GetProgress() + progressStep);
        this->InvokeEvent(vtkCommand::ProgressEvent);
    } // for each fiber

    // Set output
    output->SetPoints(outputPoints);
    output->SetLines(outputLines);
    output->GetPointData()->SetScalars(outputColors);
}

//----------------------------------------------------------------------------------------

inline void PolyDataColorLinesByOrientation::ColorByLocalOrientation(vtkSmartPointer< vtkPolyData > input, vtkSmartPointer< vtkPolyData > output, vtkSmartPointer< vtkCellArray > inputLines, vtkSmartPointer< vtkPoints > outputPoints, vtkSmartPointer< vtkCellArray > outputLines, vtkSmartPointer< vtkFloatArray > outputColors)
{
    // Initialize the colors
    outputColors->SetNumberOfComponents(1);
    outputColors->SetName("LocalOrientation");


    // Initialize the progress bar
    this->SetProgress(0.0);
    this->InvokeEvent(vtkCommand::ProgressEvent);
    unsigned int numberOfTracts = inputLines->GetNumberOfCells();
    double         progressStep = 1.0 / numberOfTracts;


    vtkIdType numberOfPoints, *pointIds;

    // Compute colors for each fiber
    while(inputLines->GetNextCell(numberOfPoints, pointIds) != 0)
    {
        itk::Vector< float,3 > displacement;
        displacement.Fill(0);

        float color[1] = { 0 };

        // Get first point
        double firstPoint[3];
        input->GetPoint(pointIds[0], firstPoint);

        // Initialize current line and set first point
        vtkSmartPointer< vtkPolyLine > line = vtkSmartPointer< vtkPolyLine >::New();
        line->GetPointIds()->SetNumberOfIds(numberOfPoints);
        line->GetPointIds()->SetId(0, outputPoints->InsertNextPoint(firstPoint));

        // Compute mean direction based on each points
        for(unsigned int i = 1; i < numberOfPoints; i++)
        {
            // Get points' coordinates
            double previousPoint[3], currentPoint[3];
            input->GetPoint(pointIds[i-1], previousPoint);
            input->GetPoint(pointIds[i], currentPoint);

            // Add current point to output
            line->GetPointIds()->SetId(i, outputPoints->InsertNextPoint(currentPoint));

            // Add current displacement vector
            displacement[0] = std::abs(-currentPoint[0]+previousPoint[0]); displacement[1] = std::abs(-currentPoint[1]+previousPoint[1]); displacement[2] = std::abs(currentPoint[2]-previousPoint[2]);

            // Compute rgb color index

            displacement.Normalize();
            displacement[0] *= 255; displacement[1] *= 255; displacement[2] *= 255;
            color[0] = btk::RGBtoIndex(displacement[0], displacement[1], displacement[2]);

            // Add color index to current previous point
            outputColors->InsertNextTupleValue(color);
        } // for each point

        // Add color of before last point to the last point
        outputColors->InsertNextTupleValue(color);

        // Add line to polydata
        outputLines->InsertNextCell(line);

        // Update progress
        this->SetProgress(this->GetProgress() + progressStep);
        this->InvokeEvent(vtkCommand::ProgressEvent);
    } // for each fiber

    // Set output
    output->SetPoints(outputPoints);
    output->SetLines(outputLines);
    output->GetPointData()->SetScalars(outputColors);
}

//----------------------------------------------------------------------------------------

int PolyDataColorLinesByOrientation::FillInputPortInformation(int port, vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
}

//----------------------------------------------------------------------------------------

void PolyDataColorLinesByOrientation::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os,indent);
}

} // namespace btk
