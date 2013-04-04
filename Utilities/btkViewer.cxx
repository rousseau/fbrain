/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 02/04/2013
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

// STL includes
#include "vector"
#include "algorithm"
#include "cmath"
#include "cfloat"
#include "ctime"

// TCLAP includes
#include "tclap/CmdLine.h"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkActor.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolygon.h"
#include "vtkColorTransferFunction.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPointData.h"
#include "vtkSphereSource.h"
#include "vtkProperty.h"
#include "vtkLineSource.h"
#include "vtkTubeFilter.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkCubeAxesActor.h"
#include "vtkImagePlaneWidget.h"
#include "vtkImageData.h"
#include "vtkInteractorStyleImage.h"
#include "vtkInteractorStyleTrackballActor.h"
#include "vtkCamera.h"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageMaskSpatialObject.h"
#include "itkCastImageFilter.h"

// Local includes
#include "btkMacro.h"
#include "btkGradientDirection.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceToDiffusionSignalFilter.h"
#include "btkDiffusionSignal.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkSphericalHarmonicsDiffusionDecompositionFilter.h"
#include "btkDiffusionTensorReconstructionFilter.h"
#include "btkDiffusionModel.h"
#include "btkOrientationDiffusionFunctionModel.h"
#include "btkTensorModel.h"
#include "btkImageHelper.h"


typedef double PrecisionType;

typedef itk::Image< short,3 >                                          ImageMask;
typedef itk::ImageMaskSpatialObject< 3 >                               MaskSpatialObject;
typedef itk::CastImageFilter< ImageMask,MaskSpatialObject::ImageType > MaskSpatialObjectCaster;


vtkSmartPointer< vtkActor > surfaceReconstruction(std::vector< btk::GradientDirection > directions, std::vector< float > response, bool normalize=true);
vtkSmartPointer< vtkActor > densityReconstruction(std::vector< btk::GradientDirection > directions, std::vector< PrecisionType > response);
vtkSmartPointer< vtkActor > displayRawSignal(btk::DiffusionSignal::Pointer signal, btk::DiffusionModel::ContinuousIndex cindex);
void displayStats(btk::DiffusionSequence::Pointer dwi, btk::DiffusionModel::ContinuousIndex cindex, btk::DiffusionModel::Pointer model);
std::vector< vtkSmartPointer< vtkActor > > displayDirections(std::vector< btk::GradientDirection > directions);
vtkSmartPointer< vtkImagePlaneWidget > displaySequence(btk::DiffusionSequence::Pointer sequence);
std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelImage(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask);
std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelDirections(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask);


class btkInteractorStyle : public vtkInteractorStyleTrackballCamera
{
    public:
        static btkInteractorStyle *New();
        vtkTypeMacro(btkInteractorStyle, vtkInteractorStyleTrackballCamera);

        virtual void OnKeyPress()
        {
            std::string key = this->Interactor->GetKeySym();

            if(key == "Up")
            {
                unsigned int currentZ = m_ImagePlaneWidget->GetSliceIndex();
                unsigned int newZ = currentZ+1;

                m_ImagePlaneWidget->SetSliceIndex(newZ);

                if(m_ImagePlaneWidget->GetSliceIndex() == newZ)
                {
                    if(m_DisplayModels)
                    {
                        for(unsigned int i = 0; i < m_Models[currentZ].size(); i++)
                        {
                            m_Models[currentZ][i]->VisibilityOff();
                        }

                        for(unsigned int i = 0; i < m_Models[newZ].size(); i++)
                        {
                            m_Models[newZ][i]->VisibilityOn();
                        }
                    }

                    if(m_DisplayMeanDirections)
                    {
                        for(unsigned int i = 0; i < m_MeanDirections[currentZ].size(); i++)
                        {
                            m_MeanDirections[currentZ][i]->VisibilityOff();
                        }

                        for(unsigned int i = 0; i < m_MeanDirections[newZ].size(); i++)
                        {
                            m_MeanDirections[newZ][i]->VisibilityOn();
                        }
                    }
                }

                m_RendererWindowInteractor->Render();
            }

            if(key == "Down")
            {
                unsigned int currentZ = m_ImagePlaneWidget->GetSliceIndex();
                unsigned int newZ = currentZ-1;

                m_ImagePlaneWidget->SetSliceIndex(newZ);

                if(m_ImagePlaneWidget->GetSliceIndex() == newZ)
                {
                    if(m_DisplayModels)
                    {
                        for(unsigned int i = 0; i < m_Models[currentZ].size(); i++)
                        {
                            m_Models[currentZ][i]->VisibilityOff();
                        }

                        for(unsigned int i = 0; i < m_Models[newZ].size(); i++)
                        {
                            m_Models[newZ][i]->VisibilityOn();
                        }
                    }

                    if(m_DisplayMeanDirections)
                    {
                        for(unsigned int i = 0; i < m_MeanDirections[currentZ].size(); i++)
                        {
                            m_MeanDirections[currentZ][i]->VisibilityOff();
                        }

                        for(unsigned int i = 0; i < m_MeanDirections[newZ].size(); i++)
                        {
                            m_MeanDirections[newZ][i]->VisibilityOn();
                        }
                    }
                }

                m_RendererWindowInteractor->Render();
            }

            if(key == "m")
            {
                unsigned int currentZ = m_ImagePlaneWidget->GetSliceIndex();

                if(m_DisplayModels)
                {
                    for(unsigned int i = 0; i < m_Models[currentZ].size(); i++)
                    {
                        m_Models[currentZ][i]->VisibilityOff();
                    }

                    m_DisplayModels = false;
                }
                else // !m_DisplayModels
                {
                    for(unsigned int i = 0; i < m_Models[currentZ].size(); i++)
                    {
                        m_Models[currentZ][i]->VisibilityOn();
                    }

                    m_DisplayModels = true;
                }

                m_RendererWindowInteractor->Render();
            }

            if(key == "d")
            {
                unsigned int currentZ = m_ImagePlaneWidget->GetSliceIndex();

                if(m_DisplayMeanDirections)
                {
                    for(unsigned int i = 0; i < m_MeanDirections[currentZ].size(); i++)
                    {
                        m_MeanDirections[currentZ][i]->VisibilityOff();
                    }

                    m_DisplayMeanDirections = false;
                }
                else // !m_DisplayModels
                {
                    for(unsigned int i = 0; i < m_MeanDirections[currentZ].size(); i++)
                    {
                        m_MeanDirections[currentZ][i]->VisibilityOn();
                    }

                    m_DisplayMeanDirections = true;
                }

                m_RendererWindowInteractor->Render();
            }

            vtkInteractorStyleTrackballCamera::OnKeyPress();
        }

        btkSetMacro(RendererWindowInteractor, vtkSmartPointer< vtkRenderWindowInteractor >);
        btkGetMacro(RendererWindowInteractor, vtkSmartPointer< vtkRenderWindowInteractor >);

        btkSetMacro(ImagePlaneWidget, vtkSmartPointer< vtkImagePlaneWidget >);
        btkGetMacro(ImagePlaneWidget, vtkSmartPointer< vtkImagePlaneWidget >);

        btkSetMacro(Models, std::vector< std::vector< vtkSmartPointer< vtkActor > > >);
        btkGetMacro(Models, std::vector< std::vector< vtkSmartPointer< vtkActor > > >);

        btkSetMacro(MeanDirections, std::vector< std::vector< vtkSmartPointer< vtkActor > > >);
        btkGetMacro(MeanDirections, std::vector< std::vector< vtkSmartPointer< vtkActor > > >);

    protected:
        btkInteractorStyle()
        {
            m_DisplayModels = false;
            m_DisplayMeanDirections = false;
        }

    private:
        vtkSmartPointer< vtkRenderWindowInteractor > m_RendererWindowInteractor;
        vtkSmartPointer< vtkImagePlaneWidget > m_ImagePlaneWidget;
        std::vector< std::vector< vtkSmartPointer< vtkActor > > > m_Models;
        std::vector< std::vector< vtkSmartPointer< vtkActor > > > m_MeanDirections;
        bool m_DisplayModels;
        bool m_DisplayMeanDirections;
};
vtkStandardNewMacro(btkInteractorStyle);


int main(int argc, char *argv[])
{
    //
    // Command line
    //

    // Defines the command line parser
    TCLAP::CmdLine cmd("DWI data viewer", ' ', "2.0");

    // Define mandatory arguments
    TCLAP::ValueArg< std::string > inputSequenceFileNameArg("i", "input_sequence", "Input diffusion sequence filename", true, "", "string", cmd);
    TCLAP::ValueArg< std::string > inputMaskFileNameArg("m", "mask", "Mask filename", false, "", "string", cmd);

    // Options
    TCLAP::ValueArg< double > percentageOfSamplePointsArg("", "percent_of_sample_points", "Percentage of sample points used for modeling display", false, 0.1, "real between 0 and 1", cmd);
    TCLAP::ValueArg< unsigned int > SphericalResolutionArg("", "spherical_resolution", "Spherical resolution used for modeling display", false, 100, "positive integer", cmd);
    TCLAP::SwitchArg meanDirectionsOnlyArg("", "mean_directions_only", "Precompute only mean directions of the modeling and not the all functions (faster)", cmd);

    // Parse command line
    cmd.parse(argc, argv);

    // Get back arguments' values
    std::string inputSequenceFileName = inputSequenceFileNameArg.getValue();
    std::string     inputMaskFileName = inputMaskFileNameArg.getValue();

    double  percentageOfSamplePoints = percentageOfSamplePointsArg.getValue();
    unsigned int SphericalResolution = SphericalResolutionArg.getValue();
    bool          meanDirectionsOnly = meanDirectionsOnlyArg.getValue();


    //
    // Read files
    //

    // Diffusion sequence
    btk::DiffusionSequence::Pointer inputSequence = btk::DiffusionSequenceHelper::ReadSequence(inputSequenceFileName);

    // VTK works in image space coordinates
    inputSequence->ConvertGradientTableToImageCoordinates();

    // Mask (if any)
    ImageMask::Pointer inputMask = NULL;

    if(!inputMaskFileName.empty())
    {
        // Read mask image
        inputMask = btk::ImageHelper< ImageMask >::ReadImage(inputMaskFileName);

        // Set the region of interest
        MaskSpatialObjectCaster::Pointer inputMaskObjectCaster = MaskSpatialObjectCaster::New();
        inputMaskObjectCaster->SetInput(inputMask);
        inputMaskObjectCaster->Update();

        MaskSpatialObject::Pointer inputMaskObject = MaskSpatialObject::New();
        inputMaskObject->SetImage(inputMaskObjectCaster->GetOutput());
        ImageMask::RegionType inputMaskRegion = inputMaskObject->GetAxisAlignedBoundingBoxRegion();

        // Set requested region
        inputMask->SetRequestedRegion(inputMaskRegion);
    }


    //
    // Diffusion modeling
    //

    btkCoutMacro("Estimating modeling ...");

    // Spherical harmonics modeling
    btk::SphericalHarmonicsDiffusionDecompositionFilter::Pointer shModelEstimateFilter = btk::SphericalHarmonicsDiffusionDecompositionFilter::New();
    shModelEstimateFilter->SetInput(inputSequence);
    shModelEstimateFilter->Update();

    // Create model function
    btk::DiffusionModel::Pointer model = NULL;

    // Use fODF model function
    btk::OrientationDiffusionFunctionModel::Pointer fODFModel = btk::OrientationDiffusionFunctionModel::New();
    fODFModel->SetInputModelImage(shModelEstimateFilter->GetOutput());
    fODFModel->SetBValue(inputSequence->GetBValues()[1]);
    fODFModel->UseSharpModelOn();
    fODFModel->SetSphericalResolution(SphericalResolution);
    fODFModel->Update();

    model = fODFModel;

    btkCoutMacro("done.");

    btkCoutMacro("Sampling model volume ...");
    btkTicTocInit();
    std::vector< std::vector< vtkSmartPointer< vtkActor > > > models;
    if(!meanDirectionsOnly)
    {
        btkTic();
        // Build model image
        models = buildModelImage(model, inputSequence, percentageOfSamplePoints, inputMask);
        btkToc();
    }
    btkTic();
    // Build model mean directions
    std::vector< std::vector< vtkSmartPointer< vtkActor > > > meanDirections = buildModelDirections(model, inputSequence, percentageOfSamplePoints, inputMask);
    btkToc();
    btkCoutMacro("done.");


    //
    // Build widgets
    //

    btkCoutMacro("Creating widgets ...");

    // Create image plane widget
    vtkSmartPointer< vtkImagePlaneWidget > imagePlaneWidget = displaySequence(inputSequence);

    btkCoutMacro("done.");


    //
    // Display
    //

    btkCoutMacro("Rendering ...");

    // Renderer
    vtkSmartPointer< vtkRenderer > renderer = vtkSmartPointer< vtkRenderer >::New();
    renderer->SetBackground(0,0,0);

    for(unsigned int z = 0; z < models.size(); z++)
    {
        for(unsigned int i = 0; i < models[z].size(); i++)
        {
            renderer->AddActor(models[z][i]);
            models[z][i]->VisibilityOff();
        }
    }

    for(unsigned int z = 0; z < meanDirections.size(); z++)
    {
        for(unsigned int i = 0; i < meanDirections[z].size(); i++)
        {
            renderer->AddActor(meanDirections[z][i]);
            meanDirections[z][i]->VisibilityOff();
        }
    }

    // Render window
    vtkSmartPointer< vtkRenderWindow > renderWindow = vtkSmartPointer< vtkRenderWindow >::New();
    renderWindow->AddRenderer(renderer);

    // Window interactor
    vtkSmartPointer< vtkRenderWindowInteractor > renderWindowInteractor = vtkSmartPointer< vtkRenderWindowInteractor >::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);

    // Interactor style
    vtkSmartPointer< btkInteractorStyle > interactorStyle = vtkSmartPointer< btkInteractorStyle >::New();
    renderWindowInteractor->SetInteractorStyle(interactorStyle);
    interactorStyle->SetImagePlaneWidget(imagePlaneWidget);
    interactorStyle->SetRendererWindowInteractor(renderWindowInteractor);
    interactorStyle->SetModels(models);
    interactorStyle->SetMeanDirections(meanDirections);

    // Add window interactor to widgets
    imagePlaneWidget->SetInteractor(renderWindowInteractor);
    imagePlaneWidget->On();

    btkCoutMacro("done.");

    // Display
    renderWindow->Render();
    renderWindowInteractor->Initialize();
    renderWindowInteractor->Start();


    return EXIT_SUCCESS;
}

vtkSmartPointer< vtkActor > surfaceReconstruction(std::vector< btk::GradientDirection > directions, std::vector< float > response, bool normalize)
{
    //
    // Parameters
    //

    // Polygonal data and colors
    vtkSmartPointer< vtkPoints >      points = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray > polygons = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkDoubleArray > colors = vtkSmartPointer< vtkDoubleArray >::New();

    // Number of samples per coordinates
    unsigned int nbSamples = directions.size();

    // Compute the regular step on the unit sphere
    unsigned int elevationResolution = static_cast< unsigned int >(std::ceil( std::sqrt(static_cast< PrecisionType >(nbSamples - 2) / 2.0f) ));
    unsigned int   azimuthResolution = 2 * elevationResolution;

    // Structure for storing points during reconstruction
    vtkIdType pid[1];

    // Normalize
    PrecisionType max = 0.f;
    for(unsigned int i = 0; i < response.size(); i++)
    {
        if(max < response[i])
            max = response[i];
    }

    for(unsigned int i = 0; i < directions.size(); i++)
    {
        if(normalize)
        {
            directions[i][0] *= response[i]/max;
            directions[i][1] *= response[i]/max;
            directions[i][2] *= response[i]/max;
        }
        else
        {
            directions[i][0] *= response[i];
            directions[i][1] *= response[i];
            directions[i][2] *= response[i];
        }
    }


    //
    // Surface reconstruction
    //

    // Set first point (north pole)
    points->InsertNextPoint(directions[0][0], directions[0][1], directions[0][2]);
    points->InsertNextPoint(directions[1][0], directions[1][1], directions[1][2]);

    // First circle
    for(unsigned int i = 2; i < 1 + azimuthResolution; i++)
    {
        // Add point
        points->InsertNextPoint(directions[i][0], directions[i][1], directions[i][2]);

        // Define triangle
        vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
        polygon->GetPointIds()->SetNumberOfIds(3);
        polygon->GetPointIds()->SetId(0,0);
        polygon->GetPointIds()->SetId(1,i-1);
        polygon->GetPointIds()->SetId(2,i);

        // Add triangle to polygons
        polygons->InsertNextCell(polygon);
    } // for each point on first circle

    // Define last triangle for first circle
    vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
    polygon->GetPointIds()->SetNumberOfIds(3);
    polygon->GetPointIds()->SetId(0,0);
    polygon->GetPointIds()->SetId(1,azimuthResolution);
    polygon->GetPointIds()->SetId(2,1);

    // Add triangle to polygons
    polygons->InsertNextCell(polygon);

    // Set general points
    for(unsigned int j = 1; j < elevationResolution-1; j++)
    {
        points->InsertNextPoint(directions[1 + j*azimuthResolution][0], directions[1 + j*azimuthResolution][1], directions[1 + j*azimuthResolution][2]);

        for(unsigned int i = 2 + j*azimuthResolution; i < 1 + (j+1) * azimuthResolution; i++)
        {
            points->InsertNextPoint(directions[i][0], directions[i][1], directions[i][2]);

            // Define rectangle
            vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
            polygon->GetPointIds()->SetNumberOfIds(4);
            polygon->GetPointIds()->SetId(0,i-1-azimuthResolution);
            polygon->GetPointIds()->SetId(1,i-1);
            polygon->GetPointIds()->SetId(2,i);
            polygon->GetPointIds()->SetId(3,i-azimuthResolution);

            // Add rectangle to polygons
            polygons->InsertNextCell(polygon);
        } // for each points

        // Define last rectangle for general circle
        vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
        polygon->GetPointIds()->SetNumberOfIds(4);
        polygon->GetPointIds()->SetId(0,j*azimuthResolution);
        polygon->GetPointIds()->SetId(1,(j+1)*azimuthResolution);
        polygon->GetPointIds()->SetId(2,1+j*azimuthResolution);
        polygon->GetPointIds()->SetId(3,1+(j-1)*azimuthResolution);

        // Add rectangle to polygons
        polygons->InsertNextCell(polygon);
    }

    // Last circle
    pid[0] = points->InsertNextPoint(directions[directions.size()-1][0], directions[directions.size()-1][1], directions[directions.size()-1][2]);
    points->InsertNextPoint(directions[1 + (elevationResolution-1)*azimuthResolution][0], directions[1 + (elevationResolution-1)*azimuthResolution][1], directions[1 + (elevationResolution-1)*azimuthResolution][2]);

    for(unsigned int i = 2 + (elevationResolution-2)*azimuthResolution; i < 1 + (elevationResolution-1)*azimuthResolution; i++)
    {
        // Add point
        points->InsertNextPoint(directions[i][0], directions[i][1], directions[i][2]);

        // Define triangle
        vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
        polygon->GetPointIds()->SetNumberOfIds(3);
        polygon->GetPointIds()->SetId(0,pid[0]);
        polygon->GetPointIds()->SetId(1,i-1);
        polygon->GetPointIds()->SetId(2,i);

        // Add triangle to polygons
        polygons->InsertNextCell(polygon);
    } // for each point on first circle

    // Define last triangle for first circle
    vtkSmartPointer< vtkPolygon > polygon2 = vtkSmartPointer< vtkPolygon >::New();
    polygon2->GetPointIds()->SetNumberOfIds(3);
    polygon2->GetPointIds()->SetId(0,pid[0]);
    polygon2->GetPointIds()->SetId(1,(elevationResolution-1)*azimuthResolution);
    polygon2->GetPointIds()->SetId(2,1 + (elevationResolution-2)*azimuthResolution);

    // Add triangle to polygons
    polygons->InsertNextCell(polygon2);


    // Set color map and add it to object
    vtkSmartPointer< vtkColorTransferFunction > map = vtkSmartPointer< vtkColorTransferFunction >::New();
    PrecisionType vmin = *std::min_element(response.begin(), response.end());
    PrecisionType vmax = *std::max_element(response.begin(), response.end());
    map->AddRGBPoint(vmin, 1, 0, 0);
    map->AddRGBPoint(vmax, 0, 0, 1);

    std::vector< float >::iterator itS;

    for(itS = response.begin(); itS != response.end(); itS++)
    {
        double color[3];
        map->GetColor((double)(*itS), color);
        colors->InsertNextTupleValue(color);
    }


    // Add data to vtk structure
    vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer< vtkPolyData >::New();
    polydata->SetPoints(points);
    polydata->SetPolys(polygons);
    polydata->GetPointData()->SetScalars(colors);


    //
    // Prepare VTK for display
    //


    // Mapper
    vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
    mapper->SetInput(polydata);

    // Actor
    vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
    actor->SetMapper(mapper);

    return actor;
}

vtkSmartPointer< vtkActor > densityReconstruction(std::vector< btk::GradientDirection > directions, std::vector< PrecisionType > response)
{
    //
    // Parameters
    //

    // Polygonal data and colors
    vtkSmartPointer< vtkPoints >      points = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray > polygons = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkDoubleArray > colors = vtkSmartPointer< vtkDoubleArray >::New();

    // Number of samples per coordinates
    unsigned int nbSamples = directions.size();

    // Compute the regular step on the unit sphere
    unsigned int elevationResolution = static_cast< unsigned int >(std::ceil( std::sqrt(static_cast< PrecisionType >(nbSamples - 2) / 2.0f) ));
    unsigned int   azimuthResolution = 2 * elevationResolution;

    // Structure for storing points during reconstruction
    vtkIdType pid[1];


    //
    // Surface reconstruction
    //

    // Set first point (north pole)
    points->InsertNextPoint(directions[0][0], directions[0][1], directions[0][2]);
    points->InsertNextPoint(directions[1][0], directions[1][1], directions[1][2]);

    // First circle
    for(unsigned int i = 2; i < 1 + azimuthResolution; i++)
    {
        // Add point
        points->InsertNextPoint(directions[i][0], directions[i][1], directions[i][2]);

        // Define triangle
        vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
        polygon->GetPointIds()->SetNumberOfIds(3);
        polygon->GetPointIds()->SetId(0,0);
        polygon->GetPointIds()->SetId(1,i-1);
        polygon->GetPointIds()->SetId(2,i);

        // Add triangle to polygons
        polygons->InsertNextCell(polygon);
    } // for each point on first circle

    // Define last triangle for first circle
    vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
    polygon->GetPointIds()->SetNumberOfIds(3);
    polygon->GetPointIds()->SetId(0,0);
    polygon->GetPointIds()->SetId(1,azimuthResolution);
    polygon->GetPointIds()->SetId(2,1);

    // Add triangle to polygons
    polygons->InsertNextCell(polygon);

    // Set general points
    for(unsigned int j = 1; j < elevationResolution-1; j++)
    {
        points->InsertNextPoint(directions[1 + j*azimuthResolution][0], directions[1 + j*azimuthResolution][1], directions[1 + j*azimuthResolution][2]);

        for(unsigned int i = 2 + j*azimuthResolution; i < 1 + (j+1) * azimuthResolution; i++)
        {
            points->InsertNextPoint(directions[i][0], directions[i][1], directions[i][2]);

            // Define rectangle
            vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
            polygon->GetPointIds()->SetNumberOfIds(4);
            polygon->GetPointIds()->SetId(0,i-1-azimuthResolution);
            polygon->GetPointIds()->SetId(1,i-1);
            polygon->GetPointIds()->SetId(2,i);
            polygon->GetPointIds()->SetId(3,i-azimuthResolution);

            // Add rectangle to polygons
            polygons->InsertNextCell(polygon);
        } // for each points

        // Define last rectangle for general circle
        vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
        polygon->GetPointIds()->SetNumberOfIds(4);
        polygon->GetPointIds()->SetId(0,j*azimuthResolution);
        polygon->GetPointIds()->SetId(1,(j+1)*azimuthResolution);
        polygon->GetPointIds()->SetId(2,1+j*azimuthResolution);
        polygon->GetPointIds()->SetId(3,1+(j-1)*azimuthResolution);

        // Add rectangle to polygons
        polygons->InsertNextCell(polygon);
    }

    // Last circle
    pid[0] = points->InsertNextPoint(directions[directions.size()-1][0], directions[directions.size()-1][1], directions[directions.size()-1][2]);
    points->InsertNextPoint(directions[1 + (elevationResolution-1)*azimuthResolution][0], directions[1 + (elevationResolution-1)*azimuthResolution][1], directions[1 + (elevationResolution-1)*azimuthResolution][2]);

    for(unsigned int i = 2 + (elevationResolution-2)*azimuthResolution; i < 1 + (elevationResolution-1)*azimuthResolution; i++)
    {
        // Add point
        points->InsertNextPoint(directions[i][0], directions[i][1], directions[i][2]);

        // Define triangle
        vtkSmartPointer< vtkPolygon > polygon = vtkSmartPointer< vtkPolygon >::New();
        polygon->GetPointIds()->SetNumberOfIds(3);
        polygon->GetPointIds()->SetId(0,pid[0]);
        polygon->GetPointIds()->SetId(1,i-1);
        polygon->GetPointIds()->SetId(2,i);

        // Add triangle to polygons
        polygons->InsertNextCell(polygon);
    } // for each point on first circle

    // Define last triangle for first circle
    vtkSmartPointer< vtkPolygon > polygon2 = vtkSmartPointer< vtkPolygon >::New();
    polygon2->GetPointIds()->SetNumberOfIds(3);
    polygon2->GetPointIds()->SetId(0,pid[0]);
    polygon2->GetPointIds()->SetId(1,(elevationResolution-1)*azimuthResolution);
    polygon2->GetPointIds()->SetId(2,1 + (elevationResolution-2)*azimuthResolution);

    // Add triangle to polygons
    polygons->InsertNextCell(polygon2);


    // Set color map and add it to object
    vtkSmartPointer< vtkColorTransferFunction > map = vtkSmartPointer< vtkColorTransferFunction >::New();
    PrecisionType vmin = *std::min_element(response.begin(), response.end());
    PrecisionType vmax = *std::max_element(response.begin(), response.end());
    map->AddRGBPoint(vmin, 1, 0, 0);
    map->AddRGBPoint(vmax, 0, 0, 1);
    btkCoutVariable(vmin);btkCoutVariable(vmax);
    std::vector< PrecisionType >::iterator itS;

    for(itS = response.begin(); itS != response.end(); itS++)
    {
        double color[3];
        map->GetColor((double)(*itS), color);
        colors->InsertNextTupleValue(color);
    }


    // Add data to vtk structure
    vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer< vtkPolyData >::New();
    polydata->SetPoints(points);
    polydata->SetPolys(polygons);
    polydata->GetPointData()->SetScalars(colors);


    //
    // Prepare VTK for display
    //


    // Mapper
    vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
    mapper->SetInput(polydata);

    // Actor
    vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
    actor->SetMapper(mapper);

    return actor;
}

vtkSmartPointer< vtkActor > displayRawSignal(btk::DiffusionSignal::Pointer signal, btk::DiffusionModel::ContinuousIndex cindex)
{
    // Polydata structures
    vtkSmartPointer< vtkPolyData > pointsData = vtkSmartPointer< vtkPolyData >::New();
    vtkSmartPointer< vtkCellArray >    vertex = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkPoints >       points = vtkSmartPointer< vtkPoints >::New();

    // Get gradient table
    std::vector< btk::GradientDirection > gradientTable = signal->GetGradientTable();

    // Get current index
    btk::DiffusionSignal::IndexType currentIndex;
    currentIndex[0] = cindex[0]; currentIndex[1] = cindex[1]; currentIndex[2] = cindex[2];

    btk::DiffusionSignal::PixelType currentResponsePixel = signal->GetPixel(currentIndex);

    // Build spherical samples from signal response
    for(unsigned int i = 0; i < gradientTable.size(); i++)
    {
        btk::GradientDirection currentDirection = gradientTable[i];
        PrecisionType currentResponse = currentResponsePixel[i];

        double coord[3] = { 0,0,0 };
        coord[0] = currentDirection[0] * currentResponse;
        coord[1] = currentDirection[1] * currentResponse;
        coord[2] = currentDirection[2] * currentResponse;

        vtkIdType id = points->InsertNextPoint(coord[0], coord[1], coord[2]);
        vertex->InsertNextCell(1,&id);
    }

    pointsData->SetPoints(points);
    pointsData->SetVerts(vertex);

    // Create mapper
    vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
    mapper->SetInput(pointsData);

    // Create actor
    vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetPointSize(3);
    actor->GetProperty()->SetColor(1,1,1);


    return actor;
}

void displayStats(btk::DiffusionSequence::Pointer dwi, btk::DiffusionModel::ContinuousIndex cindex, btk::DiffusionModel::Pointer model)
{
    // Get gradient table
    std::vector< btk::GradientDirection > gradientTable = dwi->GetGradientTable();

    // Get current index
    btk::DiffusionSequence::IndexType currentIndex;
    currentIndex[0] = cindex[0]; currentIndex[1] = cindex[1]; currentIndex[2] = cindex[2];

    btk::DiffusionSequence::IndexType baselineIndex = currentIndex;
    baselineIndex[3] = 0;

    // Compute stats
    double meanError = 0;

    for(unsigned int i = 1; i < gradientTable.size(); i++)
    {
        btk::GradientDirection currentDirection = gradientTable[i];

        currentIndex[3] = i;
        double currentResponse = static_cast< double >(dwi->GetPixel(currentIndex)) / static_cast< double >(dwi->GetPixel(baselineIndex));
        double   modelResponse = model->SignalAt(cindex, currentDirection);

        // Compute stats
        double error = std::abs(currentResponse - modelResponse);

        meanError += error;

        // Display in command line
//        btkCoutVariable(i);
//        btkCoutVariable(error);
    }

    meanError /= gradientTable.size();
    btkCoutVariable(meanError);
}

std::vector< vtkSmartPointer< vtkActor > > displayDirections(std::vector< btk::GradientDirection > directions)
{
    double p0[3] = { 0, 0, 0 };
    double p1[3];
    std::vector< vtkSmartPointer< vtkActor > > actors;

//    btkCoutVariable(directions.size());

    for(unsigned int i = 0; i < directions.size(); i++)
    {
//        btkCoutVariable(directions[i]);

        p1[0] = 1.2*directions[i][0]; p1[1] = 1.2*directions[i][1]; p1[2] = 1.2*directions[i][2];

        vtkSmartPointer< vtkLineSource > line = vtkSmartPointer< vtkLineSource >::New();
        line->SetPoint1(p0);
        line->SetPoint2(p1);
        line->Update();

        vtkSmartPointer< vtkTubeFilter > tube = vtkSmartPointer< vtkTubeFilter >::New();
        tube->SetInputConnection(line->GetOutputPort());
        tube->SetRadius(0.05);
        tube->SetNumberOfSides(10);

        vtkSmartPointer< vtkPolyDataMapper > mapper = vtkSmartPointer< vtkPolyDataMapper >::New();
        mapper->SetInputConnection(tube->GetOutputPort());

        vtkSmartPointer< vtkActor > actor = vtkSmartPointer< vtkActor >::New();
        actor->SetMapper(mapper);

        actors.push_back(actor);
    }

    return actors;
}

vtkSmartPointer< vtkImagePlaneWidget > displaySequence(btk::DiffusionSequence::Pointer sequence)
{
    // Get image information
    btk::DiffusionSequence::RegionType   region = sequence->GetLargestPossibleRegion();
    btk::DiffusionSequence::SizeType       size = region.GetSize();
    btk::DiffusionSequence::SpacingType spacing = sequence->GetSpacing();

    // Create vtk image of baseline image
    vtkSmartPointer< vtkImageData > image = vtkSmartPointer< vtkImageData >::New();
    image->SetDimensions(size[0], size[1], size[2]);
    image->SetSpacing(spacing[0], spacing[1], spacing[2]);
    image->SetScalarTypeToUnsignedChar();
    image->SetNumberOfScalarComponents(1);
    image->AllocateScalars();

    // Fill vtk image with intensity scaled baseline image
    region.SetSize(3,1);
    itk::ImageRegionConstIterator< btk::DiffusionSequence > it(sequence, region);

    btk::DiffusionSequence::PixelType min = itk::NumericTraits< btk::DiffusionSequence::PixelType >::max(), max = itk::NumericTraits< btk::DiffusionSequence::PixelType >::min();

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        if(max < it.Get())
        {
            max = it.Get();
        }

        if(min > it.Get())
        {
            min = it.Get();
        }
    }

    unsigned char *ptr = static_cast< unsigned char * >(image->GetScalarPointer());

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        double value = 255.0 * (it.Get() - min) / (max - min);
        *ptr++ = static_cast< unsigned char >( value );
    }


    // Create image plane widget
    vtkSmartPointer< vtkImagePlaneWidget > widget = vtkSmartPointer< vtkImagePlaneWidget >::New();
    widget->TextureInterpolateOff();
    widget->SetInput(image);
    widget->RestrictPlaneToVolumeOn();
    widget->SetResliceInterpolateToNearestNeighbour();
    widget->DisplayTextOn();
    widget->SetPlaneOrientationToZAxes();
    widget->PlaceWidget();

    return widget;
}

std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelImage(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask)
{
    // Create vector of actors
    std::vector< std::vector< vtkSmartPointer< vtkActor > > > *actors = new std::vector< std::vector< vtkSmartPointer< vtkActor > > >;

    // Set size and spacing
    btk::DiffusionSequence::SizeType size = sequence->GetLargestPossibleRegion().GetSize();
    size[0] *= factor; size[1] *= factor;

    btk::DiffusionSequence::SpacingType spacing = sequence->GetSpacing();
    spacing[0] /= factor; spacing[1] /= factor;

    // Build all surfaces
    for(unsigned int z = 0; z < size[2]; z++)
    {
        actors->push_back(std::vector< vtkSmartPointer< vtkActor > >());

        for(unsigned int x = 0; x < size[0]; x++)
        {
            for(unsigned int y = 0; y < size[1]; y++)
            {
                btk::DiffusionModel::ContinuousIndex cindex;
                cindex[0] = x/factor; cindex[1] = y/factor; cindex[2] = z;

                ImageMask::IndexType index;
                index[0] = cindex[0]; index[1] = cindex[1]; index[2] = cindex[2];

                if(mask.IsNull() || mask->GetPixel(index) > 0)
                {
                    vtkSmartPointer< vtkActor > actor = surfaceReconstruction(model->GetDirections(), model->ModelAt(cindex));
                    actor->SetPosition(x*spacing[0], y*spacing[1], z*spacing[2]);
                    (*actors)[z].push_back(actor);
                }
            }
        }
    }

    return *actors;
}

std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelDirections(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask)
{
    // Create vector of actors
    std::vector< std::vector< vtkSmartPointer< vtkActor > > > *actors = new std::vector< std::vector< vtkSmartPointer< vtkActor > > >;

    // Set size and spacing
    btk::DiffusionSequence::SizeType size = sequence->GetLargestPossibleRegion().GetSize();
    size[0] *= factor; size[1] *= factor;

    btk::DiffusionSequence::SpacingType spacing = sequence->GetSpacing();
    spacing[0] /= factor; spacing[1] /= factor;

    // Build all surfaces
    for(unsigned int z = 0; z < size[2]; z++)
    {
        actors->push_back(std::vector< vtkSmartPointer< vtkActor > >());

        for(unsigned int x = 0; x < size[0]; x++)
        {
            for(unsigned int y = 0; y < size[1]; y++)
            {
                btk::DiffusionModel::ContinuousIndex cindex;
                cindex[0] = x/factor; cindex[1] = y/factor; cindex[2] = z;

                ImageMask::IndexType index;
                index[0] = cindex[0]; index[1] = cindex[1]; index[2] = cindex[2];

                if(mask.IsNull() || mask->GetPixel(index) > 0)
                {
                    std::vector< btk::GradientDirection >   meanDirection = model->MeanDirectionsAt(cindex);
                    std::vector< vtkSmartPointer< vtkActor > > meanActors = displayDirections(meanDirection);

                    for(unsigned int i = 0; i < meanActors.size(); i++)
                    {
                        meanActors[i]->SetPosition(x*spacing[0], y*spacing[1], z*spacing[2]);
                        meanActors[i]->GetProperty()->SetColor(std::abs(meanDirection[i][0]), std::abs(meanDirection[i][1]), std::abs(meanDirection[i][2]));
                        (*actors)[z].push_back(meanActors[i]);
                    }
                }
            }
        }
    }

    return *actors;
}
