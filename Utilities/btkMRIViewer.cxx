/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/04/2013
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
#include "vtkArrowSource.h"
#include "vtkGlyph2D.h"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageMaskSpatialObject.h"
#include "itkCastImageFilter.h"
#include "itkDisplacementFieldTransform.h"

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

typedef itk::Image< double,3 >                                         Image;
typedef itk::Image< short,3 >                                          ImageMask;
typedef itk::ImageMaskSpatialObject< 3 >                               MaskSpatialObject;
typedef itk::CastImageFilter< ImageMask,MaskSpatialObject::ImageType > MaskSpatialObjectCaster;

typedef itk::DisplacementFieldTransform< double,3 >::DisplacementFieldType DisplacementField;


vtkSmartPointer< vtkImagePlaneWidget > displayImage(Image::Pointer image);
//std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelImage(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask);
//std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelDirections(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask);


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
                    if(m_DisplayDisplacementField)
                    {
                        for(unsigned int i = 0; i < m_DisplacementField[currentZ].size(); i++)
                        {
                            m_DisplacementField[currentZ][i]->VisibilityOff();
                        }

                        for(unsigned int i = 0; i < m_DisplacementField[newZ].size(); i++)
                        {
                            m_DisplacementField[newZ][i]->VisibilityOn();
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
                    if(m_DisplayDisplacementField)
                    {
                        for(unsigned int i = 0; i < m_DisplacementField[currentZ].size(); i++)
                        {
                            m_DisplacementField[currentZ][i]->VisibilityOff();
                        }

                        for(unsigned int i = 0; i < m_DisplacementField[newZ].size(); i++)
                        {
                            m_DisplacementField[newZ][i]->VisibilityOn();
                        }
                    }
                }

                m_RendererWindowInteractor->Render();
            }

            if(key == "f")
            {
                unsigned int currentZ = m_ImagePlaneWidget->GetSliceIndex();

                if(m_DisplayDisplacementField)
                {
                    for(unsigned int i = 0; i < m_DisplacementField[currentZ].size(); i++)
                    {
                        m_DisplacementField[currentZ][i]->VisibilityOff();
                    }

                    m_DisplayDisplacementField = false;
                }
                else // !m_DisplayDisplacementField
                {
                    for(unsigned int i = 0; i < m_DisplacementField[currentZ].size(); i++)
                    {
                        m_DisplacementField[currentZ][i]->VisibilityOn();
                    }

                    m_DisplayDisplacementField = true;
                }

                m_RendererWindowInteractor->Render();
            }

            vtkInteractorStyleTrackballCamera::OnKeyPress();
        }

        btkSetMacro(RendererWindowInteractor, vtkSmartPointer< vtkRenderWindowInteractor >);
        btkGetMacro(RendererWindowInteractor, vtkSmartPointer< vtkRenderWindowInteractor >);

        btkSetMacro(ImagePlaneWidget, vtkSmartPointer< vtkImagePlaneWidget >);
        btkGetMacro(ImagePlaneWidget, vtkSmartPointer< vtkImagePlaneWidget >);

        btkSetMacro(DisplacementField, std::vector< std::vector< vtkSmartPointer< vtkActor > > >);
        btkGetMacro(DisplacementField, std::vector< std::vector< vtkSmartPointer< vtkActor > > >);

    protected:
        btkInteractorStyle()
        {
            m_DisplayDisplacementField = false;
        }

    private:
        vtkSmartPointer< vtkRenderWindowInteractor > m_RendererWindowInteractor;
        vtkSmartPointer< vtkImagePlaneWidget > m_ImagePlaneWidget;
        std::vector< std::vector< vtkSmartPointer< vtkActor > > > m_DisplacementField;
        bool m_DisplayDisplacementField;
};
vtkStandardNewMacro(btkInteractorStyle);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#include <vtkVersion.h>
//#include <vtkArrowSource.h>
//#include <vtkCellArray.h>
//#include <vtkGlyph2D.h>
//#include <vtkPointData.h>
//#include <vtkImageData.h>
//#include <vtkInteractorStyleImage.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkPolyData.h>
//#include <vtkPoints.h>
//#include <vtkRenderer.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkSmartPointer.h>
//#include <vtkXMLPolyDataWriter.h>

//int main(int, char *[])
//{
//  // Create an image
//  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

//  // Specify the size of the image data
//  image->SetDimensions(50,50,1);

//#if VTK_MAJOR_VERSION <= 5
//  image->SetNumberOfScalarComponents(3);
//  image->SetScalarTypeToFloat();
//  image->AllocateScalars();
//#else
//  image->AllocateScalars(VTK_FLOAT,3);
//#endif
//  int* dims = image->GetDimensions();

//  // Zero the image
//  for (int y = 0; y < dims[1]; y++)
//    {
//    for (int x = 0; x < dims[0]; x++)
//      {
//      float* pixel = static_cast<float*>(image->GetScalarPointer(x,y,0));
//      pixel[0] = 0.0;
//      pixel[1] = 0.0;
//      pixel[2] = 0.0;
//      }
//    }

//  {
//  float* pixel = static_cast<float*>(image->GetScalarPointer(20,20,0));
//  pixel[0] = 0.1;
//  pixel[1] = 0.1;
//  pixel[2] = 100.0;
//  }

//  {
//  float* pixel = static_cast<float*>(image->GetScalarPointer(30,30,0));
//  pixel[0] = 10.0;
//  pixel[1] = 10.0;
//  pixel[2] = -20.0;
//  }

//  image->GetPointData()->SetActiveVectors("ImageScalars");

//  // Setup the arrows
//  vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
//  arrowSource->Update();

//  vtkSmartPointer<vtkGlyph3D> glyphFilter = vtkSmartPointer<vtkGlyph3D>::New();
//  glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
//  glyphFilter->OrientOn();
//  glyphFilter->SetVectorModeToUseVector();
//  glyphFilter->SetScaleModeToDataScalingOff();
//#if VTK_MAJOR_VERSION <= 5
//  glyphFilter->SetInputConnection(image->GetProducerPort());
//#else
//  glyphFilter->SetInputData(image);
//#endif
//  glyphFilter->Update();

//  // Create actors

//  vtkSmartPointer<vtkPolyDataMapper> vectorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
//  vectorMapper->SetInputConnection(glyphFilter->GetOutputPort());
//  vtkSmartPointer<vtkActor> vectorActor = vtkSmartPointer<vtkActor>::New();
//  vectorActor->SetMapper(vectorMapper);

//  // Setup renderer
//  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//  renderer->AddViewProp(vectorActor);
//  renderer->ResetCamera();

//  // Setup render window
//  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//  renderWindow->AddRenderer(renderer);

//  // Setup render window interactor
//  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//  vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
//  //renderWindowInteractor->SetInteractorStyle(style);

//  // Render and start interaction
//  renderWindowInteractor->SetRenderWindow(renderWindow);
//  renderWindowInteractor->Initialize();

//  renderWindowInteractor->Start();

//  return EXIT_SUCCESS;
//}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    //
    // Command line
    //

    // Defines the command line parser
    TCLAP::CmdLine cmd("DWI data viewer", ' ', "2.0");

    // Define mandatory arguments
    TCLAP::ValueArg< std::string > inputImageFileNameArg("i", "input_sequence", "Input diffusion sequence filename", true, "", "string", cmd);
    TCLAP::ValueArg< std::string > inputMaskFileNameArg("m", "mask", "Mask filename", false, "", "string", cmd);

    // Options
    TCLAP::ValueArg< double > percentageOfSamplePointsArg("", "percent_of_sample_points", "Percentage of sample points used for modeling display", false, 0.1, "real between 0 and 1", cmd);
    TCLAP::ValueArg< std::string >   displacementFieldArg("f", "displacement_field", "Displacement field to display", false, "", "string", cmd);

    // Parse command line
    cmd.parse(argc, argv);

    // Get back arguments' values
    std::string inputImageFileName = inputImageFileNameArg.getValue();
    std::string     inputMaskFileName = inputMaskFileNameArg.getValue();

    double  percentageOfSamplePoints = percentageOfSamplePointsArg.getValue();
    std::string displacementField = displacementFieldArg.getValue();


    //
    // Read files
    //

    // Diffusion sequence
    Image::Pointer inputImage = btk::ImageHelper< Image >::ReadImage(inputImageFileName);

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

    // Displacement field (if any)
    DisplacementField::Pointer field = NULL;

    if(!displacementField.empty())
    {
        field = btk::ImageHelper< DisplacementField >::ReadImage(displacementField);
    }


    //
    // Field reconstruction (glyph)
    //

    vtkSmartPointer< vtkActor > vectorActor = NULL;

    if(!displacementField.empty())
    {
        // Get image information
        DisplacementField::RegionType       region = field->GetLargestPossibleRegion();
        DisplacementField::SizeType           size = region.GetSize();
        DisplacementField::SpacingType     spacing = field->GetSpacing();
        DisplacementField::DirectionType direction = field->GetDirection();

        // Create vtk image of baseline image
        vtkSmartPointer< vtkImageData > vtkfield = vtkSmartPointer< vtkImageData >::New();
        vtkfield->SetDimensions(size[0], size[1], size[2]);
        vtkfield->SetSpacing(spacing[0], spacing[1], spacing[2]);
        vtkfield->SetScalarTypeToFloat();
        vtkfield->SetNumberOfScalarComponents(3);
        vtkfield->AllocateScalars();

        unsigned int voxelIndex = 0;
        unsigned int frequency = static_cast< unsigned int >(1.0/percentageOfSamplePoints);btkCoutVariable(frequency);

        // Fill vtk field
        for(unsigned int z = 0; z < size[2]; z++)
        {
            for(unsigned int y = 0; y < size[1]; y++)
            {
                for(unsigned int x = 0; x < size[0]; x++)
                {
                    ImageMask::IndexType index;
                    index[0] = x; index[1] = y; index[2] = z;

                    if(inputMask->GetPixel(index) > 0)
                    {
                        float *outputPixel = static_cast< float * >(vtkfield->GetScalarPointer(x,y,z));

                        if(voxelIndex % frequency == 0)
                        {
                            DisplacementField::IndexType index;
                            index[0] = x; index[1] = y; index[2] = z;
                            DisplacementField::PixelType inputPixel = field->GetPixel(index);

                            outputPixel[0] = direction(0,0) * inputPixel[0];
                            outputPixel[1] = direction(1,1) * inputPixel[1];
                            outputPixel[2] = direction(2,2) * inputPixel[2];
                        }
                        else
                        {
                            outputPixel[0] = 0.0;
                            outputPixel[1] = 0.0;
                            outputPixel[2] = 0.0;
                        }

                        voxelIndex++;
                    }
                }
            }
        }

        vtkfield->GetPointData()->SetActiveVectors("ImageScalars");

        // Setup the arrows
        vtkSmartPointer< vtkArrowSource > arrowSource = vtkSmartPointer< vtkArrowSource >::New();
        arrowSource->Update();

        vtkSmartPointer< vtkGlyph3D > glyphFilter = vtkSmartPointer< vtkGlyph3D >::New();
        glyphFilter->SetSourceConnection(arrowSource->GetOutputPort());
        glyphFilter->SetInputConnection(vtkfield->GetProducerPort());
        glyphFilter->OrientOn();
        glyphFilter->SetVectorModeToUseVector();
        glyphFilter->SetScaleModeToScaleByVector();
        glyphFilter->SetColorModeToColorByScale();
        glyphFilter->ClampingOn();
        glyphFilter->Update();

        // Setup actor
        vtkSmartPointer< vtkPolyDataMapper > vectorMapper = vtkSmartPointer< vtkPolyDataMapper >::New();
        vectorMapper->SetInputConnection(glyphFilter->GetOutputPort());

        vectorActor = vtkSmartPointer< vtkActor >::New();
        vectorActor->SetMapper(vectorMapper);
    }



    //
    // Build widgets
    //

    btkCoutMacro("Creating widgets ...");

    // Create image plane widget
    vtkSmartPointer< vtkImagePlaneWidget > imagePlaneWidget = displayImage(inputImage);

    btkCoutMacro("done.");


    //
    // Display
    //

    btkCoutMacro("Rendering ...");

    // Renderer
    vtkSmartPointer< vtkRenderer > renderer = vtkSmartPointer< vtkRenderer >::New();
    renderer->SetBackground(0,0,0);

    if(!displacementField.empty())
    {
        renderer->AddViewProp(vectorActor);
    }

//    for(unsigned int z = 0; z < models.size(); z++)
//    {
//        for(unsigned int i = 0; i < models[z].size(); i++)
//        {
//            renderer->AddActor(models[z][i]);
//            models[z][i]->VisibilityOff();
//        }
//    }

//    for(unsigned int z = 0; z < meanDirections.size(); z++)
//    {
//        for(unsigned int i = 0; i < meanDirections[z].size(); i++)
//        {
//            renderer->AddActor(meanDirections[z][i]);
//            meanDirections[z][i]->VisibilityOff();
//        }
//    }

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
//    interactorStyle->SetModels(models);
//    interactorStyle->SetMeanDirections(meanDirections);

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


vtkSmartPointer< vtkImagePlaneWidget > displayImage(Image::Pointer image)
{
    // Get image information
    Image::RegionType   region = image->GetLargestPossibleRegion();
    Image::SizeType       size = region.GetSize();
    Image::SpacingType spacing = image->GetSpacing();

    // Create vtk image of baseline image
    vtkSmartPointer< vtkImageData > vtkimage = vtkSmartPointer< vtkImageData >::New();
    vtkimage->SetDimensions(size[0], size[1], size[2]);
    vtkimage->SetSpacing(spacing[0], spacing[1], spacing[2]);
    vtkimage->SetScalarTypeToUnsignedChar();
    vtkimage->SetNumberOfScalarComponents(1);
    vtkimage->AllocateScalars();

    // Fill vtk image with intensity scaled baseline image
    itk::ImageRegionConstIterator< Image > it(image, region);

    Image::PixelType min = itk::NumericTraits< Image::PixelType >::max(), max = itk::NumericTraits< Image::PixelType >::min();

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

    unsigned char *ptr = static_cast< unsigned char * >(vtkimage->GetScalarPointer());

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        double value = 255.0 * (it.Get() - min) / (max - min);
        *ptr++ = static_cast< unsigned char >( value );
    }


    // Create image plane widget
    vtkSmartPointer< vtkImagePlaneWidget > widget = vtkSmartPointer< vtkImagePlaneWidget >::New();
    widget->TextureInterpolateOff();
    widget->SetInput(vtkimage);
    widget->RestrictPlaneToVolumeOn();
    widget->SetResliceInterpolateToNearestNeighbour();
    widget->DisplayTextOn();
    widget->SetPlaneOrientationToZAxes();
    widget->PlaceWidget();

    return widget;
}

//std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelImage(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask)
//{
//    // Create vector of actors
//    std::vector< std::vector< vtkSmartPointer< vtkActor > > > *actors = new std::vector< std::vector< vtkSmartPointer< vtkActor > > >;

//    // Set size and spacing
//    btk::DiffusionSequence::SizeType size = sequence->GetLargestPossibleRegion().GetSize();
//    size[0] *= factor; size[1] *= factor;

//    btk::DiffusionSequence::SpacingType spacing = sequence->GetSpacing();
//    spacing[0] /= factor; spacing[1] /= factor;

//    // Build all surfaces
//    for(unsigned int z = 0; z < size[2]; z++)
//    {
//        actors->push_back(std::vector< vtkSmartPointer< vtkActor > >());

//        for(unsigned int x = 0; x < size[0]; x++)
//        {
//            for(unsigned int y = 0; y < size[1]; y++)
//            {
//                btk::DiffusionModel::ContinuousIndex cindex;
//                cindex[0] = x/factor; cindex[1] = y/factor; cindex[2] = z;

//                ImageMask::IndexType index;
//                index[0] = cindex[0]; index[1] = cindex[1]; index[2] = cindex[2];

//                if(mask.IsNull() || mask->GetPixel(index) > 0)
//                {
//                    vtkSmartPointer< vtkActor > actor = surfaceReconstruction(model->GetDirections(), model->ModelAt(cindex));
//                    actor->SetPosition(x*spacing[0], y*spacing[1], z*spacing[2]);
//                    (*actors)[z].push_back(actor);
//                }
//            }
//        }
//    }

//    return *actors;
//}

//std::vector< std::vector< vtkSmartPointer< vtkActor > > > &buildModelDirections(btk::DiffusionModel::Pointer model, btk::DiffusionSequence::Pointer sequence, double factor, ImageMask::Pointer mask)
//{
//    // Create vector of actors
//    std::vector< std::vector< vtkSmartPointer< vtkActor > > > *actors = new std::vector< std::vector< vtkSmartPointer< vtkActor > > >;

//    // Set size and spacing
//    btk::DiffusionSequence::SizeType size = sequence->GetLargestPossibleRegion().GetSize();
//    size[0] *= factor; size[1] *= factor;

//    btk::DiffusionSequence::SpacingType spacing = sequence->GetSpacing();
//    spacing[0] /= factor; spacing[1] /= factor;

//    // Build all surfaces
//    for(unsigned int z = 0; z < size[2]; z++)
//    {
//        actors->push_back(std::vector< vtkSmartPointer< vtkActor > >());

//        for(unsigned int x = 0; x < size[0]; x++)
//        {
//            for(unsigned int y = 0; y < size[1]; y++)
//            {
//                btk::DiffusionModel::ContinuousIndex cindex;
//                cindex[0] = x/factor; cindex[1] = y/factor; cindex[2] = z;

//                ImageMask::IndexType index;
//                index[0] = cindex[0]; index[1] = cindex[1]; index[2] = cindex[2];

//                if(mask.IsNull() || mask->GetPixel(index) > 0)
//                {
//                    std::vector< btk::GradientDirection >   meanDirection = model->MeanDirectionsAt(cindex);
//                    std::vector< vtkSmartPointer< vtkActor > > meanActors = displayDirections(meanDirection);

//                    for(unsigned int i = 0; i < meanActors.size(); i++)
//                    {
//                        meanActors[i]->SetPosition(x*spacing[0], y*spacing[1], z*spacing[2]);
//                        meanActors[i]->GetProperty()->SetColor(std::abs(meanDirection[i][0]), std::abs(meanDirection[i][1]), std::abs(meanDirection[i][2]));
//                        (*actors)[z].push_back(meanActors[i]);
//                    }
//                }
//            }
//        }
//    }

//    return *actors;
//}
