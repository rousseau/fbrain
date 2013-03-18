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

#include "btkTractographyAlgorithm.h"


// STL includes
#include "algorithm"

// ITK includes
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageMaskSpatialObject.h"
#include "itkCastImageFilter.h"

// VTK includes
#include "vtkPolyData.h"

// Types used by this filter
typedef itk::ResampleImageFilter< btk::TractographyAlgorithm::LabelImage,btk::TractographyAlgorithm::LabelImage > LabelResampler;
typedef itk::NearestNeighborInterpolateImageFunction< btk::TractographyAlgorithm::LabelImage,double > LabelInterpolator;
typedef itk::ImageRegionIteratorWithIndex< btk::TractographyAlgorithm::LabelImage > LabelIterator;
typedef itk::ImageMaskSpatialObject< btk::TractographyAlgorithm::LabelImage::ImageDimension > LabelSpatialObject;
typedef itk::CastImageFilter< btk::TractographyAlgorithm::LabelImage,LabelSpatialObject::ImageType > LabelSpatialObjectCaster;

// Mutex used by this filter in multi-threaded processing
static itk::SimpleFastMutexLock mutex;


namespace btk
{

TractographyAlgorithm::TractographyAlgorithm() : m_RegionsOfInterest(NULL), m_SeedSpacing(1), m_SeedLabels(), m_OutputFibers(), m_OutputIndicesOfLabels(), m_DiffusionSignal(NULL), m_DiffusionModel(NULL), m_Mask(NULL), Superclass()
{
    // ----
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::Update()
{
    // Initialize filter
    this->Initialize();

    // Resample label image to seed spacing
    this->ResampleLabelImage();

    // Select the region on which will be applied the algorithm
    LabelSpatialObjectCaster::Pointer objectCaster = LabelSpatialObjectCaster::New();
    objectCaster->SetInput(m_RegionsOfInterest);
    objectCaster->Update();
    LabelSpatialObject::Pointer object = LabelSpatialObject::New();
    object->SetImage(objectCaster->GetOutput());
    m_RegionsOfInterest->SetRequestedRegion(object->GetAxisAlignedBoundingBoxRegion());

    // Initialize the progress bar
    this->SetProgress(0);
    m_ProgressStep = 1.0 / static_cast< double >(m_RegionsOfInterest->GetRequestedRegion().GetNumberOfPixels());

    // Set up the multithreaded processing
    ThreadStruct str;
    str.Filter = this;

    this->GetMultiThreader()->SetNumberOfThreads(this->GetNumberOfThreads());
    this->GetMultiThreader()->SetSingleMethod(this->ThreaderCallback, &str);

    // Multithread the execution
    this->GetMultiThreader()->SingleMethodExecute();

    // Update progress to 1
    this->SetProgress(1);
    this->InvokeEvent(itk::ProgressEvent());
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::Initialize()
{
    // ----
}

//----------------------------------------------------------------------------------------

ITK_THREAD_RETURN_TYPE TractographyAlgorithm::ThreaderCallback(void *arg)
{
    ThreadStruct *str;
    itk::ThreadIdType  total, threadId, threadCount;

    threadId    = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->ThreadID;
    threadCount = ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->NumberOfThreads;

    str = (ThreadStruct *)( ( (itk::MultiThreader::ThreadInfoStruct *)( arg ) )->UserData );

    // execute the actual method with appropriate roi region
    // first find out how many pieces extent can be split into.
    LabelImage::RegionType splitRegion;
    total = str->Filter->SplitRequestedRegion(threadId, threadCount, splitRegion);

    if (threadId < total)
    {
        str->Filter->ThreadedGenerateData(splitRegion, threadId);
    }
    // else
    //   {
    //   otherwise don't use this thread. Sometimes the threads dont
    //   break up very well and it is just as efficient to leave a
    //   few threads idle.
    //   }

    return ITK_THREAD_RETURN_VALUE;
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::ThreadedGenerateData(const LabelImage::RegionType &region, itk::ThreadIdType threadId)
{
    LabelIterator it(m_RegionsOfInterest, region);

    unsigned int numberOfSelectedLabels = m_SeedLabels.size();

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        unsigned int label = it.Get();

        // If the label of the current voxel has to bee processed
        if(label > 0 && (numberOfSelectedLabels == 0 || *std::find(m_SeedLabels.begin(), m_SeedLabels.end(), label) == label))
        {
            // Get the physical point
            Self::LabelImage::IndexType labelIndex = it.GetIndex();

            Self::PhysicalPoint point;
            m_RegionsOfInterest->TransformIndexToPhysicalPoint(labelIndex, point);

            // Check if the physical point is in the mask
            Self::MaskImage::IndexType maskIndex;
            m_Mask->TransformPhysicalPointToIndex(point, maskIndex);

            if(m_Mask->GetPixel(maskIndex)  != 0)
            {
                // Initialize vector
                vtkSmartPointer< vtkPolyData > currentFiber = NULL;

                // Start tractography from seed point
                currentFiber = PropagateSeed(point);

                // Save data
                mutex.Lock();

                for(int l = m_OutputFibers.size(); l < label; l++)
                {
                    m_OutputIndicesOfLabels.push_back(-1);
                }

                if(m_OutputIndicesOfLabels[label-1] == -1)
                {
                    m_OutputFibers.push_back(vtkSmartPointer< vtkAppendPolyData >::New());
                    m_OutputIndicesOfLabels[label-1] = m_OutputFibers.size()-1;
                }

                if(currentFiber != NULL)
                {
                    m_OutputFibers[m_OutputIndicesOfLabels[label-1]]->AddInput(currentFiber);
                }

                mutex.Unlock();
            }
        }

        // Update progress
        mutex.Lock();
        this->SetProgress(this->GetProgress() + m_ProgressStep);
        this->InvokeEvent(itk::ProgressEvent());
        mutex.Unlock();
    } // for each seed in region
}

//----------------------------------------------------------------------------------------

void TractographyAlgorithm::ResampleLabelImage()
{
    LabelResampler::Pointer resampler = LabelResampler::New();
    resampler->SetInput(m_RegionsOfInterest);
    resampler->SetInterpolator((LabelInterpolator::New()).GetPointer());

    Self::LabelImage::SpacingType spacing = m_RegionsOfInterest->GetSpacing();
    Self::LabelImage::SizeType       size = m_RegionsOfInterest->GetLargestPossibleRegion().GetSize();
    size[0] *= spacing[0]/m_SeedSpacing; size[1] *= spacing[1]/m_SeedSpacing; size[2] *= spacing[2]/m_SeedSpacing;
    spacing[0] = m_SeedSpacing; spacing[1] = m_SeedSpacing; spacing[2] = m_SeedSpacing;

    resampler->SetOutputSpacing(spacing);
    resampler->SetSize(size);
    resampler->SetOutputDirection(m_RegionsOfInterest->GetDirection());
    resampler->SetOutputOrigin(m_RegionsOfInterest->GetOrigin());

    resampler->Update();

    m_RegionsOfInterest = resampler->GetOutput();
}

//----------------------------------------------------------------------------------------

unsigned int TractographyAlgorithm::SplitRequestedRegion(unsigned int i, unsigned int num, LabelImage::RegionType &splitRegion)
{
    // Get the output pointer
    LabelImage *labelPtr = m_RegionsOfInterest.GetPointer();

    const LabelImage::SizeType &requestedRegionSize = labelPtr->GetRequestedRegion().GetSize();

    int splitAxis;
    LabelImage::IndexType splitIndex;
    LabelImage::SizeType splitSize;

    // Initialize the splitRegion to the output requested region
    splitRegion = labelPtr->GetRequestedRegion();
    splitIndex = splitRegion.GetIndex();
    splitSize = splitRegion.GetSize();

    // split on the outermost dimension available
    splitAxis = labelPtr->GetImageDimension() - 1;
    while ( requestedRegionSize[splitAxis] == 1 )
    {
        --splitAxis;
        if ( splitAxis < 0 )
        { // cannot split
            itkDebugMacro("  Cannot Split");
            return 1;
        }
    }

    // determine the actual number of pieces that will be generated
    LabelImage::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
    unsigned int valuesPerThread = itk::Math::Ceil< unsigned int >(range / (double)num);
    unsigned int maxThreadIdUsed = itk::Math::Ceil< unsigned int >(range / (double)valuesPerThread) - 1;

    // Split the region
    if ( i < maxThreadIdUsed )
    {
        splitIndex[splitAxis] += i * valuesPerThread;
        splitSize[splitAxis] = valuesPerThread;
    }
    if ( i == maxThreadIdUsed )
    {
        splitIndex[splitAxis] += i * valuesPerThread;
        // last thread needs to process the "rest" dimension being split
        splitSize[splitAxis] = splitSize[splitAxis] - i * valuesPerThread;
    }

    // set the split region ivars
    splitRegion.SetIndex(splitIndex);
    splitRegion.SetSize(splitSize);

    itkDebugMacro("  Split Piece: " << splitRegion);

    return maxThreadIdUsed + 1;
}

//----------------------------------------------------------------------------------------

std::vector< vtkSmartPointer< vtkAppendPolyData > > TractographyAlgorithm::GetOutputFiber() const
{
    std::vector< vtkSmartPointer< vtkAppendPolyData > > labeledFibers;

    for(unsigned int label = 1; label <= m_OutputIndicesOfLabels.size(); label++)
    {
        // If label is a valid and has been processed
        if(m_OutputIndicesOfLabels[label-1] == -1)
        {
            labeledFibers.push_back(NULL);
        }
        else // m_OutputIndicesOfLabels[label-1] >= 0
        {
            labeledFibers.push_back(m_OutputFibers[m_OutputIndicesOfLabels[label-1]]);
        }
    }

    return labeledFibers;
}

} // namespace btk
