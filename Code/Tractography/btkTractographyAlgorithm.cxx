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

// VTK includes
#include "vtkPolyData.h"


typedef itk::ResampleImageFilter< btk::TractographyAlgorithm::LabelImage,btk::TractographyAlgorithm::LabelImage > LabelResampler;
typedef itk::NearestNeighborInterpolateImageFunction< btk::TractographyAlgorithm::LabelImage,double > LabelInterpolator;
typedef itk::ImageRegionIteratorWithIndex< btk::TractographyAlgorithm::LabelImage > LabelIterator;


namespace btk
{

TractographyAlgorithm::TractographyAlgorithm() : m_SeedSpacing(1)
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
    btkCoutMacro("TractographyAlgorithm::Update");

    Self::ResampleLabelImage();

    LabelIterator it(m_RegionsOfInterest, m_RegionsOfInterest->GetLargestPossibleRegion());

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
                // Initialize output
                m_CurrentFiber = NULL;

                // Start tractography from seed point
                PropagateSeed(point);

                // Save data
                for(int l = m_OutputFibers.size(); l < label; l++)
                {
                    m_OutputIndicesOfLabels.push_back(-1);
                }

                if(m_OutputIndicesOfLabels[label-1] == -1)
                {
                    m_OutputFibers.push_back(vtkSmartPointer< vtkAppendPolyData >::New());
                    m_OutputIndicesOfLabels[label-1] = m_OutputFibers.size()-1;
                }

                if(m_CurrentFiber != NULL)
                {
                    m_OutputFibers[m_OutputIndicesOfLabels[label-1]]->AddInput(m_CurrentFiber);
                }
            }
        }
    } // for each seed
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
