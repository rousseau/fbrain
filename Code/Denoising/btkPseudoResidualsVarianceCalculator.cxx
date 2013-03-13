/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 25/01/2011
  Author(s): François Rousseau (rousseau@unistra.fr)
             Julien Pontabry   (pontabry@unistra.fr)
             Marc Schweitzer   (marc.schweitzer@unistra.fr)

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

#include "btkPseudoResidualsVarianceCalculator.h"

// ITK includes
#include "itkImageMaskSpatialObject.h"
#include "itkCastImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionConstIterator.h"

// Local includes
#include "btkImageHelper.h"


namespace btk
{

const double PseudoResidualsVarianceCalculator::m_1over6;
const double PseudoResidualsVarianceCalculator::m_sqrt6over7;

//-------------------------------------------------------------------------------------------------

PseudoResidualsVarianceCalculator::PseudoResidualsVarianceCalculator() : Superclass(), m_InputImage(NULL), m_MaskImage(NULL)
{
    // ----
}

//-------------------------------------------------------------------------------------------------

PseudoResidualsVarianceCalculator::~PseudoResidualsVarianceCalculator()
{
    // ----
}

//-------------------------------------------------------------------------------------------------

void PseudoResidualsVarianceCalculator::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//-------------------------------------------------------------------------------------------------

void PseudoResidualsVarianceCalculator::SetInputImage(DiffusionSignal::Pointer signal){
    m_InputImage = dynamic_cast< InternalImage *>(signal.GetPointer());
}

//-------------------------------------------------------------------------------------------------

void PseudoResidualsVarianceCalculator::Compute()
{
    // Test for input image
    if(m_InputImage.IsNull())
    {
        btkException("Missing input image of pseudo residuals variance calculator !");
    }

    // Define region for calculator
    DiffusionSignal::RegionType region;
    if(m_MaskImage.IsNotNull())
    {
        // Test if image and mask are in same physical space
        if(!btk::ImageHelper< InternalImage,MaskImage >::IsInSamePhysicalSpace(m_InputImage, m_MaskImage))
        {
            btkException("Input image and mask image are not in same physical space !");
        }

        typedef itk::ImageMaskSpatialObject< 3 >                               MaskSpatialObject;
        typedef itk::CastImageFilter< MaskImage,MaskSpatialObject::ImageType > MaskSpatialObjectCaster;

        // Set the region of interest
        MaskSpatialObjectCaster::Pointer objectCaster = MaskSpatialObjectCaster::New();
        objectCaster->SetInput(m_MaskImage);
        objectCaster->Update();
        MaskSpatialObject::Pointer object = MaskSpatialObject::New();
        object->SetImage(objectCaster->GetOutput());
        region = object->GetAxisAlignedBoundingBoxRegion();
    }
    else // m_MaskImage.IsNull()
    {
        region = m_InputImage->GetLargestPossibleRegion();

        m_MaskImage = btk::ImageHelper< InternalImage,MaskImage >::CreateNewImageFromPhysicalSpaceOf(m_InputImage, itk::NumericTraits< MaskImage::PixelType >::OneValue());
    }

    // Set requested region for input image
    m_InputImage->SetRequestedRegion(region);

    // Set requested region for mask image
    m_MaskImage->SetRequestedRegion(region);

    // Initialize vectors
    m_PseudoResidualsVariance     = std::vector< double >(m_InputImage->GetNumberOfComponentsPerPixel(), 0.0);
    m_PseudoResidualsStdDeviation = std::vector< double >(m_InputImage->GetNumberOfComponentsPerPixel(), 0.0);

    // Generate data (compute the variance of the square residuals)
    this->GenerateData();
}

//-------------------------------------------------------------------------------------------------

void PseudoResidualsVarianceCalculator::GenerateData()
{
    typedef itk::ConstNeighborhoodIterator< InternalImage > ConstNeighborhoodImageIterator;
    typedef itk::ImageRegionConstIterator< MaskImage >      ConstImageRegionIterator;

    // Define neighborhood iterator for input image
    ConstNeighborhoodImageIterator::RadiusType radius;
    radius[0] = radius[1] = radius[2] = 1;

    ConstNeighborhoodImageIterator neighborInputImageIterator(radius, m_InputImage, m_InputImage->GetRequestedRegion());

    ConstNeighborhoodImageIterator::OffsetType offset1 = { {-1, 0, 0} };
    ConstNeighborhoodImageIterator::OffsetType offset2 = { { 1, 0, 0} };
    ConstNeighborhoodImageIterator::OffsetType offset3 = { { 0,-1, 0} };
    ConstNeighborhoodImageIterator::OffsetType offset4 = { { 0, 1, 0} };
    ConstNeighborhoodImageIterator::OffsetType offset5 = { { 0, 0,-1} };
    ConstNeighborhoodImageIterator::OffsetType offset6 = { { 0, 0, 1} };

    // Define mask image iterator
    ConstImageRegionIterator maskIterator(m_MaskImage, m_MaskImage->GetRequestedRegion());


    // Loop on both images and compute the sum of square residuals
    unsigned int numberOfVoxels = 0;

    for(neighborInputImageIterator.GoToBegin(), maskIterator.GoToBegin(); !neighborInputImageIterator.IsAtEnd() && !maskIterator.IsAtEnd(); ++neighborInputImageIterator, ++maskIterator)
    {
        if(maskIterator.Get() > 0)
        {
            // Compute residual
            DiffusionSignal::PixelType residual = m_sqrt6over7 * (
                                  neighborInputImageIterator.GetCenterPixel() - m_1over6 * (
                                      neighborInputImageIterator.GetPixel(offset1) + neighborInputImageIterator.GetPixel(offset2) + neighborInputImageIterator.GetPixel(offset3) + neighborInputImageIterator.GetPixel(offset4) + neighborInputImageIterator.GetPixel(offset5) + neighborInputImageIterator.GetPixel(offset6)
                                      )
                                  );

            // Sum of square residuals
            for(unsigned int i = 0; i < m_InputImage->GetNumberOfComponentsPerPixel(); i++)
            {
                m_PseudoResidualsVariance[i] += residual[i]*residual[i];
            }

            // Compute the number of visited voxels
            numberOfVoxels++;
        }
    }

    // Normalize by the number of visited voxels and compute standard deviation
    for(unsigned int i = 0; i < m_InputImage->GetNumberOfComponentsPerPixel(); i++)
    {
        m_PseudoResidualsVariance[i]    /= numberOfVoxels;
        m_PseudoResidualsStdDeviation[i] = std::sqrt(m_PseudoResidualsVariance[i]);
    }
}

} // namespace btk
