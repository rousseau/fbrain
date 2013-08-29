/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 12/07/2012
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

#include "btkDiffusionTensorReconstructionFilter.h"


// ITK includes
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

// Definitions
typedef itk::ExtractImageFilter< btk::DiffusionSequence,itk::Image< short,3 > > ExtractImageFilter;
typedef itk::VectorImage< short,3 > VectorImage;
typedef itk::ImageRegionIterator< VectorImage > VectorImageIterator;


namespace btk
{

DiffusionTensorReconstructionFilter::DiffusionTensorReconstructionFilter() : Superclass()
{
    // ----
}

//----------------------------------------------------------------------------------------

DiffusionTensorReconstructionFilter::~DiffusionTensorReconstructionFilter()
{
    // ----
}

//----------------------------------------------------------------------------------------

void DiffusionTensorReconstructionFilter::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void DiffusionTensorReconstructionFilter::SetInput(DiffusionSequence::Pointer sequence)
{
    this->m_InputDiffusionSequence = sequence;
}

//----------------------------------------------------------------------------------------

void DiffusionTensorReconstructionFilter::Update()
{
    // Extract reference and gradient images and give it to filter.
    // This process assume that the given sequence is normalized (one reference image at first).
    std::vector< btk::GradientDirection > gradientTable = this->m_InputDiffusionSequence->GetGradientTable();
    btk::DiffusionSequence::RegionType          region = this->m_InputDiffusionSequence->GetLargestPossibleRegion();
    ExtractImageFilter::Pointer                extract = ExtractImageFilter::New();

    if(gradientTable.size() != region.GetSize(3))
    {
        throw(std::string("There is a mismatch between the number of gradient directions and gradient images ! Cannot estimate tensors !"));
    }

    if(gradientTable.size() < 6)
    {
        throw(std::string("There are less than 6 gradient directions ! Cannot estimate tensors !"));
    }

    // Get the reference image
    region.SetSize(3,0);
    extract->SetInput(this->m_InputDiffusionSequence);
    extract->SetExtractionRegion(region);
    extract->Update();

    // Build a vector image from the diffusion sequence
    // TODO : store diffusion sequence as vector image
    // But test first if linear interpolation is ok
    VectorImage::Pointer vectorImage = VectorImage::New();
    vectorImage->SetRegions(extract->GetOutput()->GetLargestPossibleRegion());
    vectorImage->SetSpacing(extract->GetOutput()->GetSpacing());
    vectorImage->SetOrigin(extract->GetOutput()->GetOrigin());
    vectorImage->SetDirection(extract->GetOutput()->GetDirection());
    vectorImage->SetVectorLength(gradientTable.size());
    vectorImage->Allocate();
    vectorImage->FillBuffer(itk::VariableLengthVector< short >(0));

    VectorImageIterator it(vectorImage, vectorImage->GetLargestPossibleRegion());

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        VectorImage::IndexType vindex = it.GetIndex();
        btk::DiffusionSequence::IndexType sindex;
        sindex[0] = vindex[0]; sindex[1] = vindex[1]; sindex[2] = vindex[2];

        VectorImage::PixelType pixelValue;
        pixelValue.SetSize(vectorImage->GetVectorLength());

        for(unsigned int k = 0; k < vectorImage->GetVectorLength(); k++)
        {
            sindex[3] = k;
            pixelValue[k] = this->m_InputDiffusionSequence->GetPixel(sindex);
        }

        it.Set(pixelValue);
    }

    // Get the gradient table
    Self::GradientDirectionContainerType::Pointer directions = Self::GradientDirectionContainerType::New();

    for(unsigned int i = 0; i < gradientTable.size(); i++)
    {
        directions->InsertElement(i, gradientTable[i].GetVnlVectorFixed());
    }

    // Set gradient table and images to filter
    this->SetGradientImage(directions, vectorImage);

    // B-value.
    // ITK consider that there is only on b-value for the whole gradient table.
    this->SetBValue(this->m_InputDiffusionSequence->GetBValues()[1]);

    // WARNING: check if there is no more bugs with netlib/dsvdc.c
//    // Needed until netlib/dsvdc.c has been fixed.
//    this->SetNumberOfThreads(1);

    // Evaluate tensors
    Superclass::Update();
}

} // namespace btk
