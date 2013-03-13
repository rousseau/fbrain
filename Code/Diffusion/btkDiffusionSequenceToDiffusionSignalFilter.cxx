/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 11/03/2013
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

#include "btkDiffusionSequenceToDiffusionSignalFilter.h"


// ITK includes
#include "itkImageRegionIterator.h"


namespace btk
{

DiffusionSequenceToDiffusionSignalFilter::DiffusionSequenceToDiffusionSignalFilter() : Superclass()
{
    // ----
}

//----------------------------------------------------------------------------------------

DiffusionSequenceToDiffusionSignalFilter::~DiffusionSequenceToDiffusionSignalFilter()
{
    // ----
}

//----------------------------------------------------------------------------------------

void DiffusionSequenceToDiffusionSignalFilter::SetInput(const DiffusionSequence *input)
{
    // FIX : We have to do this turnover to overcome the problem of dimension 4 of diffusion sequence.
    // The solution is to go to vector image for diffusion sequence storage.
    m_InputDiffusionSequence = input;
    Self::SetNumberOfRequiredInputs(0);
}

//----------------------------------------------------------------------------------------

void DiffusionSequenceToDiffusionSignalFilter::AllocateOutputs()
{
    // Get informations of input image
    InputImageType::SizeType           sequenceSize = m_InputDiffusionSequence->GetLargestPossibleRegion().GetSize();
    InputImageType::SpacingType     sequenceSpacing = m_InputDiffusionSequence->GetSpacing();
    InputImageType::PointType        sequenceOrigin = m_InputDiffusionSequence->GetOrigin();
    InputImageType::DirectionType sequenceDirection = m_InputDiffusionSequence->GetDirection();

    // Define informations for output image
    OutputImageType::SizeType outputSize;
    outputSize[0] = sequenceSize[0]; outputSize[1] = sequenceSize[1]; outputSize[2] = sequenceSize[2];

    OutputImageType::SpacingType outputSpacing;
    outputSpacing[0] = sequenceSpacing[0]; outputSpacing[1] = sequenceSpacing[1]; outputSpacing[2] = sequenceSpacing[2];

    OutputImageType::PointType outputOrigin;
    outputOrigin[0] = sequenceOrigin[0]; outputOrigin[1] = sequenceOrigin[1]; outputOrigin[2] = sequenceOrigin[2];

    OutputImageType::DirectionType outputDirection;
    outputDirection(0,0) = sequenceDirection(0,0); outputDirection(0,1) = sequenceDirection(0,1); outputDirection(0,2) = sequenceDirection(0,2);
    outputDirection(1,0) = sequenceDirection(1,0); outputDirection(1,1) = sequenceDirection(1,1); outputDirection(1,2) = sequenceDirection(1,2);
    outputDirection(2,0) = sequenceDirection(2,0); outputDirection(2,1) = sequenceDirection(2,1); outputDirection(2,2) = sequenceDirection(2,2);

    // Define new gradient table (remove first entry)
    std::vector< GradientDirection > gradientTable = m_InputDiffusionSequence->GetGradientTable();
    gradientTable.erase(gradientTable.begin());

    // Allocate output image
    OutputImageType::Pointer output = this->GetOutput();
    output->SetGradientTable(gradientTable);
    output->SetRegions(outputSize);
    output->SetSpacing(outputSpacing);
    output->SetOrigin(outputOrigin);
    output->SetDirection(outputDirection);
    output->SetVectorLength(gradientTable.size());
    output->Allocate();
    output->FillBuffer(OutputImagePixelType(0));
}

//----------------------------------------------------------------------------------------

void DiffusionSequenceToDiffusionSignalFilter::GenerateData()
{
    // Allocate outputs
    this->AllocateOutputs();


    typedef itk::ImageRegionIterator< OutputImageType > OutputImageIterator;

    // Extract reference and gradient images and give it to filter.
    // This process assume that the given sequence is normalized (one reference image at first).

    // Build a vector image from the diffusion sequence
    OutputImageType::Pointer output = this->GetOutput();
    OutputImageIterator it(output, output->GetLargestPossibleRegion());

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        OutputImageType::IndexType vindex = it.GetIndex();
        InputImageType::IndexType sindex;
        sindex[0] = vindex[0]; sindex[1] = vindex[1]; sindex[2] = vindex[2];

        OutputImageType::PixelType pixelValue;
        pixelValue.SetSize(output->GetVectorLength());

        sindex[3] = 0;
        InputImageType::PixelType B0 = m_InputDiffusionSequence->GetPixel(sindex);

        for(unsigned int k = 0; k < output->GetVectorLength(); k++)
        {
            sindex[3] = k+1;

            if(B0 != 0)
            {
                pixelValue[k] = static_cast< double >(m_InputDiffusionSequence->GetPixel(sindex)) / static_cast< double >(B0);
            }
            else // B0 == 0
            {
                pixelValue[k] = 0.0;
            }
        }

        it.Set(pixelValue);
    }
}

} // namespace btk
