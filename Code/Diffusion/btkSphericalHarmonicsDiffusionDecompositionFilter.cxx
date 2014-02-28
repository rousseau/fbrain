/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 18/07/2012
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

#include "btkSphericalHarmonicsDiffusionDecompositionFilter.h"


// Local includes
#include "btkSphericalHarmonics.h"

// ITK includes
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorImage.h"

// Definitions
typedef itk::ExtractImageFilter< btk::DiffusionSequence,itk::Image< short,3 > > ExtractImageFilter;
typedef itk::VectorImage< short,3 > VectorImage;
typedef itk::ImageRegionIterator< VectorImage > VectorImageIterator;


namespace btk
{

SphericalHarmonicsDiffusionDecompositionFilter::SphericalHarmonicsDiffusionDecompositionFilter() : Superclass()
{
    m_SphericalHarmonicsOrder = 4;
    m_RegularizationParameter = 0.006;
    m_EstimationType          = Self::DIFFUSION_SIGNAL;
}

//----------------------------------------------------------------------------------------

SphericalHarmonicsDiffusionDecompositionFilter::~SphericalHarmonicsDiffusionDecompositionFilter()
{
    // ----
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::SetInput(const DiffusionSequence::Pointer input)
{
    // FIX : We have to do this turnover to overcome the problem of dimension 4 of diffusion sequence.
    // The solution is to go to vector image for diffusion sequence storage.
    m_InputDiffusionSequence = input;
    Self::SetNumberOfRequiredInputs(0);
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::ComputeSphericalHarmonicsBasisMatrix()
{
    // This assume that the sequence has only one B0 image at first.
    btk::DiffusionSequence::GradientTable gradientTable = m_InputDiffusionSequence->GetGradientTable();
    unsigned int             numberOfGradientDirections = gradientTable.size() - 1;

    // Resize the matrix (rows: number of gradient directions, columns: number of SH coefficients).
    m_SphericalHarmonicsBasisMatrix = Self::Matrix(numberOfGradientDirections, m_NumberOfSHCoefficients);

    // Compute the basis
    for(unsigned int u = 0; u < numberOfGradientDirections; u++)
    {
        unsigned int j = 0;

        for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
        {
            for(int m = -(int)l; m <= (int)l; m++)
            {
                m_SphericalHarmonicsBasisMatrix(u,j++) = btk::SphericalHarmonics::ComputeBasis(gradientTable[u+1].GetSphericalDirection(), l, m);
            } // for each m
        } // for each even order
    } // for each gradient direction
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::ComputeRegularizationMatrix()
{
    // Resize the matrix (rows: number of SH coefficients, columns: number of SH coefficients).
    m_RegularizationMatrix = Self::Matrix(m_NumberOfSHCoefficients, m_NumberOfSHCoefficients);

    // This matrix is diagonal, so we fill the other coefficients with 0.
    m_RegularizationMatrix.Fill(0);

    // Compute the diagonal coefficients.
    unsigned int i = 0;

    for(unsigned int l = 0; l <= m_SphericalHarmonicsOrder; l += 2)
    {
        float   lp1 = l+1;
        float value = l*l * lp1*lp1;

        for(int m = -(int)l; m <= (int)l; m++)
        {
            m_RegularizationMatrix(i,i) = value;
            i++;
        } // for each m
    } // for each even order
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::ComputeTransitionMatrix()
{
    // This assume that the sequence has only one B0 image at first.
    btk::DiffusionSequence::GradientTable gradientTable = m_InputDiffusionSequence->GetGradientTable();
    unsigned int             numberOfGradientDirections = gradientTable.size() - 1;

    // Resize  the matrix (rows: number of SH coefficients, columns: number of gradient directions).
    m_TransitionMatrix = Matrix(m_NumberOfSHCoefficients, numberOfGradientDirections);

    // Compute the transpose of the spherical harmonics basis matrix.
    Self::Matrix SphericalHarmonicsBasisTranspose;
    SphericalHarmonicsBasisTranspose = m_SphericalHarmonicsBasisMatrix.GetTranspose();

    // Compute the inverse.
    Self::Matrix Inverse;
    Inverse = ((SphericalHarmonicsBasisTranspose * m_SphericalHarmonicsBasisMatrix) + (m_RegularizationMatrix * m_RegularizationParameter)).GetInverse();

    // Compute the transition matrix.
    m_TransitionMatrix = Inverse * SphericalHarmonicsBasisTranspose;
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::AllocateOutputs()
{
    // Compute the number of coefficients
    m_NumberOfSHCoefficients = 0.5 * (m_SphericalHarmonicsOrder+1) * (m_SphericalHarmonicsOrder+2);

    // Allocate memory for output image
    btk::DiffusionSequence::SizeType           sequenceSize = m_InputDiffusionSequence->GetLargestPossibleRegion().GetSize();
    btk::DiffusionSequence::SpacingType     sequenceSpacing = m_InputDiffusionSequence->GetSpacing();
    btk::DiffusionSequence::PointType        sequenceOrigin = m_InputDiffusionSequence->GetOrigin();
    btk::DiffusionSequence::DirectionType sequenceDirection = m_InputDiffusionSequence->GetDirection();

    Self::OutputImageType::SizeType outputSize;
    outputSize[0] = sequenceSize[0]; outputSize[1] = sequenceSize[1]; outputSize[2] = sequenceSize[2];

    Self::OutputImageType::SpacingType outputSpacing;
    outputSpacing[0] = sequenceSpacing[0]; outputSpacing[1] = sequenceSpacing[1]; outputSpacing[2] = sequenceSpacing[2];

    Self::OutputImageType::PointType outputOrigin;
    outputOrigin[0] = sequenceOrigin[0]; outputOrigin[1] = sequenceOrigin[1]; outputOrigin[2] = sequenceOrigin[2];

    Self::OutputImageType::DirectionType outputDirection;
    outputDirection(0,0) = sequenceDirection(0,0); outputDirection(0,1) = sequenceDirection(0,1); outputDirection(0,2) = sequenceDirection(0,2);
    outputDirection(1,0) = sequenceDirection(1,0); outputDirection(1,1) = sequenceDirection(1,1); outputDirection(1,2) = sequenceDirection(1,2);
    outputDirection(2,0) = sequenceDirection(2,0); outputDirection(2,1) = sequenceDirection(2,1); outputDirection(2,2) = sequenceDirection(2,2);

    Self::OutputImageType::Pointer output = Self::GetOutput();
    output->SetRegions(outputSize);
    output->SetSpacing(outputSpacing);
    output->SetOrigin(outputOrigin);
    output->SetDirection(outputDirection);
    output->SetVectorLength(m_NumberOfSHCoefficients);
    output->Allocate();
    output->FillBuffer(Self::OutputImagePixelType(0));
}

//----------------------------------------------------------------------------------------

//void SphericalHarmonicsDiffusionDecompositionFilter::GenerateData()
void SphericalHarmonicsDiffusionDecompositionFilter::BeforeThreadedGenerateData()
{
    // Compute needed matrices
    this->ComputeSphericalHarmonicsBasisMatrix();
    this->ComputeRegularizationMatrix();
    this->ComputeTransitionMatrix();


    // Extract reference and gradient images and give it to filter.
    // This process assume that the given sequence is normalized (one reference image at first).
    std::vector< btk::GradientDirection > gradientTable = m_InputDiffusionSequence->GetGradientTable();
    btk::DiffusionSequence::RegionType           region = m_InputDiffusionSequence->GetLargestPossibleRegion();
    ExtractImageFilter::Pointer                 extract = ExtractImageFilter::New();

    // Get the reference image
    region.SetSize(3,0);
    extract->SetInput(m_InputDiffusionSequence);
    extract->SetExtractionRegion(region);
    extract->SetDirectionCollapseToSubmatrix();
    extract->Update();

    // Build a vector image from the diffusion sequence
    // TODO : store diffusion sequence as vector image
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
            pixelValue[k] = m_InputDiffusionSequence->GetPixel(sindex);
        }

        it.Set(pixelValue);
    }

    m_VectorImage = vectorImage;
}

//----------------------------------------------------------------------------------------

void SphericalHarmonicsDiffusionDecompositionFilter::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{
    // Estimate model
    VectorImageIterator                                  it(m_VectorImage, outputRegionForThread);
    itk::ImageRegionIterator< Self::OutputImageType > outIt(this->GetOutput(), outputRegionForThread);

    unsigned int gradienTableSize = m_VectorImage->GetVectorLength();

    for(it.GoToBegin(), outIt.GoToBegin(); !it.IsAtEnd() && !outIt.IsAtEnd(); ++it, ++outIt)
    {
        VectorImage::PixelType signal = it.Get();

        // Get signal matrix
        Self::Matrix signalMatrix(gradienTableSize-1, 1);

        for(unsigned int i = 0; i < signalMatrix.Rows(); i++)
        {
            if(m_EstimationType == Self::APPARENT_DIFFUSION_PROFILE)
            {
                double b0 = static_cast< float >(signal[0]);

                if(!btkFloatingEqual(b0, 0.0))
                {
                    signalMatrix(i,0) = -(1.0/m_InputDiffusionSequence->GetBValues()[i+1]) * std::log(static_cast< float >(signal[i+1]) / b0);
                }
                else
                {
                    signalMatrix(i,0) = 0.0;
                }
            }
            else // m_EstimationType == Self::DIFFUSION_SIGNAL
            {
                double b0 = static_cast< float >(signal[0]);

                if(!btkFloatingEqual(b0, 0.0))
                {
                    signalMatrix(i,0) = static_cast< float >(signal[i+1]) / static_cast< float >(signal[0]);
                }
                else
                {
                    signalMatrix(i,0) = 0.0;
                }
            }
        }

        // Compute coefficients matrix
        Self::Matrix coefficientsMatrix(m_NumberOfSHCoefficients, 1);
        coefficientsMatrix = m_TransitionMatrix * signalMatrix;

        // Set coefficients in pixel
        Self::OutputImagePixelType coefficients(m_NumberOfSHCoefficients);

        for(unsigned int i = 0; i < coefficientsMatrix.Rows(); i++)
        {
            coefficients[i] = coefficientsMatrix(i,0);
        }

        outIt.Set(coefficients);
    }
}

} // namespace btk
