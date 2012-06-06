/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 02/12/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr)

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

#ifndef __btkSuperResolutionRigidImageFilter_txx
#define __btkSuperResolutionRigidImageFilter_txx

#include "btkSuperResolutionRigidImageFilter.h"

namespace btk
{

/**
 * Initialize new instance
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
SuperResolutionRigidImageFilter<TInputImage, TOutputImage   ,TInterpolatorPrecisionType>
::SuperResolutionRigidImageFilter()
{
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);
  m_OutputDirection.SetIdentity();

  m_UseReferenceImage = false;

  m_Size.Fill( 0 );
  m_OutputStartIndex.Fill( 0 );

  m_DefaultPixelValue = 0;

  m_Iterations = 30;
  m_Lambda = 0.1;
  m_PSF = GAUSSIAN;
}

/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage, TOutputImage   ,TInterpolatorPrecisionType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "DefaultPixelValue: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(m_DefaultPixelValue)
     << std::endl;
  os << indent << "Size: " << m_Size << std::endl;
  os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
  os << indent << "UseReferenceImage: " << (m_UseReferenceImage ? "On" : "Off") << std::endl;
  return;
}

template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::SetOutputSpacing(
  const double* spacing)
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}

template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::AddInput(InputImageType* _arg)
{
  m_ImageArray.push_back(_arg);

  this -> SetInput(_arg);

  // Add transforms for this image
  m_Transform.resize( m_Transform.size() + 1 );
  SizeType _argSize = _arg -> GetLargestPossibleRegion().GetSize();
  m_Transform[m_Transform.size()-1].resize(_argSize[2]);

  // Initialize transforms
  for (unsigned int i=0; i<_argSize[2]; i++)
    m_Transform[m_Transform.size()-1][i] = TransformType::New();
}

template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
typename SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>::IndexType
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::LinearToAbsoluteIndex(
  unsigned int linearIndex, InputImageRegionType region)
{
  IndexType absIndex;

  IndexType start = region.GetIndex();
  SizeType  size  = region.GetSize();
  IndexType diffIndex;

  diffIndex[2] = linearIndex / (size[0]*size[1]);

  diffIndex[1] = linearIndex - diffIndex[2]*size[0]*size[1];
  diffIndex[1] = diffIndex[1] / size[0];

  diffIndex[0] = linearIndex - diffIndex[2]*size[0]*size[1] - diffIndex[1]*size[0];

  absIndex[0] = diffIndex[0] + start[0];
  absIndex[1] = diffIndex[1] + start[1];
  absIndex[2] = diffIndex[2] + start[2];

  return absIndex;
}

template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::OptimizeByLeastSquares()
{
  // Fill x
  m_OutputImageRegion = this -> GetReferenceImage() -> GetLargestPossibleRegion();
  unsigned int ncols = m_OutputImageRegion.GetNumberOfPixels();

  m_x.set_size( ncols );
  OutputIteratorType hrIt( this -> GetReferenceImage(), m_OutputImageRegion );
  unsigned int linearIndex = 0;

  for (hrIt.GoToBegin(); !hrIt.IsAtEnd(); ++hrIt, linearIndex++)
    m_x[linearIndex] = hrIt.Get();

  SizeType size = m_OutputImageRegion.GetSize();
  vnl_vector<int>  x_size(3);
  x_size[0] = size[0]; x_size[1] = size[1]; x_size[2] = size[2];

  // Setup cost function

  LeastSquaresVnlCostFunction<InputImageType> f(m_x.size());

  for(unsigned int im = 0; im < m_ImageArray.size(); im++)
  {
    f.AddImage(m_ImageArray[im]);
    f.AddRegion(m_InputImageRegion[im]);

    if ( m_MaskArray.size() > 0)
      f.AddMask( m_MaskArray[im] );

    for(unsigned int i=0; i<m_Transform[im].size(); i++)
      f.SetTransform(im,i,m_Transform[im][i]);
  }
  f.SetReferenceImage(this -> GetReferenceImage());
  f.SetLambda( m_Lambda );
  f.SetPSF( m_PSF );
  f.Initialize();

  // Setup optimizer

  vnl_conjugate_gradient cg(f);
  cg.set_max_function_evals(m_Iterations);

  // Start minimization

  cg.minimize(m_x);
  cg.diagnose_outcome();

}

template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::SetOutputOrigin(
  const double* origin)
{
  PointType p(origin);
  this->SetOutputOrigin( p );
}


template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::GenerateData()
{
  OptimizeByLeastSquares();

  // Get the output pointers
  OutputImagePointer outputPtr = this->GetOutput();

  // Allocate data
  IndexType outputStart;
  outputStart[0] = 0; outputStart[1] = 0; outputStart[2] = 0;

  const OutputImageType * referenceImage = this->GetReferenceImage();

  SizeType outputSize = referenceImage -> GetLargestPossibleRegion().GetSize();

  OutputImageRegionType outputRegion;
  outputRegion.SetIndex(outputStart);
  outputRegion.SetSize(outputSize);

  outputPtr -> SetRegions(outputRegion);
  outputPtr -> Allocate();
  outputPtr -> FillBuffer(0);

  outputPtr -> SetOrigin( referenceImage -> GetOrigin() );
  outputPtr -> SetSpacing( referenceImage -> GetSpacing() );
  outputPtr -> SetDirection( referenceImage -> GetDirection() );

  IndexType hrIndex;
  IndexType hrStart = m_OutputImageRegion.GetIndex();
  SizeType  hrSize  = m_OutputImageRegion.GetSize();
  IndexType hrDiffIndex;


  for (unsigned int i = 0; i<m_x.size(); i++)
  {
    hrDiffIndex[2] = i / (hrSize[0]*hrSize[1]);

    hrDiffIndex[1] = i - hrDiffIndex[2]*hrSize[0]*hrSize[1];
    hrDiffIndex[1] = hrDiffIndex[1] / hrSize[0];

    hrDiffIndex[0] = i - hrDiffIndex[2]*hrSize[0]*hrSize[1] - hrDiffIndex[1]*hrSize[0];


    hrIndex[0] = hrDiffIndex[0] + hrStart[0];
    hrIndex[1] = hrDiffIndex[1] + hrStart[1];
    hrIndex[2] = hrDiffIndex[2] + hrStart[2];

    outputPtr -> SetPixel(hrIndex, m_x[i] );

  }

  return;
}

// TODO: We are not requiring any image region since we are using several inputs.
// We should check if this creates some problems at level of pipeline execution.
// We should also consider the implementation the class derivating from ProcessObject,
// perhaps it's a more logical choice (it does not assume a single input image)
// Same changes should be applied to ResampleImageByInjection and ResampleLabelByInjection
// classes.

/**
 * Inform pipeline of necessary input image region
 *
 * Determining the actual input region is non-trivial, especially
 * when we cannot assume anything about the transform being used.
 * So we do the easy thing and request the entire input image.
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
    return;
    }

  // get pointers to the input and output
  InputImagePointer  inputPtr  =
    const_cast< TInputImage *>( this->GetInput() );

  // Request the entire input image
  InputImageRegionType inputRegion;
  inputRegion = inputPtr->GetLargestPossibleRegion();
  inputPtr->SetRequestedRegion(inputRegion);

  return;
}

/**
 * Get the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
const typename SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>::OutputImageType *
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::GetReferenceImage() const
{
  Self * surrogate = const_cast< Self * >( this );
  const OutputImageType * referenceImage =
    static_cast<const OutputImageType *>(surrogate->ProcessObject::GetInput(1));
  return referenceImage;
}

/**
 * Set the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::SetReferenceImage( const TOutputImage *image )
{
  itkDebugMacro("setting input ReferenceImage to " << image);
  if( image != static_cast<const TOutputImage *>(this->ProcessObject::GetInput( 1 )) )
    {
    this->ProcessObject::SetNthInput(1, const_cast< TOutputImage *>( image ) );
    this->Modified();
    }
}

/** Helper method to set the output parameters based on this image */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::SetOutputParametersFromImage ( const ImageBaseType * image )
{
  this->SetOutputOrigin ( image->GetOrigin() );
  this->SetOutputSpacing ( image->GetSpacing() );
  this->SetOutputDirection ( image->GetDirection() );
  this->SetOutputStartIndex ( image->GetLargestPossibleRegion().GetIndex() );
  this->SetSize ( image->GetLargestPossibleRegion().GetSize() );
}

/**
 * Inform pipeline of required output region
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
void
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  OutputImagePointer outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }

  const OutputImageType * referenceImage = this->GetReferenceImage();

  // Set the size of the output region
  if( m_UseReferenceImage && referenceImage )
    {
    outputPtr->SetLargestPossibleRegion( referenceImage->GetLargestPossibleRegion() );
    }
  else
    {
    typename TOutputImage::RegionType outputLargestPossibleRegion;
    outputLargestPossibleRegion.SetSize( m_Size );
    outputLargestPossibleRegion.SetIndex( m_OutputStartIndex );
    outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
    }

  // Set spacing and origin
  if (m_UseReferenceImage && referenceImage)
    {
    outputPtr->SetOrigin( referenceImage->GetOrigin() );
    outputPtr->SetSpacing( referenceImage->GetSpacing() );
    outputPtr->SetDirection( referenceImage->GetDirection() );
    }
  else
    {
    outputPtr->SetOrigin( m_OutputOrigin );
    outputPtr->SetSpacing( m_OutputSpacing );
    outputPtr->SetDirection( m_OutputDirection );
    }
  return;
}

/**
 * Verify if any of the components has been modified.
 */
template <class TInputImage, class TOutputImage  , class TInterpolatorPrecisionType>
unsigned long
SuperResolutionRigidImageFilter<TInputImage,TOutputImage   ,TInterpolatorPrecisionType>
::GetMTime( void ) const
{
  unsigned long latestTime = Object::GetMTime();

  if( m_Transform.size()!=0 )
    {
      for(unsigned int i=0; i<m_Transform.size(); i++)
      {
        for(unsigned int j=0; j<m_Transform[i].size(); j++)
        {
          if( latestTime < m_Transform[i][j]->GetMTime() )
          {
            latestTime = m_Transform[i][j]->GetMTime();
          }
        }
      }
    }

  return latestTime;
}

} // end namespace btk

#endif
