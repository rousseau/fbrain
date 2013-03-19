/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 18/11/2010
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

#ifndef __btkResampleLabelsByInjection_txx
#define __btkResampleLabelsByInjection_txx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "btkResampleLabelsByInjectionFilter.h"
#include "itkObjectFactory.h"
#include "itkIdentityTransform.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkSpecialCoordinatesImage.h"
#include "itkGaussianSpatialFunction.h"

#include "vnl/vnl_inverse.h"

namespace btk
{

/**
 * Initialize new instance
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
ResampleLabelsByInjectionFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::ResampleLabelsByInjectionFilter()
{
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);
  m_OutputDirection.SetIdentity();

  m_UseReferenceImage = false;

  m_Size.Fill( 0 );
  m_OutputStartIndex.Fill( 0 );

  m_DefaultPixelValue = 0;
  m_NumberOfClasses = 6; // 0=Background, 1=CSF, 2=Cortex, 3=White Matter, 4=Brain stem, 5=Cerebellum

}


/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
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
//  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "UseReferenceImage: " << (m_UseReferenceImage ? "On" : "Off") << std::endl;
  return;
}

/**
 * Set the output image spacing.
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputSpacing(
  const double* spacing)
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}


/**
 * Set the output image origin.
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputOrigin(
  const double* origin)
{
  OriginPointType p(origin);
  this->SetOutputOrigin( p );
}


template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GenerateData()
{
  // Get the output pointers
  OutputImagePointer      outputPtr = this->GetOutput();

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
  outputPtr -> FillBuffer(0.0);

  outputPtr -> SetOrigin( referenceImage -> GetOrigin() );
  outputPtr -> SetSpacing( referenceImage -> GetSpacing() );
  outputPtr -> SetDirection( referenceImage -> GetDirection() );

  m_OutputSpacing = referenceImage -> GetSpacing();

  // Image of weights
  wtImage.resize( m_NumberOfClasses + 1 );

  unsigned int im;

  // FIXME The last image ( wtImage[m_NumberOfClasses] ) is used to put the sum
  // (for normalization)

  for(unsigned int c = 0; c <= m_NumberOfClasses; c++)
  {
    wtImage[c] = OutputImageType::New();
    wtImage[c] -> SetRegions( outputRegion );
    wtImage[c] -> Allocate();
    wtImage[c] -> FillBuffer(0.0);

    wtImage[c] -> SetOrigin( referenceImage -> GetOrigin() );
    wtImage[c] -> SetSpacing( referenceImage -> GetSpacing() );
    wtImage[c] -> SetDirection( referenceImage -> GetDirection() );

  }

  // Create gaussian (psf)

  typedef GaussianSpatialFunction< double, ImageDimension, PointType> GaussianFunctionType;
  typename GaussianFunctionType::Pointer gaussian = GaussianFunctionType::New();

  gaussian -> SetNormalized(false);

  typedef FixedArray<double,ImageDimension> ArrayType;

  ArrayType mean;
  mean[0] = 0; mean[1] = 0; mean[2] = 0;
  gaussian -> SetMean( mean  );

  ArrayType sigma;
  double cst = 2*sqrt(2*log(2.0));

  #pragma omp parallel for private(im) schedule(dynamic)

  for(im = 0; im < m_ImageArray.size(); im++)
  {
    // ijk directions for gaussian orientation
    // TODO We should create a nex class, sth like oriented gaussian (orientations given by a reference
    // image) to put all these calculations togheter and make the code more clear

    DirectionType inputDirection = m_ImageArray[im] -> GetDirection();

    PointType idir;
    idir[0] =  inputDirection(0,0); idir[1] =  inputDirection(1,0); idir[2] =  inputDirection(2,0);

    PointType jdir;
    jdir[0] =  inputDirection(0,1); jdir[1] =  inputDirection(1,1); jdir[2] =  inputDirection(2,1);

    PointType kdir;
    kdir[0] =  inputDirection(0,2); kdir[1] =  inputDirection(1,2); kdir[2] =  inputDirection(2,2);

    // Neighborhood iterators on output and weights

    typename NeighborhoodIteratorType::RadiusType radius;

    SpacingType inputSpacing = m_ImageArray[im] -> GetSpacing();

    radius[0] = ceil(inputSpacing[2] / m_OutputSpacing[0]);
    radius[1] = ceil(inputSpacing[2] / m_OutputSpacing[1]);
    radius[2] = ceil(inputSpacing[2] / m_OutputSpacing[2]);

    NeighborhoodIteratorType nbIt( radius, outputPtr, outputPtr -> GetLargestPossibleRegion() );
    nbIt.NeedToUseBoundaryConditionOff();

    NeighborhoodIteratorType nbWtIt( radius, wtImage[m_Labels[im]], wtImage[m_Labels[im]] -> GetLargestPossibleRegion() );
    nbWtIt.NeedToUseBoundaryConditionOff();

    // Division of outputImage into regions (for iteration over neighborhoods)
    FaceCalculatorType faceCalculator;
    FaceListType faceList;

    faceList = faceCalculator( outputPtr, outputPtr -> GetLargestPossibleRegion(), radius);
    typename FaceCalculatorType::FaceListType::iterator fit;
    fit=faceList.begin();

    // Change Gaussian parameters (in case of inputs with different spaces)
    sigma[0] = inputSpacing[0]/cst;
    sigma[1] = inputSpacing[1]/cst;
    sigma[2] = inputSpacing[2]/cst;

    gaussian -> SetSigma( sigma );

    IndexType inputIndex = m_InputImageRegion[im].GetIndex();
    SizeType  inputSize  = m_InputImageRegion[im].GetSize();

    for ( unsigned int i=inputIndex[2]; i < inputIndex[2] + inputSize[2]; i++ )
    {
      // Extract rotation from affine metric to orient the PDF
      VnlMatrixType NQd = m_Transform[im] -> GetSliceTransform(i) -> GetMatrix().GetVnlMatrix();

      VnlVectorType idirTransformed = NQd*idir.GetVnlVector();
      VnlVectorType jdirTransformed = NQd*jdir.GetVnlVector();
      VnlVectorType kdirTransformed = NQd*kdir.GetVnlVector();

      InputImageRegionType wholeSliceRegion;
 //     wholeSliceRegion = m_ImageArray[im] -> GetLargestPossibleRegion();
      wholeSliceRegion = m_InputImageRegion[im];

      IndexType  wholeSliceRegionIndex = wholeSliceRegion.GetIndex();
      SizeType   wholeSliceRegionSize = wholeSliceRegion.GetSize();

      wholeSliceRegionIndex[2]= i;
      wholeSliceRegionSize[2] = 1;

      wholeSliceRegion.SetIndex(wholeSliceRegionIndex);
      wholeSliceRegion.SetSize(wholeSliceRegionSize);

      ConstIteratorType fixedIt( m_ImageArray[im], wholeSliceRegion);

      IndexType fixedIndex;
      IndexType outputIndex;
      IndexType nbIndex;
      PointType physicalPoint;
      PointType nbPoint;
      PointType diffPoint;
      PointType rotPoint;
      PointType transformedPoint;

      double value;

      for(fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt)
      {
        fixedIndex = fixedIt.GetIndex();
        m_ImageArray[im] -> TransformIndexToPhysicalPoint( fixedIndex, physicalPoint );

        transformedPoint = m_Transform[im] -> TransformPoint( physicalPoint);
        outputPtr -> TransformPhysicalPointToIndex( transformedPoint, outputIndex);

        nbIt.SetLocation(outputIndex);
        nbWtIt.SetLocation(outputIndex);

        if ( (*fit).IsInside(outputIndex) )
        {
          for(unsigned int i=0; i<nbIt.Size(); i++)
          {
            nbIndex = nbIt.GetIndex(i);

            outputPtr -> TransformIndexToPhysicalPoint( nbIndex, nbPoint );

            VnlVectorType diffPoint = nbPoint.GetVnlVector() - transformedPoint.GetVnlVector();
            rotPoint[0] = dot_product(diffPoint,idirTransformed);
            rotPoint[1] = dot_product(diffPoint,jdirTransformed);
            rotPoint[2] = dot_product(diffPoint,kdirTransformed);

            value = gaussian->Evaluate( rotPoint );

            nbIt.SetPixel(i, fixedIt.Get() *value + nbIt.GetPixel(i) );
            nbWtIt.SetPixel(i, fixedIt.Get() *value + nbWtIt.GetPixel(i) );

          }
        }
        else
        {
          for(unsigned int i=0; i<nbIt.Size(); i++)
          {
            nbIndex = nbIt.GetIndex(i);

            if ( outputRegion.IsInside(nbIndex) )
            {

              outputPtr -> TransformIndexToPhysicalPoint( nbIndex, nbPoint );

              VnlVectorType diffPoint = nbPoint.GetVnlVector() - transformedPoint.GetVnlVector();
              rotPoint[0] = dot_product(diffPoint,idirTransformed);
              rotPoint[1] = dot_product(diffPoint,jdirTransformed);
              rotPoint[2] = dot_product(diffPoint,kdirTransformed);

              value = gaussian->Evaluate( rotPoint );

              nbIt.SetPixel(i, fixedIt.Get() *value + nbIt.GetPixel(i) );
              nbWtIt.SetPixel(i, fixedIt.Get() *value + nbWtIt.GetPixel(i) );

            }

          }

        }


      }

    }

  }

  // Normalization

  IteratorType sumIt(wtImage[m_NumberOfClasses],wtImage[m_NumberOfClasses] -> GetLargestPossibleRegion() );

  double value;

  for( sumIt.GoToBegin(); !sumIt.IsAtEnd(); ++sumIt  )
  {
    value = 0.0;
    for (unsigned int c = 0; c < m_NumberOfClasses; c++)
    {
      value = value + wtImage[c] -> GetPixel( sumIt.GetIndex() );
    }

    sumIt.Set(value);

  }

  for (unsigned int c = 0; c < m_NumberOfClasses; c++)
  {
    IteratorType wtIt(wtImage[c], wtImage[c] -> GetLargestPossibleRegion() );

    for(wtIt.GoToBegin(),sumIt.GoToBegin(); !wtIt.IsAtEnd(); ++wtIt, ++sumIt)
    {
      if (sumIt.Get() > 0.0)
        wtIt.Set( wtIt.Get()/sumIt.Get() );
      else if (c==0)
        wtIt.Set( 1.0 );
    }


  }


  outputPtr -> FillBuffer(0.0);
  IteratorType outputIt(outputPtr,outputRegion);

  for( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt  )
  {
    value = wtImage[0] -> GetPixel( outputIt.GetIndex() );
    outputIt.Set(0);

    for (unsigned int c = 1; c < m_NumberOfClasses; c++)
    {

     if ( wtImage[c] -> GetPixel( outputIt.GetIndex() ) > value)
     {
       value = wtImage[c] -> GetPixel( outputIt.GetIndex() );
       outputIt.Set( c );
     }
    }

  }


  return;
}


/**
 * Inform pipeline of necessary input image region
 *
 * Determining the actual input region is non-trivial, especially
 * when we cannot assume anything about the transform being used.
 * So we do the easy thing and request the entire input image.
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
 * Set the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
const typename ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>::OutputImageType *
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
unsigned long
ResampleLabelsByInjectionFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetMTime( void ) const
{
  unsigned long latestTime = Object::GetMTime();

  if( m_Transform.size()!=0 )
    {
      for(unsigned int i=0; i<m_Transform.size(); i++)
      {
        if( latestTime < m_Transform[i]->GetMTime() )
        {
          latestTime = m_Transform[i]->GetMTime();
        }
      }
    }

  return latestTime;
}

} // end namespace btk

#endif
