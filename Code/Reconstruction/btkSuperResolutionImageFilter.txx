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

#ifndef __SuperResolutionImageFilter_txx
#define __SuperResolutionImageFilter_txx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "btkSuperResolutionImageFilter.h"
#include "itkObjectFactory.h"
#include "itkIdentityTransform.h"
#include "itkProgressReporter.h"
#include "itkNeighborhoodIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkSpecialCoordinatesImage.h"
#include "btkOrientedSpatialFunction.h"

#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matops.h"

namespace btk
{

/**
 * Initialize new instance
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
SuperResolutionImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::SuperResolutionImageFilter()
{
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);
  m_OutputDirection.SetIdentity();

  m_UseReferenceImage = false;

  m_Size.Fill( 0 );
  m_OutputStartIndex.Fill( 0 );

  m_DefaultPixelValue = 0;
  m_SimulatedImagesUpdated = false;

  m_Iterations = 100;
  m_OptimizationMethod = MSE;
}


/**
 * Print out a description of self
 *
 * \todo Add details about this class
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
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
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputSpacing(
  const double* spacing)
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::UpdateSimulatedImages()
{
  unsigned int offset = 0;

  for (unsigned int im = 0; im < m_ImageArray.size(); im++)
  {

    std::cout << "Updating image " << im << std::endl;

    IndexType absIndex;

    IndexType start = m_InputImageRegion[im].GetIndex();
    SizeType  size  = m_InputImageRegion[im].GetSize();
    unsigned int nvoxels = m_InputImageRegion[im].GetNumberOfPixels();
    IndexType diffIndex;

    for( unsigned int i=0; i<nvoxels; i++)
    {

      diffIndex[2] = i / (size[0]*size[1]);

      diffIndex[1] = i - diffIndex[2]*size[0]*size[1];
      diffIndex[1] = diffIndex[1] / size[0];

      diffIndex[0] = i - diffIndex[2]*size[0]*size[1] - diffIndex[1]*size[0];

      absIndex[0] = diffIndex[0] + start[0];
      absIndex[1] = diffIndex[1] + start[1];
      absIndex[2] = diffIndex[2] + start[2];

      m_SimulatedImages[im]->SetPixel(absIndex,m_ysim[i + offset]);

    }

    offset = offset + nvoxels;

  }

}

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
typename SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>::IndexType
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::OptimizeByLeastSquares()
{

  vnl_vector<double>  x;
  x = vnl_matops::f2d(m_x);

  LeastSquaresVnlCostFunction f(x.size());
  f.SetParameters(m_H,m_y,x);

  std::cout << "after setting parameters" << std::endl; std::cout.flush();

  vnl_conjugate_gradient cg(f);

  std::cout << "after conjugate gradient" << std::endl; std::cout.flush();

  cg.minimize(x);

  std::cout << "after minimize" << std::endl; std::cout.flush();

  cg.diagnose_outcome();
}

template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::CreateH()
{

  // FIXME: replace the following code with an user-provided region
  // Just for testing purposes we have set it to the region covered
  // by the axial image mask

  IndexType start_lr = m_InputImageRegion[0].GetIndex();
  SizeType  size_lr  = m_InputImageRegion[0].GetSize();

  std::cout << "LR region = " << std::endl;
  std::cout << m_InputImageRegion[0] << std::endl;

  SpacingType spacing_lr = m_ImageArray[0] -> GetSpacing();
  SpacingType spacing_hr = this -> GetReferenceImage() -> GetSpacing();

  PointType start_lr_pt;
  m_ImageArray[0] -> TransformIndexToPhysicalPoint(start_lr,start_lr_pt);
  start_lr_pt = start_lr_pt - (spacing_lr - spacing_hr) / 2.0;

  IndexType start_hr;
  this -> GetReferenceImage() -> TransformPhysicalPointToIndex(start_lr_pt, start_hr);

  SizeType  size_hr;
  size_hr[0] = (int)(size_lr[0] * spacing_lr[0] / spacing_hr[0]);
  size_hr[1] = (int)(size_lr[1] * spacing_lr[1] / spacing_hr[1]);
  size_hr[2] = (int)(size_lr[2] * spacing_lr[2] / spacing_hr[2]);

  IndexType end_hr;
  end_hr[0] = start_hr[0] + size_hr[0] - 1 ;
  end_hr[1] = start_hr[1] + size_hr[1] - 1 ;
  end_hr[2] = start_hr[2] + size_hr[2] - 1 ;

  m_OutputImageRegion.SetIndex( start_hr );
  m_OutputImageRegion.SetSize(  size_hr );

  // matrix resizing

  unsigned int ncols = m_OutputImageRegion.GetNumberOfPixels();

//  m_ysim.resize(m_ImageArray.size());

  m_x.set_size( ncols );

  std::cout << "Number of voxels SR image = " << m_x.size() << std::endl;

  // Fills x vector, since it does not change during H contruction
  OutputIteratorType hrIt( this -> GetReferenceImage(), m_OutputImageRegion );

  unsigned int linearIndex = 0;

  for (hrIt.GoToBegin(); !hrIt.IsAtEnd(); ++hrIt, linearIndex++)
    m_x[linearIndex] = hrIt.Get();

  // Interpolator for HR image

  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator -> SetInputImage( this -> GetReferenceImage() );

  // Differential continuous indexes to perform the neighborhood iteration

  std::vector<ContinuousIndexType> deltaIndexes;
  int npoints =  spacing_lr[2] / (2.0 * spacing_hr[2]) ;
  ContinuousIndexType delta;
  delta[0] = 0.0; delta[1] = 0.0;

  for (int i = -npoints ; i <= npoints; i++ )
  {
    delta[2] = i * 0.5 / (double) npoints;
    deltaIndexes.push_back(delta);
    std::cout << delta << std::endl;
  }

  // Set size of matrices

  unsigned int nrows = 0;
  for(unsigned int im = 0; im < m_ImageArray.size(); im++)
    nrows += m_InputImageRegion[im].GetNumberOfPixels();

  std::cout << "rows in h = " << nrows << std::endl;

  m_H.set_size(nrows, ncols);
//  m_Ht.set_size(ncols, nrows);
//  m_Hbp.set_size(ncols,nrows);
  m_y.set_size(nrows);
//  m_ysim[im].set_size(nrows);

  // H is different for each input image
  unsigned int offset = 0;

  for(unsigned int im = 0; im < m_ImageArray.size(); im++)
  {
    std::cout << "Creating H matrix for image " << im << std::endl; std::cout.flush();

    // Neighborhood iterator

    typename NeighborhoodIteratorType::RadiusType radius;

    SpacingType inputSpacing = m_ImageArray[im] -> GetSpacing();

    // FIXME: the radius should depend on the specific PSF
    // One posibility is to create a member function of the
    // PSF to provide the support region of the function ...
    radius[0] = 1; radius[1] = 1; radius[2] = 1;

    NeighborhoodIteratorType nbIt( radius, m_ImageArray[im], m_ImageArray[im] -> GetLargestPossibleRegion() );
    nbIt.NeedToUseBoundaryConditionOff();


    // PSF definition

    typename FunctionType::Pointer function = FunctionType::New();
    function -> SetDirection( m_ImageArray[im]  -> GetDirection() );
    function -> SetSpacing( m_ImageArray[im] -> GetSpacing() );

    // Iteration over slices

    IndexType inputIndex = m_InputImageRegion[im].GetIndex();
    SizeType  inputSize  = m_InputImageRegion[im].GetSize();

    IndexType lrIndex;
    IndexType lrDiffIndex;
    unsigned int lrLinearIndex;

    IndexType hrIndex;
    IndexType hrDiffIndex;
    ContinuousIndexType hrContIndex;
    unsigned int hrLinearIndex;

    ContinuousIndexType nbIndex;

    PointType lrPoint;
    PointType nbPoint;
    PointType transformedPoint;

    for ( unsigned int i=inputIndex[2]; i < inputIndex[2] + inputSize[2]; i++ )
    {

      InputImageRegionType wholeSliceRegion;
      wholeSliceRegion = m_InputImageRegion[im];

      IndexType  wholeSliceRegionIndex = wholeSliceRegion.GetIndex();
      SizeType   wholeSliceRegionSize  = wholeSliceRegion.GetSize();

      wholeSliceRegionIndex[2]= i;
      wholeSliceRegionSize[2] = 1;

      wholeSliceRegion.SetIndex(wholeSliceRegionIndex);
      wholeSliceRegion.SetSize(wholeSliceRegionSize);

      ConstIteratorType fixedIt( m_ImageArray[im], wholeSliceRegion);

      double lrValue;
      double hrValue;

      for(fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt)
      {
        lrIndex = fixedIt.GetIndex();

        lrDiffIndex[0] = lrIndex[0] - inputIndex[0];
        lrDiffIndex[1] = lrIndex[1] - inputIndex[1];
        lrDiffIndex[2] = lrIndex[2] - inputIndex[2];

        lrLinearIndex = lrDiffIndex[0] + lrDiffIndex[1]*inputSize[0] + lrDiffIndex[2]*inputSize[0]*inputSize[1];

        m_y[lrLinearIndex + offset] = fixedIt.Get();

        m_ImageArray[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );
        nbIt.SetLocation( lrIndex);
        function -> SetCenter( lrPoint );

        for(unsigned int k=0; k<deltaIndexes.size(); k++)
        {
          nbIndex[0] = deltaIndexes[k][0] + lrIndex[0];
          nbIndex[1] = deltaIndexes[k][1] + lrIndex[1];
          nbIndex[2] = deltaIndexes[k][2] + lrIndex[2];

          m_ImageArray[im] -> TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );
          lrValue = function -> Evaluate(nbPoint);

          if ( lrValue > 0)
          {

// FIXME Exchange lines after testing with simulated images
//          transformedPoint = m_Transform[im][i] -> TransformPoint( nbPoint);
            transformedPoint = nbPoint;

            this->GetReferenceImage() -> TransformPhysicalPointToContinuousIndex( transformedPoint, hrContIndex );

            bool isInsideHR = true;

            if ( (hrContIndex[0] < start_hr[0]) || (hrContIndex[0] > end_hr[0]) ||
                 (hrContIndex[1] < start_hr[1]) || (hrContIndex[1] > end_hr[1]) ||
                 (hrContIndex[2] < start_hr[2]) || (hrContIndex[2] > end_hr[2]) )
               isInsideHR = false;

            if ( isInsideHR )
            {
              hrValue = interpolator -> Evaluate( transformedPoint );

              for(unsigned int n=0; n<interpolator -> GetContributingNeighbors(); n++)
              {
                hrIndex = interpolator -> GetIndex(n);

                hrDiffIndex[0] = hrIndex[0] - start_hr[0];
                hrDiffIndex[1] = hrIndex[1] - start_hr[1];
                hrDiffIndex[2] = hrIndex[2] - start_hr[2];

                hrLinearIndex = hrDiffIndex[0] + hrDiffIndex[1]*size_hr[0] + hrDiffIndex[2]*size_hr[0]*size_hr[1];
                m_H(lrLinearIndex + offset,hrLinearIndex) = interpolator -> GetOverlap(n)* lrValue;
//                m_Hbp(hrLinearIndex,lrLinearIndex + offset) = lrValue;

              }

            }

          }

        }

      }

    }

    offset += m_InputImageRegion[im].GetNumberOfPixels();

  }

  // Normalize H

  for (unsigned int i = 0; i < m_H.rows(); i++)
  {
    double sum = m_H.sum_row(i);

    VnlSparseMatrixType::row & r = m_H.get_row(i);
    VnlSparseMatrixType::row::iterator col_iter = r.begin();

    for ( ;col_iter != r.end(); ++col_iter)
      (*col_iter).second = (*col_iter).second / sum;
  }

  // Create simulated images

  m_H.mult(m_x,m_ysim);
  UpdateSimulatedImages();
}


/**
 * Set the output image origin.
 */
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputOrigin(
  const double* origin)
{
  OriginPointType p(origin);
  this->SetOutputOrigin( p );
}


template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GenerateData()
{
  CreateH();
  OptimizeByLeastSquares();
  UpdateSimulatedImages();

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
  SizeType hrSize  = m_OutputImageRegion.GetSize();
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

//    std::cout << hrIndex << std::endl;

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
template <class TInputImage, class TOutputImage, class TInterpolatorPrecisionType>
void
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
const typename SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>::OutputImageType *
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
SuperResolutionImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetMTime( void ) const
{
  unsigned long latestTime = Object::GetMTime();

// FIXME Uncomment after testing with simulated images
/*  if( m_Transform.size()!=0 )
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
*/
  return latestTime;
}

} // end namespace itk

#endif
