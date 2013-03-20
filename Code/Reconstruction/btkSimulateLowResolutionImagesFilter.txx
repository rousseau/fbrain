/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/03/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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

#ifndef BTK_BTKSIMULATELOWRESOLUTIONIMAGESFILTER_TXX
#define BTK_BTKSIMULATELOWRESOLUTIONIMAGESFILTER_TXX


#include "btkSimulateLowResolutionImagesFilter.hxx"

namespace btk
{
//-------------------------------------------------------------------------------------------------
template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType>
SimulateLowResolutionImagesFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::SimulateLowResolutionImagesFilter()
{
}
//-------------------------------------------------------------------------------------------------
template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType>
void
SimulateLowResolutionImagesFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::Update()
{
    std::cout<<"Compute H..."<<std::flush;
    this->ComputeH();
    std::cout<<" Done !"<<std::endl;

    std::cout<<"Compute Simulated Images..."<<std::flush;
    this->UpdateSimulatedImages();
    std::cout<<" Done !"<<std::endl;
}
//-------------------------------------------------------------------------------------------------
template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType>
void
SimulateLowResolutionImagesFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::ComputeH()
{
    //FIXME : This is an old implementation of SuperResolution, it may be not optimized.
    // Result are not good, when we put motion image is the average of image without motion and image with motion...

    SpacingType spacing_lr = m_Images[0] -> GetSpacing();
    SpacingType spacing_hr = this -> GetReferenceImage() -> GetSpacing();
    IndexType start_hr;
    SizeType  size_hr;

    m_OutputImageRegion = this -> GetReferenceImage() -> GetLargestPossibleRegion();
    start_hr = m_OutputImageRegion.GetIndex();
    size_hr = m_OutputImageRegion.GetSize();

    IndexType end_hr;
    end_hr[0] = start_hr[0] + size_hr[0] - 1 ;
    end_hr[1] = start_hr[1] + size_hr[1] - 1 ;
    end_hr[2] = start_hr[2] + size_hr[2] - 1 ;

    // Fill x

    unsigned int ncols = m_OutputImageRegion.GetNumberOfPixels();

    m_x.set_size( ncols );

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
    }

    // Set size of matrices

    unsigned int nrows = 0;
    for(unsigned int im = 0; im < m_Images.size(); im++)
      nrows += m_InputImageRegion[im].GetNumberOfPixels();

    m_H.set_size(nrows, ncols);
    m_y.set_size(nrows);
    m_y.fill(0.0);

    // H is different for each input image
    unsigned int offset = 0;

    for(unsigned int im = 0; im < m_Images.size(); im++)
    {
      SpacingType inputSpacing = m_Images[im] -> GetSpacing();

      // PSF definition

      typename FunctionType::Pointer function = FunctionType::New();
      function -> SetPSF( m_PSF );
      function -> SetDirection( m_Images[im] -> GetDirection() );
      function -> SetSpacing( m_Images[im] -> GetSpacing() );

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

        ConstIteratorType fixedIt( m_Images[im], wholeSliceRegion);

        double lrValue;
        double hrValue;

        for(fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt)
        {
          lrIndex = fixedIt.GetIndex();
          m_Images[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );

          if ( interpolator -> IsInsideBuffer( lrPoint ) )
          {

            lrDiffIndex[0] = lrIndex[0] - inputIndex[0];
            lrDiffIndex[1] = lrIndex[1] - inputIndex[1];
            lrDiffIndex[2] = lrIndex[2] - inputIndex[2];

            lrLinearIndex = lrDiffIndex[0] + lrDiffIndex[1]*inputSize[0] + lrDiffIndex[2]*inputSize[0]*inputSize[1];

            m_y[lrLinearIndex + offset] = fixedIt.Get();

            m_Images[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );
            function -> SetCenter( lrPoint );

            for(unsigned int k=0; k<deltaIndexes.size(); k++)
            {
              nbIndex[0] = deltaIndexes[k][0] + lrIndex[0];
              nbIndex[1] = deltaIndexes[k][1] + lrIndex[1];
              nbIndex[2] = deltaIndexes[k][2] + lrIndex[2];

              m_Images[im] -> TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );
              lrValue = function -> Evaluate(nbPoint);

              if ( lrValue > 0)
              {
                //transformedPoint = m_Transform[im][i] -> TransformPoint( nbPoint);
                transformedPoint = m_Transforms[im]-> TransformPoint( nbPoint);

                this->GetReferenceImage() -> TransformPhysicalPointToContinuousIndex( transformedPoint, hrContIndex );

                bool isInsideHR = true;

                // FIXME This checking should be done for all points first, and discard the point
                // if al least one point is out of the reference image

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
                     m_H(lrLinearIndex + offset, hrLinearIndex) += interpolator -> GetOverlap(n)* lrValue;


                  }

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


    m_H.mult(m_x,m_ysim);

}
//-------------------------------------------------------------------------------------------------
template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType>
void
SimulateLowResolutionImagesFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::UpdateSimulatedImages()
{
    unsigned int offset = 0;

    for (unsigned int im = 0; im < m_Images.size(); im++)
    {
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
//-------------------------------------------------------------------------------------------------
} // namespace btk

#endif
