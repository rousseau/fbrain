/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/03/2010
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

#ifndef __btkSliceToSliceInterpolateImageFunction_txx
#define __btkSliceToSliceInterpolateImageFunction_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

#include "btkSliceToSliceInterpolateImageFunction.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_inverse.h"


namespace btk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep>
const unsigned long
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep>
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::SliceToSliceInterpolateImageFunction()
{
  m_GradientTableCartesian = 0;
  m_GradientTableSpherical = 0;
//  m_RBFInterpolator = 0;
//  m_Sigma = 1;
//  m_kdTreeSphere = 0;
//  m_kdTreeSpace = 0;
  m_TransformsAreSet = false;
  m_Rspa = 1;
}

/**
 * Destructor
 */
template<class TInputImage, class TCoordRep>
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::~SliceToSliceInterpolateImageFunction()
{

  delete m_kdTreeSpace;

  delete []m_dataVals;

  annClose();

}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}

template<class TInputImage, class TCoordRep>
void
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::SetGradientTable( const char* input )
{

  m_GradientTableCartesian.set_size( m_ImageSize[3],3 );
  m_OriginalGradientTable.set_size( m_ImageSize[3],3 );

  // Fill gradient table in cartesian coordinates

  FILE* fr;
  fr = fopen( input, "r" );

  float value;
  unsigned ncol = 1;
  unsigned nrow = 1;

  for (unsigned int j=1; j <= 3*m_ImageSize[3]; j++)
  {
    if ( ( j % m_ImageSize[3] ) == 0 )
    {
      fscanf( fr, "%f\n", &value);
      m_GradientTableCartesian(nrow-1,ncol-1)=value;
      m_OriginalGradientTable(nrow-1,ncol-1)=value;
      nrow = 0;
      ncol++;
    } else
    {
      fscanf( fr, "%f ", &value);
      m_GradientTableCartesian(nrow-1,ncol-1)=value;
      m_OriginalGradientTable(nrow-1,ncol-1)=value;
    }
    nrow++;
  }

  fclose (fr);

}

template<class TInputImage, class TCoordRep>
void
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::RotateGradients()
{

//  std::cout << "original gradient table = "<< std::endl;
//  std::cout << m_GradientTableCartesian << std::endl;

  for (unsigned int i=1; i <= m_NumberOfGradients; i++)
  {

    VnlMatrixType Md = m_TransformArray[i-1] -> GetMatrix().GetVnlMatrix();

    VnlMatrixType PQd = Md;
    VnlMatrixType NQd = Md;
    VnlMatrixType PQNQDiffd;

    for(unsigned int ni = 0; ni < 100; ni++ )
    {
      // Average current Qi with its inverse transpose
      NQd = ( PQd + vnl_inverse_transpose( PQd ) ) / 2.0;
      PQNQDiffd = NQd - PQd;
      if( PQNQDiffd.frobenius_norm() < 1e-7 )
      {
 //       std::cout << "Polar decomposition used "
 //                 << ni << " iterations " << std::endl;
        break;
      }
      else
      {
        PQd = NQd;
      }
    }

//    std::cout << NQd << std::endl;

//    MatrixType QMatrix;
//    QMatrix = NQd;

//   transform -> SetMatrix(QMatrix);
//   RigidTransformType::PointType originalGradient;
//   originalGradient[0] = m_GradientTableCartesian(nrow-1,ncol-1)
//   RigidTransformType::PointType rotatedGradient;

    VnlVectorType grad = m_GradientTableCartesian.get_row(i);
    m_GradientTableCartesian.set_row(i, NQd*grad);

  }

//  std::cout << "rotated gradient tabble = "<< std::endl;
//  std::cout << m_GradientTableCartesian << std::endl;

  // data structure for nearest neighbor search with ANN
  m_dataPtsSphere = annAllocPts(2*m_NumberOfGradients, 3);

  // Conversion to spherical coordinates

  m_GradientTableSpherical.set_size( m_ImageSize[3],2 );

  m_GradientTableSpherical(0,0) = 0.0;
  m_GradientTableSpherical(0,1) = 0.0;

  double x,y,z,r;

  for (unsigned int i=1; i <= m_NumberOfGradients; i++)
  {
    x = m_GradientTableCartesian(i,0);
    y = m_GradientTableCartesian(i,1);
    z = m_GradientTableCartesian(i,2);

    r = std::sqrt(x*x + y*y + z*z);

    // re-normalization to avoid errors due to accuracy in gradients
    x = x/r; y = y/r; z = z/r; r=1;

    m_GradientTableSpherical(i,0) = std::acos(z/r);
    m_GradientTableSpherical(i,1) = std::atan2(y,x);

    // For the nearest neighbor search, we use the negative point too
    // since they represent the same direction

    m_dataPtsSphere[i-1][0] = x;
    m_dataPtsSphere[i-1][1] = y;
    m_dataPtsSphere[i-1][2] = z;

    m_dataPtsSphere[ i-1 + m_NumberOfGradients ][0] = -x;
    m_dataPtsSphere[ i-1 + m_NumberOfGradients ][1] = -y;
    m_dataPtsSphere[ i-1 + m_NumberOfGradients ][2] = -z;

  }

  m_kdTreeSphere = new ANNkd_tree(m_dataPtsSphere, 2*m_NumberOfGradients, 3);

}

template<class TInputImage, class TCoordRep>
void
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::GetGradientDirection( unsigned int index, double &theta, double &phi )
{
  double x,y,z,r;

  x = m_OriginalGradientTable(index,0);
  y = m_OriginalGradientTable(index,1);
  z = m_OriginalGradientTable(index,2);

  r = std::sqrt(x*x + y*y + z*z);

  // re-normalization to avoid errors due to accuracy in gradients
  x = x/r; y = y/r; z = z/r; r=1;

  theta = std::acos(z/r);
  phi = std::atan2(y,x);

}


template<class TInputImage, class TCoordRep>
void
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::SetInputImage(const InputImageType *ptr)
{
  this->Superclass::SetInputImage(ptr);
}

template<class TInputImage, class TCoordRep>
void
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::Initialize( const ImageRegionType & region)
{
  // In case of leave1out
//  std::cout << "Number of piwels in region = " << region.GetNumberOfPixels() << std::endl; std::cout.flush();

  ImageIndexType start = region.GetIndex();
  ImageSizeType  size = region.GetSize();

  ImageRegionType extRegion;

  ImageIndexType  extStart;
  extStart[0]= start[0] - 5;
  extStart[1]= start[1] - 5;
  extStart[2]= start[2] - 1;

  ImageSizeType   extSize;
  extSize[0]= size[0] + 10;
  extSize[1]= size[1] + 10;
  extSize[2]= size[2] + 2;


  extRegion.SetSize(extSize);
  extRegion.SetIndex(extStart);

//  std::cout << "Number of piwels in ext region = " << extRegion.GetNumberOfPixels() << std::endl; std::cout.flush();

  IteratorType extRegionIt( this->GetInputImage(), extRegion );

  IndexType index;
  PointType point;
  PointType transformedPoint;

  m_dataPtsSpace = annAllocPts( extRegion.GetNumberOfPixels(), ImageDimension);
  m_dataVals = new InputPixelType[ extRegion.GetNumberOfPixels() ];

//  std::cout << "ImageDimension = " << ImageDimension << std::endl;
//  std::cout << "reigon = " << extRegion << std::endl;

  unsigned int i = 0;

  for (extRegionIt.GoToBegin(); !extRegionIt.IsAtEnd(); ++extRegionIt)
  {
    index = extRegionIt.GetIndex();

    this->GetInputImage() -> TransformIndexToPhysicalPoint(index,point);
    transformedPoint = m_TransformArray[ index[2] ] -> TransformPoint(point);

    m_dataPtsSpace[i][0] = transformedPoint[0];
    m_dataPtsSpace[i][1] = transformedPoint[1];
    m_dataPtsSpace[i][2] = transformedPoint[2];

    m_dataVals[i] = extRegionIt.Get();

    if ( m_dataVals[i] < 0 )
    {
      std::cout << "Valor menor que 0 !!!! @ " << index << " : " << m_dataVals[i] << " " << extRegionIt.Get() << " " << std::endl; std::cout.flush();

    }

    i++;

  }

  m_kdTreeSpace = new ANNkd_tree(m_dataPtsSpace, extRegion.GetNumberOfPixels(), ImageDimension);

}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType& index) const
{
  unsigned int dim;  // index over dimension

  /**
   * Compute base index = closet index below point
   * Compute distance from point to base index
   */
  signed long baseIndex[ImageDimension];
  double distance[ImageDimension];
  long tIndex;

  for( dim = 0; dim < ImageDimension; dim++ )
    {
    // The following "if" block is equivalent to the following line without
    // having to call floor.
    //    baseIndex[dim] = (long) vcl_floor(index[dim] );
    if (index[dim] >= 0.0)
      {
      baseIndex[dim] = (long) index[dim];
      }
    else
      {
      tIndex = (long) index[dim];
      if (double(tIndex) != index[dim])
        {
        tIndex--;
        }
      baseIndex[dim] = tIndex;
      }
    distance[dim] = index[dim] - double( baseIndex[dim] );
    }

  /**
   * Interpolated value is the weighted sum of each of the surrounding
   * neighbors. The weight for each neighbor is the fraction overlap
   * of the neighbor pixel with respect to a pixel centered on point.
   */
  RealType value = NumericTraits<RealType>::Zero;

  typedef typename NumericTraits<InputPixelType>::ScalarRealType ScalarRealType;
  ScalarRealType totalOverlap = NumericTraits<ScalarRealType>::Zero;

  for( unsigned int counter = 0; counter < m_Neighbors; counter++ )
    {

    double overlap = 1.0;          // fraction overlap
    unsigned int upper = counter;  // each bit indicates upper/lower neighbour
    IndexType neighIndex;

    // get neighbor index and overlap fraction
    for( dim = 0; dim < ImageDimension; dim++ )
      {

      if ( upper & 1 )
        {
        neighIndex[dim] = baseIndex[dim] + 1;
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
        // Take care of the case where the pixel is just
        // in the outer upper boundary of the image grid.
        if( neighIndex[dim] > this->m_EndIndex[dim] )
          {
          neighIndex[dim] = this->m_EndIndex[dim];
          }
#endif
        overlap *= distance[dim];
        }
      else
        {
        neighIndex[dim] = baseIndex[dim];
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
        // Take care of the case where the pixel is just
        // in the outer lower boundary of the image grid.
        if( neighIndex[dim] < this->m_StartIndex[dim] )
          {
          neighIndex[dim] = this->m_StartIndex[dim];
          }
#endif
        overlap *= 1.0 - distance[dim];
        }

      upper >>= 1;

      }

    // get neighbor value only if overlap is not zero
    if( overlap )
      {
      value += static_cast<RealType>( this->GetInputImage()->GetPixel( neighIndex ) ) * overlap;
      totalOverlap += overlap;
      }

    if( totalOverlap == 1.0 )
      {
      // finished
      break;
      }

    }

  return ( static_cast<OutputType>( value ) );
}

/**
 * Evaluate at image index position
 */

template<class TInputImage, class TCoordRep>
typename SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::OutputType
SliceToSliceInterpolateImageFunction< TInputImage, TCoordRep >
::Evaluate( const PointType& point ) const
{

    double r0 = m_Rspa;

    // Search k_spa spatial neighbors

    ANNpoint spatialPt;
    spatialPt = annAllocPt(3);
    spatialPt[0] = point[0];
    spatialPt[1] = point[1];
    spatialPt[2] = point[2];

//    ImageIndexType index;
//    this -> GetInputImage() -> TransformPhysicalPointToIndex(point,index);

//    if ( (index[0]==66) && (index[1]==39) && (index[2]==25))
//    {
//      std::cout << "Punto fisico = " << point << std::endl;
//    }

    // Number of spatial neighbors inside the support

    unsigned int k_spa;
    k_spa = m_kdTreeSpace -> annkFRSearch(spatialPt, 16*m_Rspa*m_Rspa, 0);

    ANNidxArray nnIdx_spa;
    ANNdistArray dists_spa;

    // Nearest neighbors search

    if (k_spa == 0 )
    {
      k_spa = 4;
      nnIdx_spa = new ANNidx[ k_spa ];
      dists_spa = new ANNdist[ k_spa ];
      m_kdTreeSpace -> annkSearch(spatialPt, k_spa, nnIdx_spa, dists_spa, 0);
      r0 = std::sqrt(dists_spa[k_spa-1])/4.0;
    } else
    {
      nnIdx_spa = new ANNidx[ k_spa ];
      dists_spa = new ANNdist[ k_spa ];
      m_kdTreeSpace -> annkFRSearch(spatialPt, 16*m_Rspa*m_Rspa, k_spa, nnIdx_spa, dists_spa, 0);
    }

    // Fill data structures for RBF interpolation

    MatDoub pts(k_spa, ImageDimension);
    VecDoub y(k_spa);

    for( unsigned int i=0; i< k_spa; ++i)
      {
        y[i] = m_dataVals[ nnIdx_spa[i] ];
        pts[i][0] = m_dataPtsSpace[ nnIdx_spa[i] ][0];
        pts[i][1] = m_dataPtsSpace[ nnIdx_spa[i] ][1];
        pts[i][2] = m_dataPtsSpace[ nnIdx_spa[i] ][2];

//        if ( (index[0]==66) && (index[1]==39) && (index[2]==25))
//        {
//          std::cout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << ";" << std::endl;
//        }

      }

    // Cleaning

    delete [] nnIdx_spa;
    delete [] dists_spa;

    double max_val = y[0];
    double min_val = y[0];

    for( unsigned int i=1; i< k_spa; ++i)
    {
      if (y[i] > max_val) max_val = y[i];
      if (y[i] < min_val) min_val = y[i];
    }


    RBF_gauss gaussian(r0);
    RBF_interp RBFInterpolator(pts,y,gaussian,1);

    RealType value = NumericTraits<RealType>::Zero;

    VecDoub pt(3);
    pt[0] = point[0]; pt[1] = point[1]; pt[2] = point[2];

    value = RBFInterpolator.interp(pt);

    if ( (value<min_val) || (value>max_val) )
    {
      value = y[0];
    }


    return ( static_cast<OutputType>( value ) );

}


} // end namespace btk

#endif

