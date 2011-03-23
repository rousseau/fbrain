/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/03/2010
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

#ifndef __btkRBFInterpolateImageFunctionS2S_txx
#define __btkRBFInterpolateImageFunctionS2S_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

#include "btkRBFInterpolateImageFunctionS2S.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_inverse.h"


namespace btk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep>
const unsigned long
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep>
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::RBFInterpolateImageFunctionS2S()
{
  m_GradientTableCartesian = 0;
  m_kdTreeSphere = 0;
  m_TransformsAreSet = false;
  m_NumberOfSlices = 0;
  m_FirstSlice = 0;
  m_LastSlice  = 0;
}

/**
 * Destructor
 */
template<class TInputImage, class TCoordRep>
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::~RBFInterpolateImageFunctionS2S()
{
  for (unsigned int i=0; i<m_NumberOfGradients*m_NumberOfSlices; i++)
  {
    delete m_kdTreeSpace[i];
  }

  if (m_kdTreeSphere != 0)
  {
    delete m_kdTreeSphere;
  }

  annClose();

}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}

template<class TInputImage, class TCoordRep>
void
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::SetTransforms( const char* tpath )
{
  for (unsigned int i=1; i<m_ImageSize[3]; i++)
  {

    char fullTrafoName[255]; strcpy ( fullTrafoName,tpath );
    char trafoName[255];

    sprintf ( trafoName, "/%d.txt", i );
    strcat ( fullTrafoName,trafoName );

    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader->SetFileName( fullTrafoName );
    transformReader->Update();

    TransformListType transforms = transformReader->GetTransformList();

    TransformReaderType::TransformListType::const_iterator titr =
    transforms->begin();

//    std::cout << "transform size in " << fullTrafoName << " " << transforms -> size() << std::endl;

    for(unsigned int j=0; j<transforms -> size(); j++,titr++)
    {
      if( !strcmp((*titr)->GetNameOfClass(),"AffineTransform") )
        {
        m_TransformArray[i-1][j] = dynamic_cast< TransformType * >( titr->GetPointer() );

        if( !m_TransformArray[i-1][j] )//  m_kdTreeSpace = 0;
          {
          std::cerr << "Error reading Transform" << std::endl;
          }

        }
      else
        {
        std::cerr << "Input file does not contain an Affinr Transform" << std::endl;
        }
    }

  }

  m_TransformsAreSet = true;

}


template<class TInputImage, class TCoordRep>
void
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
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

/*  std::cout << "Original gradient table = " << std::endl;
  for (unsigned int r=0; r < m_ImageSize[3]; r++)
  {
    for (unsigned int c=0; c < 3; c++)
    {
      std::cout << m_OriginalGradientTable(r,c) << " ";
    }
    std::cout << std::endl;

  } */

}

template<class TInputImage, class TCoordRep>
void
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::RotateGradients()
{

  // data structure for nearest neighbor search with ANN
  m_dataPtsSphere = annAllocPts(2*m_NumberOfGradients*m_NumberOfSlices, 3);

  double x,y,z,r;

  unsigned int k=0;

  for (unsigned int i=1; i <= m_NumberOfGradients; i++)
  {
    for (unsigned int j=m_FirstSlice; j<=m_LastSlice; j++)
    {

      VnlMatrixType Md = m_TransformArray[i-1][j] -> GetMatrix().GetVnlMatrix();

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
           break;
        }
        else
        {
          PQd = NQd;
        }
      }

      VnlVectorType grad = m_GradientTableCartesian.get_row(i);
      VnlVectorType rotGrad = NQd*grad;

      x = rotGrad(0);
      y = rotGrad(1);
      z = rotGrad(2);

      r = std::sqrt(x*x + y*y + z*z);

      // re-normalization to avoid errors due to accuracy in gradients
      x = x/r; y = y/r; z = z/r; r=1;

      // For the nearest neighbor search, we use the negative point too
      // since they represent the same direction

      m_dataPtsSphere[k][0] = x;
      m_dataPtsSphere[k][1] = y;
      m_dataPtsSphere[k][2] = z;

      m_dataPtsSphere[ k + m_NumberOfGradients*m_NumberOfSlices ][0] = -x;
      m_dataPtsSphere[ k + m_NumberOfGradients*m_NumberOfSlices ][1] = -y;
      m_dataPtsSphere[ k + m_NumberOfGradients*m_NumberOfSlices ][2] = -z;

      k++;

    }

  }

/*  std::cout << "[";
  for (k=0; k<2*m_NumberOfGradients*m_NumberOfSlices; k++)
  {
    std::cout << m_dataPtsSphere[k][0] << " " << m_dataPtsSphere[k][1] << " " << m_dataPtsSphere[k][2] << ";" << std::endl;
  }
  std::cout << "]" << std::endl; */

  m_kdTreeSphere = new ANNkd_tree(m_dataPtsSphere, 2*m_NumberOfGradients*m_NumberOfSlices, 3);

}

template<class TInputImage, class TCoordRep>
void
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
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
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::SetInputImage(const InputImageType *ptr)
{
  this->Superclass::SetInputImage(ptr);
  m_ImageSize = this->GetInputImage()->GetLargestPossibleRegion().GetSize();
  m_TransformArray.resize(m_ImageSize[3]-1);

  for(unsigned j=0; j<m_ImageSize[3]-1; j++)
  {
    m_TransformArray[j].resize(m_ImageSize[2]);
  }

//  std::cout << "tamanios = " << m_TransformArray.size() << " " << m_TransformArray[0].size() << std::endl;

//  m_InverseTransformArray.resize(m_ImageSize[3]-1);
  m_NumberOfGradients = m_ImageSize[3]-1;
}

template<class TInputImage, class TCoordRep>
void
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::Initialize( ImageRegionType & region)
{
  // In case of leave1out

  if (! m_TransformsAreSet )
  {
    for (unsigned int i=0; i < m_TransformArray.size(); i++)
    {
      for (unsigned int j=0; j< m_TransformArray[i].size(); j++)
      {
        m_TransformArray[i][j] = TransformType::New();
        m_TransformArray[i][j]-> SetIdentity();

      }
    }
  }

  ImageRegionType region2 = region;

  ImageSizeType size2 = region2.GetSize();
  ImageIndexType start2 = region2.GetIndex();

  size2[3] = 1;
  start2[3] = 1;

  region2.SetSize(size2);
  region2.SetIndex(start2);

  m_NumberOfSlices = size2[2];
  m_FirstSlice     = start2[2];
  m_LastSlice      = start2[2] + size2[2] - 1;

  this -> RotateGradients();

  std::cout << m_NumberOfSlices << " " << m_FirstSlice << " " << m_LastSlice << std::endl;

  IndexType index;
  PointType point;
  typename TransformType::OutputPointType transformedPoint;
  typename TransformType::InputPointType spatialPoint;

  m_dataPtsSpace.resize(m_NumberOfGradients*m_NumberOfSlices);
  m_dataVals.resize(m_NumberOfGradients*m_NumberOfSlices);

  m_kdTreeSpace.resize(m_NumberOfGradients*m_NumberOfSlices);

  unsigned int nn = 0;

  for (unsigned int n=0; n<m_NumberOfGradients; n++ )
  {
    for (unsigned int j=m_FirstSlice; j<=m_LastSlice; j++)
    {

      m_dataPtsSpace[nn] = annAllocPts( size2[0]*size2[1], ImageDimension-1);
      m_dataVals[nn] = annAllocPts( size2[0]*size2[1], 1);

      unsigned int i = 0;

      start2[2] = j;
      size2[2] = 1;

      region2.SetIndex( start2 );
      region2.SetSize( size2 );

      IteratorType regionIt(this->GetInputImage(),region2);

      for (regionIt.GoToBegin(); !regionIt.IsAtEnd(); ++regionIt)
      {

        index = regionIt.GetIndex();
        this->GetInputImage() -> TransformIndexToPhysicalPoint(index,point);

        spatialPoint[0] = point[0];
        spatialPoint[1] = point[1];
        spatialPoint[2] = point[2];

        transformedPoint = m_TransformArray[n][j]->TransformPoint(spatialPoint);

        m_dataPtsSpace[nn][i][0] = transformedPoint[0];
        m_dataPtsSpace[nn][i][1] = transformedPoint[1];
        m_dataPtsSpace[nn][i][2] = transformedPoint[2];

        m_dataVals[nn][i][0] = regionIt.Get();

        i++;

      }

      m_kdTreeSpace[nn] = new ANNkd_tree(m_dataPtsSpace[nn], size2[0]*size2[1], ImageDimension-1);

      nn++;

    }

    start2[3] = start2[3]+1;


  }

}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::OutputType
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
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
typename RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::OutputType
RBFInterpolateImageFunctionS2S< TInputImage, TCoordRep >
::Evaluate( const PointType& point,
            double theta, double phi,
            double r_spa, double r_gra, char init) const
{

    unsigned int k_gra;
    std::vector< unsigned int > k_spa;

    // Search k_spa spatial neighbors

    ANNpoint spatialPt;
    spatialPt = annAllocPt(3);
    spatialPt[0] = point[0];
    spatialPt[1] = point[1];
    spatialPt[2] = point[2];

    // Search k_gra neighbors in the sphere closer than 3*r_grad

    ANNpoint spherePt;
    spherePt = annAllocPt(3);
    spherePt[0] = std::sin(theta) * std::cos(phi);//m_GradientTableCartesian[1][0];
    spherePt[1] = std::sin(theta) * std::sin(phi);
    spherePt[2] = std::cos(theta);

    k_gra = m_kdTreeSphere -> annkFRSearch(spherePt, 9*r_gra*r_gra, 0);

    ANNidxArray nnIdx_gra;
    nnIdx_gra = new ANNidx[k_gra];

    ANNdistArray dists_gra;
    dists_gra = new ANNdist[k_gra];

    m_kdTreeSphere -> annkFRSearch(spherePt, 9*r_gra*r_gra, k_gra, nnIdx_gra, dists_gra, 0);

    // Calculate number of spatial neighbors

    k_spa.resize(k_gra);
    unsigned int numberOfPoints = 0;
    bool fixedNN = false;

//    std::cout << "Proximos al punto " << spherePt[0] << " " << spherePt[1] << " " << spherePt[2] << " : " << std::endl;

    for(unsigned int j=init; j < k_gra; j++)
    {
      unsigned int gradIndex =  nnIdx_gra[j] % (m_NumberOfGradients*m_NumberOfSlices);

//      std::cout << "gradIndex = " << gradIndex << " " <<  m_dataPtsSphere[ nnIdx_gra[j] ][0] <<
//      " " << m_dataPtsSphere[ nnIdx_gra[j] ][1] << " " << m_dataPtsSphere[ nnIdx_gra[j] ][2] << std::endl;

      k_spa[j] = m_kdTreeSpace[gradIndex] -> annkFRSearch(spatialPt, 9*r_spa*r_spa, 0);
      numberOfPoints += (k_spa[j] - init);
    }

    if (numberOfPoints == 0)
    {
      std::cout << "Warning: no neighbors found in the search area. Forcing search." << std::endl;
      for(unsigned int j=init; j < k_gra; j++)
      {
        unsigned int gradIndex =  nnIdx_gra[j] % (m_NumberOfGradients*m_NumberOfSlices);

        k_spa[j] = 4;
        numberOfPoints += (k_spa[j] - init);
      }
      fixedNN = true;
    }

    MatDoub pts(numberOfPoints,5);
    VecDoub y(numberOfPoints);

    typename IteratorType::IndexType index;
    typename InputImageType::PointType physicalPoint;

    typename TransformType::InputPointType spatialPoint;
    typename TransformType::OutputPointType transformedPoint;

    unsigned int n = 0;
    double r_g, x_g, y_g, z_g;

    for(unsigned int j=init; j < k_gra; j++)
    {

      if ( k_spa[j] == 0)
      {
        continue;
      }

      x_g = m_dataPtsSphere[ nnIdx_gra[j] ][0];
      y_g = m_dataPtsSphere[ nnIdx_gra[j] ][1];
      z_g = m_dataPtsSphere[ nnIdx_gra[j] ][2];
      r_g = std::sqrt(x_g*x_g + y_g*y_g + z_g*z_g);

      unsigned int gradIndex =  nnIdx_gra[j] % (m_NumberOfGradients*m_NumberOfSlices);
//      index[0] = 0; index[1] = 0; index[2] = 0; index[3] = gradIndex+1;
//      this->GetInputImage()->TransformIndexToPhysicalPoint( index, physicalPoint );

      ANNidxArray nnIdx_spa;
      nnIdx_spa = new ANNidx[ k_spa[j] ];

      ANNdistArray dists_spa;
      dists_spa = new ANNdist[ k_spa[j] ];

      if (fixedNN)
      {
        m_kdTreeSpace[gradIndex] -> annkSearch(spatialPt, k_spa[j], nnIdx_spa, dists_spa, 0);
      } else
       {
         m_kdTreeSpace[gradIndex] -> annkFRSearch(spatialPt, 9*r_spa*r_spa, k_spa[j], nnIdx_spa, dists_spa, 0);
       }

      for( unsigned int i=init; i< k_spa[j]; ++i)
      {

        spatialPoint[0] = m_dataPtsSpace[gradIndex][ nnIdx_spa[i] ][0];
        spatialPoint[1] = m_dataPtsSpace[gradIndex][ nnIdx_spa[i] ][1];
        spatialPoint[2] = m_dataPtsSpace[gradIndex][ nnIdx_spa[i] ][2];

//        transformedPoint = m_InverseTransformArray[gradIndex]->TransformPoint(spatialPoint);
//
//        physicalPoint[0] = transformedPoint[0];
//        physicalPoint[1] = transformedPoint[1];
//        physicalPoint[2] = transformedPoint[2];

//        this -> GetInputImage() -> TransformPhysicalPointToIndex( physicalPoint,index );
//        y[n] = this -> GetInputImage() -> GetPixel( index );

        y[n] = m_dataVals[gradIndex][nnIdx_spa[i]][0];

        // First point components are the voxel coordinates
        pts[n][0] = spatialPoint[0];
        pts[n][1] = spatialPoint[1];
        pts[n][2] = spatialPoint[2];

        // Last point components are the \phi and \theta gradient angles

        pts[n][3] = std::acos(z_g/r_g);
        pts[n][4] = std::atan2(y_g,x_g);

        n++;
      }

      delete [] nnIdx_spa;
      delete [] dists_spa;

    }

    double max_val = y[0];
    double min_val = y[0];

    for( unsigned int i=1; i< numberOfPoints; ++i)
    {
      if (y[i] > max_val) max_val = y[i];
      if (y[i] < min_val) min_val = y[i];
    }



    RBF2_gauss gaussian(r_spa,r_gra);
    RBF_interp RBFInterpolator(pts,y,gaussian,1);

    RealType value = NumericTraits<RealType>::Zero;

    VecDoub pt(5);
    pt[0] = point[0]; pt[1] = point[1]; pt[2] = point[2];
    pt[3] = theta; pt[4] = phi;

    value = RBFInterpolator.interp(pt);


    if ( (value<min_val) || (value>max_val) )
    {
      value = y[0];
    }

//    if ( fixedNN )
//    {
//      value = 200;
//    }

    // cleaning

    delete [] nnIdx_gra;
    delete [] dists_gra;

  return ( static_cast<OutputType>( value ) );

}


} // end namespace btk

#endif

