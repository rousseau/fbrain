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

#ifndef _btkSliceBySliceRegistration_txx
#define _btkSliceBySliceRegistration_txx

#include "btkSliceBySliceRegistration.h"

namespace btk
{

/*
 * Constructor
 */
template < typename ImageType >
SliceBySliceRegistration<ImageType>
::SliceBySliceRegistration()
{
  m_TransformArrayIsSet = false;
}

/*
 * Initialize by setting the interconnects between components.
 */
template < typename ImageType >
void
SliceBySliceRegistration<ImageType>
::Initialize() throw (ExceptionObject)
{

  Superclass::Initialize();

}

/*
 * PrintSelf
 */
template < typename ImageType >
void
SliceBySliceRegistration<ImageType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}

template < typename ImageType >
void
SliceBySliceRegistration<ImageType>
::StartRegistration()
{

// Initialize transformations if are not set

  if (!m_TransformArrayIsSet)
  {
    SizeType fullSize = this -> GetFixedImage() -> GetLargestPossibleRegion().GetSize();
    m_TransformArray.resize( fullSize[2] );

    for (unsigned int k=0; k<fullSize[2]; k++)
    {
      m_TransformArray[k] =  TransformType::New();
    }
  }

//  Registration objects

  rigidRegistration   = RigidRegistrationType::New();
  affineRegistration  = AffineRegistrationType::New();

  rigidRegistration->SetMovingImage(  this -> GetMovingImage() );
  rigidRegistration->SetFixedImage(   this -> GetFixedImage()  );

  affineRegistration->SetMovingImage( this -> GetMovingImage() );
  affineRegistration->SetFixedImage(  this -> GetFixedImage() );

  RigidTransformType::Pointer transform = RigidTransformType::New();

  IndexType start = this -> GetFixedImageRegion().GetIndex();
  SizeType  size  = this -> GetFixedImageRegion().GetSize();

  unsigned int x1 = start[0];
  unsigned int y1 = start[1];
  unsigned int z1 = start[2];

  unsigned int x2 = x1 + size[0] -1;
  unsigned int y2 = y1 + size[1] -1;
  unsigned int z2 = z1 + size[2] -1;

  for ( unsigned int i=z1; i<=z2; i++ )
  {
    // Fixed region for slice i

    RegionType fixedImageRegion;
    IndexType  fixedImageRegionIndex;
    SizeType   fixedImageRegionSize;

    fixedImageRegionIndex[0] = x1; fixedImageRegionIndex[1]= y1; fixedImageRegionIndex[2]= i;
    fixedImageRegionSize[0]  = x2 - x1 + 1; fixedImageRegionSize[1] = y2 - y1 + 1; fixedImageRegionSize[2] = 1;

    fixedImageRegion.SetIndex(fixedImageRegionIndex);
    fixedImageRegion.SetSize(fixedImageRegionSize);

    ParametersType initialAffineParameters( 12 );

    if (!m_TransformArrayIsSet)
    {

      rigidRegistration->SetFixedImageRegion( fixedImageRegion );

      try
      {
        rigidRegistration -> StartRegistration();
      }
      catch( itk::ExceptionObject & err )
      {
        throw err;
      }

      ParametersType finalParameters = rigidRegistration -> GetLastTransformParameters();
      PointType      transformCenter = rigidRegistration -> GetTransformCenter();

      // Initialize affine
      transform -> SetCenter( transformCenter);
      transform -> SetParameters( finalParameters );
      VnlMatrixType rigidMatrix = transform -> GetMatrix().GetVnlMatrix();

      unsigned short p = 0;
      for ( int ii=0; ii<3; ii++)
      {
        for ( int jj = 0; jj < 3; jj++)
        {
          initialAffineParameters[p] = rigidMatrix(ii,jj);
          p++;
        }
      }
      initialAffineParameters[9] = finalParameters[3];
      initialAffineParameters[10] = finalParameters[4];
      initialAffineParameters[11] = finalParameters[5];

    }
    else
    {
      initialAffineParameters = m_TransformArray[i] -> GetParameters();

    }

    affineRegistration->SetFixedImageRegion( fixedImageRegion );
    affineRegistration->SetInitialTransformParameters( initialAffineParameters );

//    std::cout << "Initial affine parameters = " << initialAffineParameters << std::endl;

    try
      {
      affineRegistration -> StartRegistration();
      }
    catch( itk::ExceptionObject & err )
      {
      throw err;
      }

//    finalParameters = affineRegistration->GetLastTransformParameters();

//    std::cout << "Final affine parameters ( z = " << i << " ) " << affineRegistration->GetLastTransformParameters() << std::endl;

    m_TransformArray[i] -> SetIdentity();
    m_TransformArray[i] -> SetCenter( affineRegistration->GetTransformCenter() );
    m_TransformArray[i] -> SetParameters( affineRegistration->GetLastTransformParameters() );

  } // end for in z

}


} // end namespace btk


#endif
