#ifndef _btkLowToHighImageResolutionMethod_txx
#define _btkLowToHighImageResolutionMethod_txx

#include "btkLowToHighImageResolutionMethod.h"

namespace btk
{

/*
 * Constructor
 */
template < typename ImageType >
LowToHighImageResolutionMethod<ImageType>
::LowToHighImageResolutionMethod()
{

  m_TransformArray.resize(0);
  m_Interpolator = 0;
  m_NumberOfImages = 0;
  m_TargetImage = 0;
  m_InitializeWithMask= false;

}

template < typename ImageType >
void LowToHighImageResolutionMethod<ImageType>
::SetNumberOfImages(int N)
{

  m_TransformArray.resize(N);
  m_ImageArray.resize(N);
  m_ResampledImageArray.resize(N);
  m_RegionArray.resize(N);
  m_ImageMaskArray.resize(N);

  for(int i=m_NumberOfImages; i<N; i++)
  {
    m_ImageArray[i]=0;
    m_ResampledImageArray[i]=0;
  }

  m_NumberOfImages  = N;

}
/*
 * Initialize by setting the interconnects between components.
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::Initialize() throw (ExceptionObject)
{

/*  m_Registration = RegistrationType::New();
  m_Registration -> SetFixedImage(  m_ImageArray[m_TargetImage]  );
  m_Registration -> SetFixedImageRegion( m_RegionArray[m_TargetImage] );

  // FIXME The following lines should be executed only if masks are set

  m_Registration -> SetFixedImageMask( m_ImageMaskArray[m_TargetImage] );
  m_Registration -> InitializeWithMask();
  m_Registration -> SetEnableObserver( false ); */

  // Initial rigid parameters computed from the center of ROIs
/*  PointType centerPointTarget = m_Registration -> GetRotationCenter();

  m_InitialRigidParameters.resize(m_NumberOfImages);

  for (unsigned int i=0; i<m_NumberOfImages; i++)
  {
    // if the user doesn't provide masks, we initialize the transformation with
    // rois, if not we use the masks
    if (m_InitializeWithMask)
    {
      m_InitialRigidParameters[i].SetSize(6);
      m_InitialRigidParameters[i].Fill(0.0);

      IndexType regionStart = m_RegionArray[i].GetIndex();
      SizeType  regionSize  = m_RegionArray[i].GetSize();

      IndexType centerIndex;
      centerIndex[0] = regionStart[0] + regionSize[0]/2;
      centerIndex[1] = regionStart[1] + regionSize[1]/2;
      centerIndex[2] = regionStart[2] + regionSize[2]/2;

      PointType centerPoint;
      m_ImageArray[i] -> TransformIndexToPhysicalPoint(centerIndex,centerPoint);

      m_InitialRigidParameters[i][3] = centerPoint[0] - centerPointTarget[0];
      m_InitialRigidParameters[i][4] = centerPoint[1] - centerPointTarget[1];
      m_InitialRigidParameters[i][5] = centerPoint[2] - centerPointTarget[2];
    } else
    {

      TransformInitializerType::Pointer initializer = TransformInitializerType::New();
      EulerTransformType::Pointer transform = EulerTransformType::New();

      initializer->SetTransform( transform );
      initializer->SetFixedImage(  m_ImageMaskArray[m_TargetImage] );
      initializer->SetMovingImage( m_ImageMaskArray[i] );
      initializer->MomentsOn();
      initializer->InitializeTransform();
      m_InitialRigidParameters[i] = transform -> GetParameters();

    }
  }

*/


   /* Create interpolator */
   m_Interpolator = InterpolatorType::New();

  /* Resampling matrix */
  SpacingType  fixedSpacing   = m_ImageArray[m_TargetImage] -> GetSpacing();
  SizeType     fixedSize      = m_ImageArray[m_TargetImage] -> GetLargestPossibleRegion().GetSize();

  m_ResampleSpacing[0] = fixedSpacing[0];
  m_ResampleSpacing[1] = fixedSpacing[0];
  m_ResampleSpacing[2] = fixedSpacing[0];

  m_ResampleSize[0] = floor(fixedSize[0]*fixedSpacing[0]/m_ResampleSpacing[0]);
  m_ResampleSize[1] = floor(fixedSize[1]*fixedSpacing[1]/m_ResampleSpacing[1]);
  m_ResampleSize[2] = floor(fixedSize[2]*fixedSpacing[2]/m_ResampleSpacing[2]);

  /* Create high resolution image */

  m_HighResolutionImage = ImageType::New();

  IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  RegionType region;
  region.SetIndex(start);
  region.SetSize(m_ResampleSize);

  m_HighResolutionImage -> SetRegions( region );
  m_HighResolutionImage -> Allocate();
  m_HighResolutionImage -> SetOrigin( m_ImageArray[m_TargetImage] -> GetOrigin() );
  m_HighResolutionImage -> SetSpacing( m_ResampleSpacing );
  m_HighResolutionImage -> SetDirection( m_ImageArray[m_TargetImage]-> GetDirection() );
  m_HighResolutionImage -> FillBuffer( 0 );

}

/*
 * Starts the Registration Process
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::StartRegistration( void )
{

  try
  {
    // initialize the interconnects between components
    this->Initialize();
  }
  catch( ExceptionObject& err )
  {
    // pass exception to caller
    throw err;
  }

  RegistrationPointer m_Registration = RegistrationType::New();
  m_Registration -> SetFixedImage(  m_ImageArray[m_TargetImage]  );
  m_Registration -> SetFixedImageRegion( m_RegionArray[m_TargetImage] );

  // FIXME The following lines should be executed only if masks are set

  m_Registration -> SetFixedImageMask( m_ImageMaskArray[m_TargetImage] );
  m_Registration -> InitializeWithMask();
  m_Registration -> SetEnableObserver( false );

  m_ResamplingStatus.resize( m_NumberOfImages );
  for (unsigned int i=0; i < m_ResamplingStatus.size(); i++)
    m_ResamplingStatus[i] = false;

  for (unsigned int i=0; i < m_NumberOfImages; i++)
  {
    m_Registration->SetMovingImage( m_ImageArray[i] );
    m_Registration->SetMovingImageMask( m_ImageMaskArray[i] );
//      m_Registration->SetInitialTransformParameters( m_InitialRigidParameters[i] );

    if (i != m_TargetImage )
    {
      std::cout << "Registering image " << i << " to " << m_TargetImage << " ... "; std::cout.flush();

      try
      {
        m_Registration->StartRegistration();
      }
      catch( itk::ExceptionObject & err )
      {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        throw err;
      }

      std::cout << "done." << std::endl; std::cout.flush();
      m_TransformArray[i] = m_Registration->GetTransform();
      }
      else
      {
        m_TransformArray[i] = TransformType::New();
        m_TransformArray[i] -> SetIdentity();
//        m_TransformArray[i] -> SetCenter( m_Registration -> GetTransformCenter() );

      }

    }

    IteratorType imageIt( m_HighResolutionImage, m_HighResolutionImage->GetLargestPossibleRegion() );

    PointType physicalPoint;
    PointType transformedPoint;

    IndexType index;
    int value;
    unsigned int counter;

//    std::cout << "Creating initial HR image ... " ; std::cout.flush();

    for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
    {
      index = imageIt.GetIndex();
      m_HighResolutionImage -> TransformIndexToPhysicalPoint( index, physicalPoint);

      value = 0;
      counter = 0;
      for (unsigned int i=0; i < m_NumberOfImages; i++)
      {
        m_Interpolator -> SetInputImage( m_ImageArray[i] );
        transformedPoint = m_TransformArray[i]->TransformPoint(physicalPoint);
        if ( m_Interpolator -> IsInsideBuffer (transformedPoint) )
        {
          value+= m_Interpolator->Evaluate(transformedPoint);
          counter++;
        }
      }
      if ( counter>0 ) imageIt.Set(value/counter);
     }
//     std::cout << "done." << std::endl; std::cout.flush();

}

/*
 * Writes transforms to file
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::WriteTransforms( const char *folder )
{

  typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();

  for (unsigned int i=0; i < m_NumberOfImages; i++)
  {

    transformWriter->SetInput( m_TransformArray[i] );

    char fullTrafoName[255]; strcpy ( fullTrafoName, folder );
    char trafoName[255];

    sprintf ( trafoName, "/%d.txt", i );
    strcat ( fullTrafoName,trafoName );

    transformWriter -> SetFileName ( fullTrafoName );

    try
    {
      std::cout << "Writing transform to " << fullTrafoName << " ... " ; std::cout.flush();
      transformWriter -> Update();
      std::cout << " done! " << std::endl;
    }
    catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while saving the transforms" << std::endl;
      std::cerr << excp << std::endl;
      std::cout << "[FAILED]" << std::endl;
      throw excp;
    }

  }

}

/*
 * Writes a specific transform to file
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::WriteTransforms( unsigned int i, const char *filename )
{

  typename TransformWriterType::Pointer transformWriter = TransformWriterType::New();
  transformWriter->SetInput( m_TransformArray[i] );
  transformWriter -> SetFileName ( filename );

  try
  {
    std::cout << "Writing transform to " << filename << " ... " ; std::cout.flush();
    transformWriter -> Update();
    std::cout << " done! " << std::endl;
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while saving the transforms" << std::endl;
    std::cerr << excp << std::endl;
    std::cout << "[FAILED]" << std::endl;
    throw excp;
  }

}

/*
 * Writes resampled images to disk
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::WriteResampledImages( const char *folder )
{

  typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();

  for (unsigned int i=0; i < m_NumberOfImages; i++)
  {
    if (!m_ResamplingStatus[i])
    {
      m_Resample =  ResampleType::New();
      m_Resample -> SetTransform( m_TransformArray[i] );
      m_Resample -> SetInput( m_ImageArray[i] );
      m_Resample -> SetReferenceImage( m_HighResolutionImage );
      m_Resample -> SetUseReferenceImage( true );
      m_Resample -> SetDefaultPixelValue( 0 );
      m_Resample -> Update();
      m_ResampledImageArray[i] = m_Resample -> GetOutput();
      m_ResamplingStatus[i] = true;
    }

    imageWriter->SetInput( m_ResampledImageArray[i] );

    char fullImageName[255]; strcpy ( fullImageName, folder );
    char imageName[255];

    sprintf ( imageName, "/%d2HR.nii.gz", i );
    strcat ( fullImageName, imageName );

    imageWriter -> SetFileName ( fullImageName );

    try
    {
      std::cout << "Writing resampled image to " << fullImageName << " ... " ; std::cout.flush();
      imageWriter -> Update();
      std::cout << " done! " << std::endl;
    }
    catch ( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while saving the resampled image" << std::endl;
      std::cerr << excp << std::endl;
      std::cout << "[FAILED]" << std::endl;
      throw excp;
    }

  }

}

/*
 * Writes a specific resampled image to disk
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::WriteResampledImages( unsigned int i, const char *filename )
{

  if (!m_ResamplingStatus[i])
  {
    m_Resample =  ResampleType::New();
    m_Resample -> SetTransform( m_TransformArray[i] );
    m_Resample -> SetInput( m_ImageArray[i] );
    m_Resample -> SetReferenceImage( m_HighResolutionImage );
    m_Resample -> SetUseReferenceImage( true );
    m_Resample -> SetDefaultPixelValue( 0 );
    m_Resample -> Update();
    m_ResampledImageArray[i] = m_Resample -> GetOutput();
    m_ResamplingStatus[i] = true;
  }

  typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();

  imageWriter -> SetInput( m_ResampledImageArray[i] );
  imageWriter -> SetFileName ( filename );

  try
  {
    std::cout << "Writing resampled image to " << filename << " ... " ; std::cout.flush();
    imageWriter -> Update();
    std::cout << " done! " << std::endl;
  }
  catch ( itk::ExceptionObject & excp )
  {
    std::cerr << "Error while saving the resampled image" << std::endl;
    std::cerr << excp << std::endl;
    std::cout << "[FAILED]" << std::endl;
    throw excp;
  }

}


/*
 * PrintSelf
 */
template < typename ImageType >
void
LowToHighImageResolutionMethod<ImageType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


} // end namespace btk


#endif
