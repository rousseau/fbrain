#ifndef _btkSliceBySliceRigidRegistration_txx
#define _btkSliceBySliceRigidRegistration_txx

#include "btkSliceBySliceRigidRegistration.h"

namespace btk
{

/*
 * Constructor
 */
template < typename ImageType >
SliceBySliceRigidRegistration<ImageType>
::SliceBySliceRigidRegistration()
{
  m_Interpolator = 0;
  m_ImageMask = 0;
  m_Transform = 0;
  m_Iterations = 200;
}

/*
 * Initialize by setting the interconnects between components.
 */
template < typename ImageType >
void
SliceBySliceRigidRegistration<ImageType>
::Initialize() throw (itk::ExceptionObject)
{
  // Configure registration

  m_Registration = RegistrationType::New();
  m_Registration -> SetFixedImage(  m_FixedImage  );
  m_Registration -> SetMovingImage( m_MovingImage );
  m_Registration -> InitializeWithTransform();
  m_Registration -> SetEnableObserver( false );
  m_Registration -> SetIterations( m_Iterations );

  // TODO We have to decide after checking the results which one is the
  // the default behavior

  if ( !m_ImageMask )
  {
//    std::cout << "image mask IS NOT defined" << std::endl;
    m_ROI  = m_FixedImage -> GetLargestPossibleRegion().GetSize();
  } else
    {
//      std::cout << "image mask IS defined" << std::endl;
      typename MaskType::Pointer mask = MaskType::New();
      mask -> SetImage( m_ImageMask );
      m_Registration -> SetFixedImageMask( m_ImageMask );
      m_ROI = mask -> GetAxisAlignedBoundingBoxRegion();
    }
    
//  if ( !m_Transform)
//  {
//    m_Transform = SliceBySliceTransformType::New();
//    m_Transform -> SetImage( m_FixedImage );
//    m_Transform -> Initialize();
//  }

}

/*
 * Starts the Registration Process
 */
template < typename ImageType >
void
SliceBySliceRigidRegistration<ImageType>
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

    unsigned int k1 =  m_ROI.GetIndex()[2];
    unsigned int k2 =  k1 + m_ROI.GetSize()[2] -1;

    for ( unsigned int i = k1; i <= k2; i++ )
    {
//      std::cout << "Registering slice " << i << std::endl;

      // Fixed region for slice i

      RegionType fixedImageRegion;
      IndexType  fixedImageRegionIndex;
      SizeType   fixedImageRegionSize;

      fixedImageRegionIndex = m_ROI.GetIndex();
      fixedImageRegionIndex[2] = i;

      fixedImageRegionSize = m_ROI.GetSize();
      fixedImageRegionSize[2] = 1;

      fixedImageRegion.SetIndex(fixedImageRegionIndex);
      fixedImageRegion.SetSize(fixedImageRegionSize);

      ParametersType initialRigidParameters( 6 );
      ParametersType finalParameters;

      m_Registration -> SetFixedImageRegion( fixedImageRegion );
      m_Registration -> SetInitialTransformParameters( m_Transform -> GetSliceTransform(i) -> GetParameters() );
      m_Registration -> SetTransformCenter( m_Transform -> GetSliceTransform(i) -> GetCenter() );

//      std::cout << "Initial registration parameters = " << m_Transform -> GetSliceTransform(i) -> GetParameters() << std::endl;

      try
        {
        //m_Registration -> StartRegistration();// FIXME : in ITK4 StartRegistration() is replaced by Update()
        m_Registration->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        // TODO: Always check in case of problems: the following lines have been commented
        // on purpose (if the registration fails we prefer keep the old parameters)

  //      std::cerr << "ExceptionObject caught !" << std::endl;
  //      std::cerr << err << std::endl;
  //      return EXIT_FAILURE;
        }

      finalParameters = m_Registration -> GetLastTransformParameters();

//      std::cout << "Final rigid parameters = " << finalParameters << std::endl;

      m_Transform -> SetSliceParameters( i, finalParameters );

    } // end for in z

}

/*
 * PrintSelf
 */
template < typename ImageType >
void
SliceBySliceRigidRegistration<ImageType>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


} // end namespace btk


#endif
