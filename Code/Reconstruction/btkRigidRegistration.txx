/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 12/11/2010
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

#ifndef _btkRigidRegistration_cxx
#define _btkRigidRegistration_cxx

#include "btkRigidRegistration.h"

namespace btk
{

/*
 * Constructor
 */
template < typename ImageType >
RigidRegistration<ImageType>
::RigidRegistration()
{
  m_Iterations = 200;
  m_EnableObserver = false;
  m_FixedImageMask = 0;
  m_MovingImageMask = 0;

  m_InitializeWithTransform = false;
  m_InitializeWithMask = false;
  m_InitializeWithRegion = false;
}

/*
 * Initialize by setting the interconnects between components.
 */
template < typename ImageType >
void
RigidRegistration<ImageType>
::Initialize() throw (ExceptionObject)
{
  // Configure transform

  m_Transform = TransformType::New();

  if (m_InitializeWithMask)
  {
    TransformInitializerType::Pointer initializer = TransformInitializerType::New();

    initializer -> SetTransform( m_Transform );

    // TODO Perhaps here I could add methods to set masks (as spatial objects)
    // and then recover images by using the GetImage() method.

    initializer -> SetFixedImage(  this -> GetFixedImageMask() );
    initializer -> SetMovingImage( this -> GetMovingImageMask() );
    initializer -> MomentsOn();
    initializer -> InitializeTransform();
  } else if (m_InitializeWithTransform)
    {
      m_Transform -> SetParameters( this -> GetInitialTransformParameters() );
      m_Transform -> SetCenter( this -> GetTransformCenter() );
    }

  this -> SetInitialTransformParameters( m_Transform -> GetParameters() );

//  std::cout << "Transform center = " << m_Transform -> GetCenter() << std::endl;
//  std::cout << "Transform parameters = " << m_Transform -> GetParameters() << std::endl;

  m_Interpolator = InterpolatorType::New();
  m_Metric = MetricType::New();
  m_Optimizer = OptimizerType::New();

  // Configure metric
  // FIXME Delete lines afterchecking results with masks
  m_FixedMask = MaskType::New();
  m_FixedMask -> SetImage( this -> GetFixedImageMask() );
  m_Metric -> SetFixedImageMask( m_FixedMask );

//  std::cout << this -> GetFixedImageRegion().GetNumberOfPixels() << std::endl;

  // FIXME Uncomment after testing NC
//  m_Metric -> SetNumberOfHistogramBins( 24 );
//  m_Metric -> UseAllPixelsOn();
  m_Metric -> SetNumberOfSpatialSamples(0.2*this -> GetFixedImageRegion().GetNumberOfPixels());

  // Configure optimizer

  m_Optimizer->MinimizeOn();

  // FIXME Use MaximumStepLength = 1 for Initialization (3D) and 0.5 for SliceBySlice
  // Perhaps I can adjust these parameters taking into account the maximal expected
  // rotations and translations (+/- 10 deg and +/- 10 mm) . Check notes in my notebook
  // on 4th November
  m_Optimizer->SetMaximumStepLength( 0.1 );
  m_Optimizer->SetMinimumStepLength( 0.001 );
  m_Optimizer->SetNumberOfIterations( m_Iterations );
  m_Optimizer->SetRelaxationFactor( 0.8 );

  OptimizerScalesType optimizerScales( m_Transform -> GetNumberOfParameters() );

  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0/100.0;
  optimizerScales[4] =  1.0/100.0;
  optimizerScales[5] =  1.0/100.0;

  m_Optimizer->SetScales( optimizerScales );

  m_Observer = CommandIterationUpdate::New();

  if (m_EnableObserver)
  {
    m_Optimizer -> AddObserver( itk::IterationEvent(), m_Observer );
  }

/*  PointType rotationCenter;
  IndexType centerIndex;

  IndexType  fixedImageRegionIndex = this -> GetFixedImageRegion().GetIndex();
  SizeType   fixedImageRegionSize  = this -> GetFixedImageRegion().GetSize();

  centerIndex[0] = fixedImageRegionIndex[0] + fixedImageRegionSize[0] / 2.0;
  centerIndex[1] = fixedImageRegionIndex[1] + fixedImageRegionSize[1] / 2.0;
  centerIndex[2] = fixedImageRegionIndex[2] + fixedImageRegionSize[2] / 2.0;

  this -> GetFixedImage() -> TransformIndexToPhysicalPoint(centerIndex, rotationCenter);

  std::cout << "rotationIndex = " << centerIndex << std::endl;
  std::cout << "rotationPoint = " << rotationCenter << std::endl;

  m_Transform -> SetIdentity();
  m_Transform -> SetCenter( rotationCenter );

  */

  // Connect components
  SetTransform( this -> m_Transform );
  SetMetric( this -> m_Metric );
  SetOptimizer( this -> m_Optimizer );
  SetInterpolator( this -> m_Interpolator );

//  m_InitialTransformParameters = m_Transform -> GetParameters();
//  m_Registration->SetInitialTransformParameters( m_InitialTransformParameters );

  Superclass::Initialize();

}

/*
 * PrintSelf
 */
template < typename ImageType >
void
RigidRegistration<ImageType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


} // end namespace btk


#endif
