/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 14/04/2010
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

#ifndef __BTK_AFFINEREGISTRATION_TXX__
#define __BTK_AFFINEREGISTRATION_TXX__

#include "btkAffineRegistration.h"

namespace btk
{

/*
 * Constructor
 */
template < typename ImageType >
AffineRegistration<ImageType>
::AffineRegistration()
{
  m_Iterations = 300;
  m_EnableObserver = false;
  m_FixedImageMask = 0;
}

/*
 * Initialize by setting the interconnects between components.
 */
template < typename ImageType >
void
AffineRegistration<ImageType>
::Initialize() throw (ExceptionObject)
{

  m_Transform = TransformType::New();
  m_Interpolator = InterpolatorType::New();
  m_Metric = MetricType::New();
  m_Optimizer = OptimizerType::New();

  // Configure metric

  if (!m_FixedImageMask)
  {
      m_FixedMask = MaskType::New();
      m_FixedMask -> SetImage(m_FixedImageMask );
      m_Metric->SetFixedImageMask(m_FixedMask);
  }
  else
  {
      m_FixedMask = MaskType::New();
      m_FixedMask -> SetImage( this -> GetFixedImageMask() );
      m_Metric -> SetFixedImageMask( m_FixedMask );
  }
  m_Metric->SetNumberOfHistogramBins( 24 );
  m_Metric->UseAllPixelsOn();

  // Configure optimizer

  m_Optimizer->MinimizeOn();
  m_Optimizer->SetMaximumStepLength( 0.2 );
  m_Optimizer->SetMinimumStepLength( 0.001 );
  m_Optimizer->SetNumberOfIterations( m_Iterations );
  m_Optimizer->SetRelaxationFactor( 0.8 );

  OptimizerScalesType optimizerScales( m_Transform -> GetNumberOfParameters() );

  optimizerScales[0] =  1.0;
  optimizerScales[1] =  1.0;
  optimizerScales[2] =  1.0;
  optimizerScales[3] =  1.0;
  optimizerScales[4] =  1.0;
  optimizerScales[5] =  1.0;
  optimizerScales[6] =  1.0;
  optimizerScales[7] =  1.0;
  optimizerScales[8] =  1.0;
  optimizerScales[9] =  1.0 / 1000.0;
  optimizerScales[10] =  1.0 / 1000.0;
  optimizerScales[11] =  1.0 / 1000.0;

  m_Optimizer->SetScales( optimizerScales );

  m_Observer = CommandIterationUpdate::New();

  if (m_EnableObserver)
  {
    m_Optimizer -> AddObserver( itk::IterationEvent(), m_Observer );
  }

  // Configure transform

  PointType rotationCenter;
  IndexType centerIndex;

  IndexType  fixedImageRegionIndex = this -> GetFixedImageRegion().GetIndex();
  SizeType   fixedImageRegionSize  = this -> GetFixedImageRegion().GetSize();

  centerIndex[0] = fixedImageRegionIndex[0] + fixedImageRegionSize[0] / 2.0;
  centerIndex[1] = fixedImageRegionIndex[1] + fixedImageRegionSize[1] / 2.0;
  centerIndex[2] = fixedImageRegionIndex[2] + fixedImageRegionSize[2] / 2.0;

  this -> GetFixedImage() -> TransformIndexToPhysicalPoint(centerIndex, rotationCenter);

  m_Transform -> SetIdentity();
  m_Transform -> SetCenter( rotationCenter );

  // Connect components

  this->SetTransform( this -> m_Transform );
  this->SetMetric( this -> m_Metric );
  this->SetOptimizer( this -> m_Optimizer );
  this->SetInterpolator( this -> m_Interpolator );

  if (this -> GetInitialTransformParameters().GetSize() == 1)
    this -> SetInitialTransformParameters( m_Transform -> GetParameters() );

  Superclass::Initialize();

}

/*
 * PrintSelf
 */
template < typename ImageType >
void
AffineRegistration<ImageType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
}


} // end namespace btk


#endif
