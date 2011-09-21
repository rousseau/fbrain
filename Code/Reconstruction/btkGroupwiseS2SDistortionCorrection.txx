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

#ifndef _btkGroupwiseS2SDistortionCorrection_txx
#define _btkGroupwiseS2SDistortionCorrection_txx

#include "btkGroupwiseS2SDistortionCorrection.h"

#include "vnl/vnl_math.h"
#include "vnl/vnl_inverse.h"

namespace btk
{

/*
 * Constructor
 */
template < typename TSequence >
GroupwiseS2SDistortionCorrection<TSequence>
::GroupwiseS2SDistortionCorrection()
{
  m_Input = 0;
  m_Output = 0;
  m_TransformArray.resize(0);
  m_Interpolator = 0;
  m_FixedImage = 0;
  m_T2epi = 0;
  m_Resample = 0;

  m_FixedImageRegionDefined = false;

  m_Iterations = 2;

}

template < typename TSequence >
void
GroupwiseS2SDistortionCorrection<TSequence>
::SetGradientTable( const char* input )
{

  SequenceSizeType imageSize = m_Input -> GetLargestPossibleRegion().GetSize();
  m_GradientTable.set_size( imageSize[3],3 );

  // Fill gradient table in cartesian coordinates

  FILE* fr;
  fr = fopen( input, "r" );

  float value;
  unsigned ncol = 1;
  unsigned nrow = 1;

  for (unsigned int j=1; j <= 3*imageSize[3]; j++)
  {
    if ( ( j % imageSize[3] ) == 0 )
    {
      fscanf( fr, "%f\n", &value);
      m_GradientTable(nrow-1,ncol-1)=value;
      nrow = 0;
      ncol++;
    } else
    {
      fscanf( fr, "%f ", &value);
      m_GradientTable(nrow-1,ncol-1)=value;
    }
    nrow++;
  }

  fclose (fr);

  std::cout << m_GradientTable << std::endl;

}


template < typename TSequence >
void
GroupwiseS2SDistortionCorrection<TSequence>
::WriteTransforms( const char* folder )
{

  for (unsigned int i=0; i<m_TransformArray.size(); i++)
  {

    TransformWriterType::Pointer transformWriter = TransformWriterType::New();

    char fullTrafoName[255]; strcpy ( fullTrafoName, folder );
    char trafoName[255];

    sprintf ( trafoName, "/%d.txt", i+1 );
    strcat ( fullTrafoName,trafoName );

    transformWriter -> SetFileName ( fullTrafoName );

    for(unsigned int j = 0; j < m_TransformArray[i].size(); j++)
    {
      transformWriter -> AddTransform( m_TransformArray[i][j] );
    }

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
    }

  }
}

/**

 * Set the region of the fixed image to be considered for registration
 */
template < typename TSequence >
void
GroupwiseS2SDistortionCorrection<TSequence>
::SetFixedImageRegion( const RegionType & region )
{
  m_FixedImageRegion = region;
  m_FixedImageRegionDefined = true;
  this->Modified();
}

template < typename TSequence >
TSequence *
GroupwiseS2SDistortionCorrection<TSequence>
::GetOutput( )
{
    // Copy T2 epi into output

    SequenceIteratorType outputIt( m_Output, m_Output->GetLargestPossibleRegion() );
    outputIt.GoToBegin();

    IteratorType t2epiIt( m_T2epi, m_T2epi -> GetLargestPossibleRegion() );

    for( t2epiIt.GoToBegin(); !t2epiIt.IsAtEnd(); ++t2epiIt, ++outputIt)
    {
      outputIt.Set( t2epiIt.Get() );
    }

    // Copy resampled gradient images to output

    SequenceRegionType movingRegion;

    SequenceSizeType movingSize = m_SequenceSize;
    movingSize[3] = 0;
    movingRegion.SetSize( movingSize );

    SequenceIndexType movingIndex = m_SequenceIndex;


    for (unsigned int i = 1; i <= m_GradientDirections; i++)
    {

      movingIndex[3] = i;
      movingRegion.SetIndex( movingIndex );

      ImageExtractorPointer movingExtractor = ImageExtractorType::New();
      movingExtractor -> SetExtractionRegion( movingRegion );
      movingExtractor -> SetInput( m_Input );
      movingExtractor -> SetDirectionCollapseToSubmatrix();
      movingExtractor -> Update();

      ImagePointer movingImage = movingExtractor -> GetOutput();

      RBFInterpolatorPointer interpolator = RBFInterpolatorType::New();
      interpolator -> SetInputImage( movingImage );
      interpolator -> SetTransformArray( m_TransformArray[i-1] );
      interpolator -> SetRspa( 1.0 );
      interpolator -> Initialize( m_FixedImageRegion );


//      std::cout << "transform parametes for image " << i << std::endl;
//      for (unsigned t=0; t<m_TransformArray[i-1].size(); t++)
//      {
//        std::cout << m_TransformArray[i-1][t] -> GetParameters() << std::endl;
//      }

//      for( meanGradientIt.GoToBegin(); !meanGradientIt.IsAtEnd(); ++meanGradientIt)
//        {
//          index = meanGradientIt.GetIndex();
//          m_MeanGradient -> TransformIndexToPhysicalPoint( index, point );
//          meanGradientIt.Set( meanGradientIt.Get() + interpolator -> Evaluate(point) );
//        }
//
//
//      ResamplePointer resample = ResampleType::New();
//      resample->SetTransform( m_TransformArray[i-1] );
//
//      resample->SetInput( movingImage );
//      resample->SetInterpolator(m_Interpolator);
//
//      resample->SetSize( m_FixedImage -> GetLargestPossibleRegion().GetSize() );
//      resample->SetOutputOrigin(  m_FixedImage -> GetOrigin() );
//      resample->SetOutputSpacing( m_FixedImage -> GetSpacing() );
//      resample->SetOutputDirection( m_FixedImage -> GetDirection() );
//      resample->SetDefaultPixelValue( 100 );
//      resample->Update();
//
      IteratorType movingIt( movingImage, movingImage -> GetLargestPossibleRegion() );
      PointType point;
      IndexType index;

      for( movingIt.GoToBegin(); !movingIt.IsAtEnd(); ++movingIt, ++outputIt)
      {
        index = movingIt.GetIndex();
        if (m_FixedImageRegion.IsInside(index))
        {
          movingImage -> TransformIndexToPhysicalPoint( index, point );
          outputIt.Set( interpolator -> Evaluate(point) );
        }


      }

  }

  return this->m_Output;
}


template < typename TSequence >
void
GroupwiseS2SDistortionCorrection<TSequence>
::Initialize()
{

  if (!m_Input)
  {
      std::cout << "Input sequence is not present" << std::endl;
  }

  if (!m_FixedImageRegionDefined)
  {
      std::cout << "Fixed image region is not present" << std::endl;
  }

  m_Interpolator  = InterpolatorType::New();

  // Sequence region information

  m_SequenceRegion 		 = m_Input -> GetLargestPossibleRegion();
  m_SequenceIndex 		 = m_SequenceRegion.GetIndex();
  m_SequenceSize 			 = m_SequenceRegion.GetSize();
  m_GradientDirections = m_SequenceSize[3] - 1;

  // Extract T2 epi to be used as final reference

  SequenceSizeType extractorSize = m_SequenceSize;
  extractorSize[3] = 0;

  SequenceIndexType extractorIndex = m_SequenceIndex;
  extractorIndex[3] = 0;

  SequenceRegionType extractorRegion;
  extractorRegion.SetSize(  extractorSize  );
  extractorRegion.SetIndex( extractorIndex );

  ImageExtractorPointer T2epiExtractor  = ImageExtractorType::New();
  T2epiExtractor -> SetExtractionRegion( extractorRegion );
  T2epiExtractor -> SetInput( m_Input );
  T2epiExtractor -> SetDirectionCollapseToSubmatrix();
  T2epiExtractor -> Update();

  m_T2epi = T2epiExtractor -> GetOutput();

  // Extract G1 to be used as initial reference

  extractorIndex[3] = 1;
  extractorRegion.SetIndex( extractorIndex );

  ImageExtractorPointer G1Extractor  = ImageExtractorType::New();
  G1Extractor -> SetExtractionRegion( extractorRegion );
  G1Extractor -> SetDirectionCollapseToSubmatrix();
  G1Extractor -> SetInput( m_Input );

  // Smooth image to get a smoother cost function

  GaussianFilterPointer fixedSmoother  = GaussianFilterType::New();
  fixedSmoother -> SetUseImageSpacingOn();
  fixedSmoother -> SetVariance( 2.0 );
  fixedSmoother -> SetInput(   G1Extractor -> GetOutput() );
  fixedSmoother -> Update();

  m_FixedImage = fixedSmoother -> GetOutput();

  // Create registration object

  m_Registration = RegistrationType::New();

  // Set fixed image region
  // Here we swap fixed and moving images to perform slice by slice registration
  // and then we compute the transform inverse. Perhaps it is better to change
  // the names.

  m_Registration -> SetMovingImage( m_FixedImage );
  m_Registration -> SetFixedImageRegion( m_FixedImageRegion );

  // Create resampler

  m_Resample = ResampleType::New();

  m_Resample -> SetSize( m_FixedImage->GetLargestPossibleRegion().GetSize() );
  m_Resample -> SetOutputOrigin(  m_FixedImage->GetOrigin() );
  m_Resample -> SetOutputSpacing( m_FixedImage->GetSpacing() );
  m_Resample -> SetOutputDirection( m_FixedImage->GetDirection() );
  m_Resample -> SetDefaultPixelValue( 100 );

  // Create empty mean gradient image

  m_MeanGradient = ImageType::New();
  m_MeanGradient -> SetRegions( m_FixedImage -> GetLargestPossibleRegion() );
  m_MeanGradient -> Allocate();
  m_MeanGradient -> FillBuffer(0);
  m_MeanGradient -> SetOrigin( m_FixedImage -> GetOrigin() );
  m_MeanGradient -> SetDirection( m_FixedImage -> GetDirection() );
  m_MeanGradient -> SetSpacing( m_FixedImage -> GetSpacing() );

  // Create empty output sequence

  m_Output = SequenceType::New();
  m_Output -> SetRegions( m_Input -> GetLargestPossibleRegion() );
  m_Output -> Allocate();
  m_Output -> FillBuffer(0);
  m_Output -> SetOrigin( m_Input -> GetOrigin() );
  m_Output -> SetDirection( m_Input -> GetDirection() );
  m_Output -> SetSpacing( m_Input -> GetSpacing() );

   // Setup transforms

  PointType centerPoint;
  IndexType centerIndex;

  IndexType fixedImageRegionIndex = m_FixedImageRegion.GetIndex();
  SizeType  fixedImageRegionSize  = m_FixedImageRegion.GetSize();

  centerIndex[0] = fixedImageRegionIndex[0] + fixedImageRegionSize[0] / 2.0;
  centerIndex[1] = fixedImageRegionIndex[1] + fixedImageRegionSize[1] / 2.0;
  centerIndex[2] = fixedImageRegionIndex[2] + fixedImageRegionSize[2] / 2.0;

  m_FixedImage -> TransformIndexToPhysicalPoint(centerIndex, centerPoint);

  m_TransformArray.resize(m_GradientDirections);

}

/*
 * Starts the Registration Process
 */
template < typename TSequence >
void
GroupwiseS2SDistortionCorrection<TSequence>
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

//  TransformType::MatrixType meanAffine;
//  TransformType::MatrixType invMeanAffine;
//
//  // previousTransforms is used for assessing convergence
//
//  TransformPointerArray previousTransforms( m_GradientDirections );
//
//  for (unsigned int i=0; i<m_GradientDirections; i++)
//  {
//    previousTransforms[i] = TransformType::New();
//    previousTransforms[i] -> SetIdentity();
//    previousTransforms[i] -> SetCenter( m_TransformArray[i] -> GetCenter() );
//  }


  unsigned short nrep = 0;
  double transformError = 0.0;
//    double currentError  = 1e+10;

  SequenceRegionType movingRegion;

  SequenceSizeType movingSize = m_SequenceSize;
  movingSize[3] = 0;
  movingRegion.SetSize( movingSize );

  SequenceIndexType movingIndex = m_SequenceIndex;

  IteratorType meanGradientIt( m_MeanGradient, m_FixedImageRegion );
  IteratorType fixedIt( m_FixedImage, m_FixedImage -> GetLargestPossibleRegion() );

  do
  {
    m_MeanGradient -> FillBuffer(0.0);

    for (unsigned int i = 1; i <= m_GradientDirections; i++)
    {

      movingIndex[3] = i;
      movingRegion.SetIndex( movingIndex );

      ImageExtractorPointer movingExtractor = ImageExtractorType::New();
      movingExtractor -> SetExtractionRegion( movingRegion );
      movingExtractor -> SetInput( m_Input );
      movingExtractor -> SetDirectionCollapseToSubmatrix();

      GaussianFilterPointer movingSmoother  = GaussianFilterType::New();
      movingSmoother -> SetUseImageSpacingOn();
      movingSmoother -> SetVariance( 2.0 );
      movingSmoother -> SetInput( movingExtractor -> GetOutput() );
      movingSmoother -> Update();

      ImagePointer movingImage = movingSmoother -> GetOutput();

      m_Registration -> SetFixedImage( movingImage );

      if (nrep>0)
      {
        m_Registration -> SetTransformArray( m_TransformArray[i-1] );
      }
      std::cout << "Registering diffusion image " << i << " (iter " << nrep
          << ") ... "; std::cout.flush();

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

      m_TransformArray[i-1] = m_Registration -> GetTransformArray();

      IndexType index;
      PointType point;
      PointType transformedPoint;

      RBFInterpolatorPointer interpolator = RBFInterpolatorType::New();
      interpolator -> SetInputImage( movingImage );
      interpolator -> SetTransformArray( m_Registration -> GetTransformArray() );
      interpolator -> SetRspa( 1.0 );
      interpolator -> Initialize( m_FixedImageRegion );

      for( meanGradientIt.GoToBegin(); !meanGradientIt.IsAtEnd(); ++meanGradientIt)
        {
          index = meanGradientIt.GetIndex();
          m_MeanGradient -> TransformIndexToPhysicalPoint( index, point );
          meanGradientIt.Set( meanGradientIt.Get() + interpolator -> Evaluate(point) );
        }
    }

    for( meanGradientIt.GoToBegin(); !meanGradientIt.IsAtEnd(); ++meanGradientIt)
      {
        meanGradientIt.Set( meanGradientIt.Get() / m_GradientDirections );
      }

    // Update transforms for initialization and compute difference in transforms

    transformError = 0.0;

//    for (unsigned int i = 1; i <= m_GradientDirections; i++)
//    {
//
//      ParametersType previousParameters = previousTransforms[i-1] -> GetParameters();
//      ParametersType currentParameters  = m_TransformArray[i-1] -> GetParameters();
//
//      for (unsigned p = 0; p < m_TransformArray[i-1] -> GetNumberOfParameters(); p++)
//      {
//        transformError += std::abs( currentParameters[p] - previousParameters[p]) * m_Optimizer -> GetScales()[p];
//      }
//
//
//      previousTransforms[i-1] -> SetParameters( m_TransformArray[i-1] -> GetParameters() );
//
//      m_TransformArray[i-1] -> SetMatrix( m_TransformArray[i-1] -> GetMatrix() * invMeanAffine );
//      m_TransformArray[i-1] -> SetOffset( m_TransformArray[i-1] -> GetOffset() );
//    }
//
//    transformError /= ( (double)m_GradientDirections*(double)m_TransformArray[0] -> GetNumberOfParameters() );
//
//    std::cout << "transformError = " << transformError << std::endl;

    // Normalized becomes the new fixed

    IteratorType meanGradientIt( m_MeanGradient,
                                 m_MeanGradient -> GetLargestPossibleRegion() );
    for( meanGradientIt.GoToBegin(), fixedIt.GoToBegin();
        !meanGradientIt.IsAtEnd();
        ++meanGradientIt, ++fixedIt)
    {
        fixedIt.Set( meanGradientIt.Get() );
    }

    // Compute errors in parameters

    nrep++;


  } while ( (nrep < m_Iterations) ); //&& ( transformError > 0.005 ) );

  // Registration with T2 epi

  std::cout << "Registering mean gradient to B0 (T:mean -> B0) ... "; std::cout.flush();
//
  m_AffineRegistration = AffineRegistrationType::New();
  m_AffineRegistration -> SetFixedImage( m_T2epi );
  m_AffineRegistration -> SetMovingImage( m_MeanGradient );
  m_AffineRegistration -> SetFixedImageRegion( m_FixedImageRegion );
  m_AffineRegistration -> SetIterations( 0 );

  ParametersType initialParameters( 12 );
  initialParameters.Fill(0.0);
  initialParameters[0] = 1.0;
  initialParameters[4] = 1.0;
  initialParameters[8] = 1.0;

  m_AffineRegistration -> SetInitialTransformParameters( initialParameters );

  try
  {
    m_AffineRegistration -> StartRegistration();
  }
  catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      throw err;
    }

  std::cout << " done." << std::endl; std::cout.flush();

  std::cout << "final parameters = " << m_AffineRegistration -> GetTransform() -> GetParameters() << std::endl;
  std::cout << "center = " << m_AffineRegistration -> GetTransform() -> GetCenter() << std::endl;

  // Check registration mean grad to dti

  typename ResampleType::Pointer resampler = ResampleType::New();

  resampler -> SetUseReferenceImage( true );
  resampler -> SetReferenceImage(  m_T2epi );
  resampler -> SetDefaultPixelValue( 100 );
  resampler  -> SetTransform( m_AffineRegistration -> GetTransform() );
  resampler  -> SetInput( m_MeanGradient );
  resampler  -> Update();

  typename ImageWriterType::Pointer imageWriter =  ImageWriterType::New();

  imageWriter->SetFileName( "/home/miv/oubel/devel/RBFInterpolation/TAB_Ar/meanGradient2dti.nii.gz" );
  imageWriter->SetInput( resampler -> GetOutput() );
  imageWriter->Update();


  // Invert affine tranform and composition with transformations with slices

  TransformType::Pointer inverseTransform = TransformType::New();
  inverseTransform -> SetIdentity();
  inverseTransform -> SetCenter( m_AffineRegistration -> GetTransform() -> GetCenter() );
  inverseTransform -> SetParameters( m_AffineRegistration -> GetTransform() -> GetParameters() );
  inverseTransform -> GetInverse( inverseTransform );

  std::cout << "final parameters (inverse) = " << inverseTransform -> GetParameters() << std::endl;
  std::cout << "center  (inverse)  = " << inverseTransform -> GetCenter() << std::endl;

  // Final combination of transform to map everything into T2 epi

  MatrixType M2 = inverseTransform -> GetMatrix();
  OffsetType O2 = inverseTransform -> GetOffset();

  for (unsigned int i = 0; i < m_TransformArray.size(); i++)
  {
    for (unsigned int j = 0; j < m_TransformArray[i].size(); j++)
    {
      MatrixType M1 = m_TransformArray[i][j] -> GetMatrix();
      OffsetType O1 = m_TransformArray[i][j] -> GetOffset();

      m_TransformArray[i][j] -> SetMatrix( M2 * M1  );
      m_TransformArray[i][j] -> SetOffset( M2 * O1 + O2 );
    }
  }

}


/*
 * PrintSelf
 */
template < typename TSequence >
void
GroupwiseS2SDistortionCorrection<TSequence>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Member function to be implemented. " << std::endl;
}


} // end namespace btk


#endif
