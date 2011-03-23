/*==========================================================================

© Université de Strasbourg - Centre National de la Recherche Scientifique

Date: 23/03/2010
Author(s): Estanislao Oubel (oubel@unistra.fr)

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.

==========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImage.h"
#include "CmdLine.h"
#include "btkRBFInterpolateImageFunctionS2S.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "time.h"

#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"


int main( int argc, char *argv[] )
{

  try {

  const char *input = NULL, *refFile = NULL, *output = NULL;
  const char *maskFile = NULL, *trefFile = NULL;
  const char *gtable = NULL, *rbf = NULL;
  const char *tpath = NULL;

  double rspa, rgra;

  TCLAP::CmdLine cmd("Resampling of dwi sequences using Radial Basis Functions", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Original sequence (.nii;.nii.gz)",true,"","string",cmd);

  TCLAP::ValueArg<std::string> outputArg("o","output","Corrected sequence (.nii;.nii.gz)",true,"","string",cmd);
  TCLAP::ValueArg<std::string> gtableArg("g","gradients","Gradient table (.bvec)",true,"","string",cmd);
  TCLAP::ValueArg<std::string> tpathArg("t","tpath","Path to transforms",true,"none","string",cmd);
  TCLAP::ValueArg<std::string> maskArg("m","mask","Image mask in the resampling space (.nii;.nii.gz)",false,"","string",cmd);

  TCLAP::ValueArg<std::string> refArg("r","ref","Reference image to specify the resampling grid (.nii;.nii.gz)",false,"","string",cmd);
  TCLAP::ValueArg<std::string> trefArg("","tref","Transform from reference to sequence (.nii;.nii.gz)",false,"","string",cmd);

  TCLAP::ValueArg<std::string> rbfArg("","rbf","Radial basis function",false,"gauss","string",cmd);
  TCLAP::ValueArg<float> rspaArg("","rspa","Spatial kernel width",false,1,"float",cmd);
  TCLAP::ValueArg<float> rgraArg("","rgra","Angular kernel width",false,0.5,"float",cmd);

  TCLAP::SwitchArg oriSwitch("","originalSpace","Original resampling space", cmd, false);
  TCLAP::SwitchArg refSwitch("","referenceSpace","Resampling space provided by image reference", cmd, false);

  // Parse the argv array.
  cmd.parse( argc, argv );

  tpath = tpathArg.getValue().c_str();
  input  = inputArg.getValue().c_str();
  refFile  = refArg.getValue().c_str();
  output = outputArg.getValue().c_str();
  gtable = gtableArg.getValue().c_str();

  maskFile = maskArg.getValue().c_str();
  trefFile = trefArg.getValue().c_str();

  // Interpolation parameters
  rbf = rbfArg.getValue().c_str();
  rspa = rspaArg.getValue();
  rgra = rgraArg.getValue();

  // Read sequence

  typedef short         PixelType;

  const   unsigned int  Dimension = 4;
  typedef itk::Image< PixelType, Dimension >   SequenceType;
  typedef itk::ImageFileReader< SequenceType >  SequenceReaderType;

  SequenceReaderType::Pointer sequenceReader = SequenceReaderType::New();
  sequenceReader -> SetFileName( input );
  sequenceReader -> Update();
  SequenceType::Pointer       sequence     = sequenceReader -> GetOutput();

  // Read mask for DW sequence

  typedef itk::Image< unsigned char, 3 >   ImageMaskType;
  typedef itk::ImageFileReader< ImageMaskType >  ImageMaskReaderType;

  ImageMaskReaderType::Pointer imageMaskReader = ImageMaskReaderType::New();
  imageMaskReader -> SetFileName( maskFile );
  imageMaskReader -> Update();

  typedef itk::ImageMaskSpatialObject< 3 >  MaskType;
  MaskType::Pointer mask = MaskType::New();
  mask -> SetImage( imageMaskReader -> GetOutput() );

  // Calcule 4D diffusion ROI for input sequence from mask

  SequenceType::RegionType seqROI = sequence -> GetLargestPossibleRegion();
  SequenceType::IndexType  seqROIIndex = seqROI.GetIndex();
  SequenceType::SizeType   seqROISize  = seqROI.GetSize();;

  MaskType::IndexType seqROI3DIndex = mask -> GetAxisAlignedBoundingBoxRegion().GetIndex();
  MaskType::SizeType  seqROI3DSize  = mask -> GetAxisAlignedBoundingBoxRegion().GetSize();

  seqROIIndex[0] = seqROI3DIndex[0];
  seqROIIndex[1] = seqROI3DIndex[1];
  seqROIIndex[2] = seqROI3DIndex[2];
  seqROIIndex[3] = 0; // The roi starts at the first diffusion image

  seqROISize[0] = seqROI3DSize[0];
  seqROISize[1] = seqROI3DSize[1];
  seqROISize[2] = seqROI3DSize[2];
//  seqROISize[3] = seqROISize[3] - 1;
  seqROISize[3] = 2;

  seqROI.SetIndex(seqROIIndex);
  seqROI.SetSize(seqROISize);

  SequenceType::RegionType recROI = seqROI; // By default we reconstruct on the original grid

   // Create reconstructed sequence

  SequenceType::Pointer recSequence = SequenceType::New();

  SequenceType::SpacingType   recSpacing   = sequence -> GetSpacing();
  SequenceType::RegionType    recRegion    = sequence -> GetLargestPossibleRegion();
  SequenceType::SizeType      recSize      = recRegion.GetSize();
  SequenceType::IndexType     recIndex     = recRegion.GetIndex();
  SequenceType::PointType     recOrigin    = sequence -> GetOrigin();
  SequenceType::DirectionType recDirection = sequence -> GetDirection();

  typedef itk::Image< PixelType, 3 >   ImageType;
  typedef itk::ImageFileReader< ImageType >  ImageReaderType;

  if ( oriSwitch.isSet() )
  {

  } else if ( refSwitch.isSet() )
  {
    ImageReaderType::Pointer refImageReader = ImageReaderType::New();
    refImageReader -> SetFileName( refFile );
    refImageReader -> Update();

    ImageType::Pointer       refImage     = refImageReader -> GetOutput();
    ImageType::SpacingType   refSpacing   = refImage -> GetSpacing();
    ImageType::SizeType      refSize      = refImage -> GetLargestPossibleRegion().GetSize();
    ImageType::IndexType     refIndex     = refImage -> GetLargestPossibleRegion().GetIndex();
    ImageType::DirectionType refDirection = refImage -> GetDirection();
    ImageType::PointType     refOrigin    = refImage -> GetOrigin();

    recSpacing[0] = refSpacing[0];  recSpacing[1] = refSpacing[1]; recSpacing[2] = refSpacing[2];

    for(unsigned int i=0; i<3; i++)
      for(unsigned int j=0; j<3; j++)
        recDirection(i,j)= refDirection(i,j);

    recSize[0] = refSize[0];  recSize[1] = refSize[1]; recSize[2] = refSize[2];
    recIndex[0] = refIndex[0];  recIndex[1] = refIndex[1]; recIndex[2] = refIndex[2];

    recRegion.SetIndex( recIndex );
    recRegion.SetSize( recSize );

    recOrigin[0] = refOrigin[0];  recOrigin[1] = refOrigin[1]; recOrigin[2] = refOrigin[2];

  } else
    {
      recSize[2]    = floor( recSize[2]*recSpacing[2] / recSpacing[0] + 0.5 );
      recSpacing[2] = recSpacing[0];
      recRegion.SetSize( recSize );

      SequenceType::IndexType recROIIndex = recROI.GetIndex();
      SequenceType::SizeType  recROISize  = recROI.GetSize();

      recROIIndex[2] = floor (recROIIndex[2] * sequence -> GetSpacing()[2] / recSpacing[2] + 0.5);
      recROISize[2]  = floor (recROISize[2] * sequence -> GetSpacing()[2] / recSpacing[2] + 0.5);

      recROI.SetIndex( recROIIndex );
      recROI.SetSize( recROISize );
    }

  recSequence -> SetRegions( recRegion );
  recSequence -> Allocate();
  recSequence -> FillBuffer(0.0);

  recSequence -> SetDirection( recDirection);
  recSequence -> SetOrigin( recOrigin );
  recSequence -> SetSpacing( recSpacing );

  // Read transform ref -> dwi

  typedef itk::AffineTransform< double, 3 > TransformType;
  TransformType::Pointer tref;

  // Read transformation file
  if ( refSwitch.isSet() )
  {
    typedef itk::TransformFileReader     TransformReaderType;
    typedef TransformReaderType::TransformListType * TransformListType;

    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader->SetFileName( trefFile );
    transformReader->Update();

    TransformListType transforms = transformReader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator titr = transforms->begin();
    tref = dynamic_cast< TransformType * >( titr->GetPointer() ) ;
  }

  // Create interpolator
  // WARNING : THE ORDER OF PARAMETER SETTING IS IMPORTANT !!!!

  typedef btk::RBFInterpolateImageFunctionS2S< SequenceType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator -> SetInputImage( sequence );
  interpolator -> SetTransforms(tpath);
  interpolator -> SetGradientTable( gtable );
  interpolator -> Initialize( seqROI );

  typedef itk::ImageRegionIteratorWithIndex<SequenceType> IteratorType;

  IteratorType recIt( recSequence, recROI);

  SequenceType::IndexType index;
  SequenceType::PointType pointSeq;
  SequenceType::PointType pointRef;
  ImageType::PointType pointRef3D;
  ImageType::PointType pointSeq3D;

  double theta;
  double phi;
  double value;

  clock_t start, finish;

  start = clock();

  // First resample the T2 image
  typedef itk::LinearInterpolateImageFunction< SequenceType,
                                               double >      LinearInterpolatorType;
  LinearInterpolatorType::Pointer linearInterpolator = LinearInterpolatorType::New();
  linearInterpolator -> SetInputImage( sequence );

  for (recIt.GoToBegin(); !recIt.IsAtEnd(); ++recIt)
  {
    index = recIt.GetIndex();
    if (index[3] ==1) break;
    recSequence -> TransformIndexToPhysicalPoint(index,pointRef);
    recIt.Set( linearInterpolator -> Evaluate(pointRef) );
  }

  unsigned int image = 0;

  for (; !recIt.IsAtEnd(); ++recIt)
  {
    index = recIt.GetIndex();

    if (index[3]!=image)
    {
      image = index[3];
      std::cout << "Resampling image " << image << " ... " << std::endl; std::cout.flush();
    }

    value = 0;

    recSequence -> TransformIndexToPhysicalPoint(index,pointRef);
    pointRef3D[0] = pointRef[0];  pointRef3D[1] = pointRef[1]; pointRef3D[2] = pointRef[2];

    if ( mask -> IsInside(pointRef3D) )
    {
      if ( refSwitch.isSet() )
        pointSeq3D = tref -> TransformPoint( pointRef3D );
      else
        pointSeq3D = pointRef3D;

      pointSeq[0] = pointSeq3D[0];
      pointSeq[1] = pointSeq3D[1];
      pointSeq[2] = pointSeq3D[2];
      pointSeq[3] = pointRef[3];

      interpolator -> GetGradientDirection(index[3],theta,phi);

      value =  interpolator -> Evaluate(pointSeq, theta, phi, rspa, rgra, 0);
    }

    recIt.Set( (short)value );

  }
  finish = clock();

  std::cout << "duracion (seg) = " << (finish - start) / CLOCKS_PER_SEC << std::endl;

  // Write interpolated sequence

  typedef itk::ImageFileWriter< SequenceType >  SequenceWriterType;
  SequenceWriterType::Pointer sequenceWriter = SequenceWriterType::New();
  sequenceWriter -> SetFileName( output );
  sequenceWriter -> SetInput( recSequence );

  if ( strcmp(output,"none") != 0 )
  {
    sequenceWriter -> Update();
  }

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
