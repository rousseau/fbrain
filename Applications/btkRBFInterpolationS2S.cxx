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


int main( int argc, char *argv[] )
{

  try {

  const char *input = NULL, *ref = NULL, *output = NULL;
  const char *gtable = NULL, *rbf = NULL;
  const char *tpath = NULL;

  double sx, sy, sz;
  int x1, y1, z1, x2, y2, z2;
  double rspa, rgra;
// int x1hr, y1hr, z1hr, x2hr, y2hr, z2hr;

  TCLAP::CmdLine cmd("Resampling of dwi sequences using Radial Basis Functions", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Original sequence (.nii;.nii.gz)",true,"none","string");
  TCLAP::ValueArg<std::string> refArg("r","ref","Reference sequence (.nii;.nii.gz)",true,"none","string");
  TCLAP::ValueArg<std::string> outputArg("o","output","Corrected sequence (.nii;.nii.gz)",false,"none","string");
  TCLAP::ValueArg<std::string> gtableArg("g","gradients","Gradient table (.bvec)",true,"none","string");
  TCLAP::ValueArg<std::string> tpathArg("t","tpath","Path to transforms",true,"none","string");

  TCLAP::ValueArg<float> sxArg("","sx","New spacing in x",false,0,"float");
  TCLAP::ValueArg<float> syArg("","sy","New spacing in y",false,0,"float");
  TCLAP::ValueArg<float> szArg("","sz","New spacing in z",false,0,"float");

  TCLAP::ValueArg<int> x1Arg("","x1","ROI's min x value",false,0,"int");
  TCLAP::ValueArg<int> y1Arg("","y1","ROI's min y value",false,0,"int");
  TCLAP::ValueArg<int> z1Arg("","z1","ROI's min z value",false,0,"int");
  TCLAP::ValueArg<int> x2Arg("","x2","ROI's max x value",false,0,"int");
  TCLAP::ValueArg<int> y2Arg("","y2","ROI's max x value",false,0,"int");
  TCLAP::ValueArg<int> z2Arg("","z2","ROI's max x value",false,0,"int");

  TCLAP::ValueArg<std::string> rbfArg("","rbf","Radial basis function",false,"gauss","string");
  TCLAP::ValueArg<float> rspaArg("","rspa","Spatial kernel width",false,1,"float");
  TCLAP::ValueArg<float> rgraArg("","rgra","Angular kernel width",false,0.5,"float");

  cmd.add( z2Arg );
  cmd.add( y2Arg );
  cmd.add( x2Arg );
  cmd.add( z1Arg );
  cmd.add( y1Arg );
  cmd.add( x1Arg );

  cmd.add( szArg );
  cmd.add( syArg );
  cmd.add( sxArg );

  cmd.add( tpathArg );
  cmd.add( gtableArg );
  cmd.add( outputArg );
  cmd.add( refArg );
  cmd.add( inputArg );

  cmd.add( rbfArg );

  cmd.add( rspaArg );
  cmd.add( rgraArg );

  // Parse the argv array.
  cmd.parse( argc, argv );

  tpath = tpathArg.getValue().c_str();
  input = inputArg.getValue().c_str();
  ref = refArg.getValue().c_str();
  output = outputArg.getValue().c_str();
  gtable = gtableArg.getValue().c_str();

  sx = sxArg.getValue();
  sy = syArg.getValue();
  sz = szArg.getValue();

  // ROI on original image
  x1 = x1Arg.getValue();
  y1 = y1Arg.getValue();
  z1 = z1Arg.getValue();

  x2 = x2Arg.getValue();
  y2 = y2Arg.getValue();
  z2 = z2Arg.getValue();

  // Interpolation parameters
  rbf = rbfArg.getValue().c_str();
  rspa = rspaArg.getValue();
  rgra = rgraArg.getValue();


  // Read sequence

  typedef short PixelType;

  const unsigned int Dimension = 4;
  typedef itk::Image< PixelType, Dimension > SequenceType;

  typedef itk::ImageFileReader< SequenceType > SequenceReaderType;
  SequenceReaderType::Pointer sequenceReader = SequenceReaderType::New();
  sequenceReader -> SetFileName( input );
  sequenceReader -> Update();

  SequenceType::Pointer sequence = sequenceReader -> GetOutput();

  // Read reference

  SequenceReaderType::Pointer refSequenceReader = SequenceReaderType::New();
  refSequenceReader -> SetFileName( ref );
  refSequenceReader -> Update();

  SequenceType::Pointer refSequence = refSequenceReader -> GetOutput();


  // Create high resolution sequence (only in space by the moment)

  // High resolution spacing

  SequenceType::SpacingType spacing = sequence -> GetSpacing();

  SequenceType::SpacingType hrSpacing;

  if (sx > 0)
  {
    hrSpacing[0] = sx;
  } else
  {
    hrSpacing[0] = spacing[0];
  }

  if (sy > 0)
  {
    hrSpacing[1] = sy;
  } else
  {
    hrSpacing[1] = spacing[1];
  }

  if (sz > 0)
  {
    hrSpacing[2] = sz;
  } else
  {
    hrSpacing[2] = spacing[2];
  }

  hrSpacing[3] = spacing[3];

  // High resolution size

  SequenceType::SizeType size = sequence -> GetLargestPossibleRegion().GetSize();

  SequenceType::SizeType hrSize;
  hrSize[0] = ceil (size[0]*spacing[0]/hrSpacing[0]);
  hrSize[1] = ceil (size[1]*spacing[1]/hrSpacing[1]);
  hrSize[2] = ceil (size[2]*spacing[2]/hrSpacing[2]);
  hrSize[3] = ceil (size[3]*spacing[3]/hrSpacing[3]);

  // Create HR image

  SequenceType::Pointer hrSequence = SequenceType::New();

  SequenceType::IndexType hrStart;
  hrStart[0] = 0;
  hrStart[1] = 0;
  hrStart[2] = 0;
  hrStart[3] = 0;

  SequenceType::RegionType hrRegion;
  hrRegion.SetSize(hrSize);
  hrRegion.SetIndex(hrStart);

  hrSequence -> SetRegions(hrRegion);
  hrSequence -> Allocate();
  hrSequence -> FillBuffer(0.0);

  hrSequence -> SetDirection( sequence -> GetDirection() );
  hrSequence -> SetOrigin( sequence -> GetOrigin() );
  hrSequence -> SetSpacing( hrSpacing );

  // Resampling in ROI (testing purposes)
  // ROI definition on HR image

  // Compute points defining the ROI on the high resolution image

  SequenceType::IndexType LRIndex1;
  LRIndex1[0] = x1; LRIndex1[1] = y1; LRIndex1[2] = z1; LRIndex1[3] = 1;

  SequenceType::PointType LRPoint1;
  sequence -> TransformIndexToPhysicalPoint(LRIndex1, LRPoint1);

  SequenceType::IndexType HRIndex1;
  hrSequence -> TransformPhysicalPointToIndex(LRPoint1, HRIndex1);


  SequenceType::IndexType LRIndex2;
  LRIndex2[0] = x2; LRIndex2[1] = y2; LRIndex2[2] = z2; LRIndex2[3] = 1;

  SequenceType::PointType LRPoint2;
  sequence -> TransformIndexToPhysicalPoint(LRIndex2, LRPoint2);

  SequenceType::IndexType HRIndex2;
  hrSequence -> TransformPhysicalPointToIndex(LRPoint2, HRIndex2);

  // Create ROI in high resolution

  SequenceType::RegionType hrROI;

  SequenceType::IndexType roiStart = HRIndex1;

  SequenceType::SizeType roiSize;
  roiSize[0] = HRIndex2[0] - HRIndex1[0] + 1;
  roiSize[1] = HRIndex2[1] - HRIndex1[1] + 1;
  roiSize[2] = HRIndex2[2] - HRIndex1[2] + 1;
  roiSize[3] = 30;

  hrROI.SetIndex(roiStart);
  hrROI.SetSize(roiSize);

  // Create interpolator
  // WARNING : THE ORDER OF PARAMETER SETTING IS IMPORTANT !!!!

  typedef btk::RBFInterpolateImageFunctionS2S< SequenceType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator -> SetInputImage( sequence );
  interpolator -> SetTransforms(tpath);
  interpolator -> SetGradientTable( gtable );
  std::cout << "Initializing interpolator ... "; std::cout.flush();
  interpolator -> Initialize( hrROI );
  std::cout << "done. " << std::endl; std::cout.flush();

  typedef itk::ImageRegionIteratorWithIndex<SequenceType> IteratorType;

  IteratorType refIt(refSequence, hrROI);

  SequenceType::IndexType index;
  SequenceType::PointType point;

  double theta;
  double phi;

  double value;
  double error;

  clock_t start, finish;

  double mse = 0;
  double time = 0;

  start = clock();

  unsigned int image = 0;

  for (refIt.GoToBegin(); !refIt.IsAtEnd(); ++refIt)
  {
    index = refIt.GetIndex();

    if (index[3]!=image)
    {
      image = index[3];
      std::cout << "Resampling image " << image << " ... " << std::endl; std::cout.flush();
    }

// std::cout << "slice nro " << index[2] << std::endl;

    refSequence -> TransformIndexToPhysicalPoint(index,point);
    interpolator -> GetGradientDirection(index[3],theta,phi);

    value = interpolator -> Evaluate(point, theta, phi, rspa, rgra, 0);

    refIt.Set( (short)value );

  }
  finish = clock();

  // Write interpolated sequence

  typedef itk::ImageFileWriter< SequenceType > SequenceWriterType;
  SequenceWriterType::Pointer sequenceWriter = SequenceWriterType::New();
  sequenceWriter -> SetFileName( output );
  sequenceWriter -> SetInput( refSequence );

  if ( strcmp(output,"none") != 0 )
  {
    sequenceWriter -> Update();
  }

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e) // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

