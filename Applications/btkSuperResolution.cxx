/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 02/12/2010
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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkEuler3DTransform.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageMaskSpatialObject.h"

#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "btkSliceBySliceTransform.h"

#include "btkSuperResolutionImageFilter.h"

#include "CmdLine.h"


int main( int argc, char *argv[] )
{

  try {

  std::vector< std::string > input;
  std::vector< std::string > mask;
  std::vector< std::string > transform;

  const char *outImage = NULL;
  const char *refImage = NULL;

  std::vector< int > x1, y1, z1, x2, y2, z2;

  unsigned int iter;
  float lambda;

  // Parse arguments

  TCLAP::CmdLine cmd("Resample a set of images using the injection method", ' ', "Unversioned");

  TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
  TCLAP::MultiArg<std::string> maskArg("m","mask","low-resolution image mask file",false,"string",cmd);
  TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image for initialization. "
      "Typically the output of btkImageReconstruction is used." ,true,"none","string",cmd);
  TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",true,"none","string",cmd);
  TCLAP::ValueArg<int> iterArg  ("","iter","Number of iterations",false, 10,"int",cmd);
  TCLAP::ValueArg<float> lambdaArg  ("","lambda","Regularization factor",false, 0.1,"float",cmd);

  TCLAP::MultiArg<std::string> transArg("t","transform","transform file",false,"string",cmd);


  // Parse the argv array.
  cmd.parse( argc, argv );

  input = inputArg.getValue();
  mask = maskArg.getValue();
  refImage = refArg.getValue().c_str();
  outImage = outArg.getValue().c_str();
  transform = transArg.getValue();
  iter = iterArg.getValue();
  lambda = lambdaArg.getValue();

  // typedefs
  const   unsigned int    Dimension = 3;
  typedef btk::SliceBySliceTransform< double, Dimension > TransformType;
  itk::TransformFactory<TransformType>::RegisterTransform();

  typedef float  PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;
  typedef std::vector<ImagePointer>           ImagePointerArray;

  typedef itk::Image< unsigned char, Dimension >  ImageMaskType;
  typedef itk::ImageFileReader< ImageMaskType > MaskReaderType;
  typedef itk::ImageMaskSpatialObject< Dimension > MaskType;

  typedef ImageType::RegionType               RegionType;
  typedef std::vector< RegionType >           RegionArrayType;

  typedef itk::ImageFileReader< ImageType >   ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >   WriterType;

  typedef itk::TransformFileReader     TransformReaderType;
  typedef TransformReaderType::TransformListType * TransformListType;

  typedef btk::SuperResolutionImageFilter< ImageType, ImageType >  ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();

  // Filter setup

  unsigned int numberOfImages = input.size();

  ImageType::IndexType  roiIndex;
  ImageType::SizeType   roiSize;

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    // add image
    ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input[i].c_str() );
    imageReader -> Update();
    resampler -> AddInput( imageReader -> GetOutput() );

    // add region
    if ( mask.size() > 0 )
    {
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader -> SetFileName( mask[i].c_str() );
      maskReader -> Update();

      MaskType::Pointer mask = MaskType::New();
      mask -> SetImage( maskReader -> GetOutput() );

      RegionType roi = mask -> GetAxisAlignedBoundingBoxRegion();
      roiIndex = roi.GetIndex();
      roiSize  = roi.GetSize();

    } else
        {
          roiSize  = imageReader -> GetOutput() -> GetLargestPossibleRegion().GetSize();
          roiIndex = imageReader -> GetOutput() -> GetLargestPossibleRegion().GetIndex();
        }

    RegionType imageRegion;
    imageRegion.SetIndex(roiIndex);
    imageRegion.SetSize (roiSize);
    resampler -> AddRegion( imageRegion );


    // Comment when testing with simulated images

    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader -> SetFileName( transform[i] );
    transformReader -> Update();

    TransformListType transforms = transformReader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator titr = transforms->begin();
    TransformType::Pointer trans = dynamic_cast< TransformType * >( titr->GetPointer() );

    for(unsigned int j=0; j< trans -> GetNumberOfSlices(); j++)
    {
      resampler -> SetTransform(i, j, trans -> GetSliceTransform(j) ) ;
      std::cout << trans -> GetSliceTransform(j) -> GetParameters() << " ; " << trans -> GetSliceTransform(j) -> GetFixedParameters() << std::endl;
    }
    std::cout << std::endl;


    // Finish comment when testing with simulated images

  }

  // Set reference image

  ImageReaderType::Pointer refReader = ImageReaderType::New();
  refReader -> SetFileName( refImage );
  refReader -> Update();

  resampler -> UseReferenceImageOn();
  resampler -> SetReferenceImage( refReader -> GetOutput() );
  resampler -> SetIterations(iter);
  resampler -> SetLambda( lambda );
  resampler -> Update();

  // Write image

  WriterType::Pointer writer =  WriterType::New();
  writer -> SetFileName( outImage );
  writer -> SetInput( resampler -> GetOutput() );

  if ( strcmp(outImage,"none") != 0)
  {
    std::cout << "Writing " << outImage << " ... ";
    writer->Update();
    std::cout << "done." << std::endl;
  }

  writer -> SetInput( resampler -> GetSimulatedImage(0) );
  writer -> SetFileName( "axial_simulated.nii.gz" );
  writer -> Update();

  writer -> SetInput( resampler -> GetSimulatedImage(1) );
  writer -> SetFileName( "coronal_simulated.nii.gz" );
  writer -> Update();

  writer -> SetInput( resampler -> GetSimulatedImage(2) );
  writer -> SetFileName( "sagital_simulated.nii.gz" );
  writer -> Update();



  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

