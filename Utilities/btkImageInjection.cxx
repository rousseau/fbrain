/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 04/08/2011
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

#include "btkLowToHighImageResolutionMethod.h"
#include "btkResampleImageByInjectionFilter.h"

#include "btkSliceBySliceTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"

#include "itkImage.h"
#include <tclap/CmdLine.h>

#include "itkImageMaskSpatialObject.h"
#include "btkImageIntersectionCalculator.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <stdio.h>

int main( int argc, char *argv[] )
{

  try {

  std::vector< std::string > input;
  std::vector< std::string > mask;
  std::vector< std::string > transform;
  std::vector< std::string > roi;

  double margin;

  const char *outImage = NULL;

  // Parse arguments

  TCLAP::CmdLine cmd("Creates a high resolution image from a set of low resolution images", ' ', "Unversioned");

  TCLAP::MultiArg<std::string> inputArg("i","input","Image file",true,"string",cmd);
  TCLAP::MultiArg<std::string> maskArg("m","","Mask file",false,"string",cmd);
  TCLAP::MultiArg<std::string> transformArg("t","transform","Transform input file",true,"string",cmd);
  TCLAP::ValueArg<std::string> outArg("o","output","High resolution image",true,"","string",cmd);

  TCLAP::ValueArg<double> marginArg("","margin","Adds a margin to the reconstructed images"
  "to compensate for a small FOV in the reference. The value must be provided in millimeters."
  " (default 0)",false, 0.0,"double",cmd);

  TCLAP::SwitchArg  boxSwitchArg("","box","Use intersections for roi calculation",false);
  TCLAP::SwitchArg  maskSwitchArg("","mask","Use masks for roi calculation",false);
  TCLAP::SwitchArg  allSwitchArg("","all","Use the whole image FOV",false);

  // xor arguments for roi assessment

  std::vector<TCLAP::Arg*>  xorlist;
  xorlist.push_back(&boxSwitchArg);
  xorlist.push_back(&maskSwitchArg);
  xorlist.push_back(&allSwitchArg);

  cmd.xorAdd( xorlist );

  // Parse the argv array.
  cmd.parse( argc, argv );

  input = inputArg.getValue();
  mask = maskArg.getValue();
  outImage = outArg.getValue().c_str();
  transform = transformArg.getValue();
  margin = marginArg.getValue();

  // typedefs

  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;

  typedef itk::Image< PixelType, Dimension >  ImageType;
  typedef ImageType::Pointer                  ImagePointer;

  typedef ImageType::RegionType               RegionType;
  typedef std::vector< RegionType >           RegionArrayType;

  typedef itk::ImageFileReader< ImageType  >  ImageReaderType;

  typedef itk::Image< unsigned char, Dimension >    ImageMaskType;
  typedef ImageMaskType::Pointer                    ImageMaskPointer;

  typedef itk::ImageFileReader< ImageMaskType >     MaskReaderType;
  typedef itk::ImageMaskSpatialObject< Dimension >  MaskType;
  typedef MaskType::Pointer  MaskPointer;

  typedef btk::LowToHighImageResolutionMethod<ImageType> LowToHighResFilterType;
  LowToHighResFilterType::Pointer lowToHighResFilter = LowToHighResFilterType::New();

  typedef btk::ResampleImageByInjectionFilter< ImageType, ImageType >  ResamplerType;

  typedef btk::SliceBySliceTransform< double, Dimension > TransformType;
  typedef TransformType::Pointer                          TransformPointer;

  typedef btk::Euler3DTransform< double > Euler3DTransformType;
  typedef Euler3DTransformType::Pointer   Euler3DTransformPointer;

  typedef itk::TransformFileReader     TransformReaderType;
  typedef TransformReaderType::TransformListType * TransformListType;

  typedef btk::ImageIntersectionCalculator<ImageType> IntersectionCalculatorType;
  IntersectionCalculatorType::Pointer intersectionCalculator = IntersectionCalculatorType::New();

  // Filter setup
  unsigned int numberOfImages = input.size();
  std::vector< ImagePointer >         		images(numberOfImages);
  std::vector< ImageMaskPointer >     		imageMasks(numberOfImages);
  std::vector< TransformPointer >     		transforms(numberOfImages);
  std::vector< Euler3DTransformPointer >  euler3DTransforms(numberOfImages);
  std::vector< MaskPointer >          		masks(numberOfImages);
  std::vector< RegionType >           		rois(numberOfImages);

  ImagePointer refImage;
  ImagePointer recImage;

  // Create reference image

  lowToHighResFilter -> SetNumberOfImages(numberOfImages);
  lowToHighResFilter -> SetTargetImage( 0 );
  lowToHighResFilter -> SetMargin( margin );

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    // Read and set image

    ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input[i].c_str() );
    imageReader -> Update();
    images[i] = imageReader -> GetOutput();

    lowToHighResFilter -> SetImageArray(i, images[i] );

    if ( boxSwitchArg.isSet() )
    {
      intersectionCalculator -> AddImage( images[i] );
    }

  }

  if ( boxSwitchArg.isSet() )
  {
    intersectionCalculator -> Compute();
  }

  // Set roi according to the provided arguments

  for (unsigned int i=0; i<numberOfImages; i++)
  {

    if ( maskSwitchArg.isSet() )
    {
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader -> SetFileName( mask[i].c_str() );
      maskReader -> Update();
      imageMasks[i] = maskReader -> GetOutput();
      lowToHighResFilter -> SetImageMaskArray( i, imageMasks[i] );

      masks[i] = MaskType::New();
      masks[i] -> SetImage( imageMasks[i] );
      rois[i] = masks[i] -> GetAxisAlignedBoundingBoxRegion();

      lowToHighResFilter -> SetRegionArray( i, rois[i] );

    } else if ( boxSwitchArg.isSet() )
      {
        imageMasks[i] = intersectionCalculator -> GetImageMask(i);
        lowToHighResFilter -> SetImageMaskArray( i, imageMasks[i] );

        masks[i] = MaskType::New();
        masks[i] -> SetImage( imageMasks[i] );
        rois[i] = masks[i] -> GetAxisAlignedBoundingBoxRegion();
        lowToHighResFilter -> SetRegionArray( i, rois[i] );

      } else if ( allSwitchArg.isSet() )
        {
          rois[i] = images[i] -> GetLargestPossibleRegion();

          lowToHighResFilter -> SetRegionArray( i, rois[i] );
        }
  }

  std::cout << "Creating isotropic image ... ";

  // Start registration
  try
    {
    lowToHighResFilter->StartRegistration();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "done. " << std::endl;

  lowToHighResFilter -> Initialize();
  refImage = lowToHighResFilter->GetHighResolutionImage();

  // Read transforms
  itk::TransformFactory<TransformType>::RegisterTransform();
  itk::TransformFactory<Euler3DTransformType>::RegisterTransform();

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    if (transform.size() > 0)
    {
      TransformReaderType::Pointer transformReader = TransformReaderType::New();
      transformReader -> SetFileName( transform[i] );
      transformReader -> Update();

      TransformListType transformList = transformReader->GetTransformList();
      TransformReaderType::TransformListType::const_iterator titr = transformList->begin();

      const char * className = titr -> GetPointer() -> GetNameOfClass();

      if (strcmp(className,"Euler3DTransform") == 0)
      {
        euler3DTransforms[i] = dynamic_cast< Euler3DTransformType * >( titr->GetPointer() );

        transforms[i] = TransformType::New();
        transforms[i] -> SetImage( images[i] );
        transforms[i] -> Initialize( euler3DTransforms[i] );
      } else if (strcmp(className,"SliceBySliceTransform") == 0)
        {
          transforms[i] = dynamic_cast< TransformType * >( titr->GetPointer() );
          transforms[i] -> SetImage( images[i]);
        } else
          {
            std::cout << className << " is not a valid transform. Exiting ..."
                << std::endl;
            return EXIT_FAILURE;
          }
    }
  }

  // Inject images

  std::cout << "Injecting images ... "; std::cout.flush();

  ResamplerType::Pointer resampler = ResamplerType::New();

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    resampler -> AddInput( images[i] );
    resampler -> AddRegion( rois[i] );
    resampler -> SetTransform(i, transforms[i] ) ;
  }

  resampler -> UseReferenceImageOn();
  resampler -> SetReferenceImage( refImage );
  resampler -> SetImageMask(lowToHighResFilter -> GetImageMaskCombination());
  resampler -> Update();

  std::cout << "done." << std::endl;

  recImage = resampler -> GetOutput() ;


  // Write HR image

  typedef itk::ImageFileWriter< ImageType >  WriterType;

  WriterType::Pointer writer =  WriterType::New();
  writer-> SetFileName( outImage );
  writer-> SetInput( recImage );

  std::cout << "Writing " << outImage << " ... ";
  writer-> Update();
  std::cout << "done." << std::endl;

  return EXIT_SUCCESS;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}

