/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 09/08/2012
  Author(s): Youssef Taleb

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


/* ITK */
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"

/* BTK */
#include "btkImageHelper.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>


int main( int argc, char *argv[] )
{

    //TCLAP Commands for arguments

    TCLAP::CmdLine cmd("Resample input image (mask) in the space of reference image", ' ', "Unversioned");
    TCLAP::ValueArg<std::string> inputArg("i","input","Mask file",true,"","string",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","Fixed Image",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outArg  ("o","output","Extended mask file",true,"","string",cmd);


    std::string  inputMaskFile;
    std::string  outputMaskFile;
    std::string  refImageFile;



    // Parse the argv array.
    cmd.parse( argc, argv );
    inputMaskFile = inputArg.getValue();
    refImageFile  = refArg.getValue();
    outputMaskFile = outArg.getValue();

  if( argc < 4 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " -i movingImageMask -r fixedImageFile -o outputImageMask"<< std::endl;
    return EXIT_FAILURE;
    }


  // type definitions

  const    unsigned int    Dimension = 3;

  typedef  float                                                        PixelType;

  typedef itk::Image< PixelType, Dimension >                            ImageType;
  typedef itk::AffineTransform< double, Dimension >                     TransformType;
  typedef itk::Image<unsigned char, Dimension>                          ImageMaskType;
  typedef itk::ResampleImageFilter< ImageType,ImageType >               ResampleFilterType;
  typedef itk::CastImageFilter< ImageType,ImageMaskType >               CastFilterType;


 // Read Reference Image

  ImageType::Pointer  refImage;
  refImage = btk::ImageHelper<ImageType>::ReadImage(refImageFile);


  // Creation of an empty fixed image with the size of reference image
  ImageType::Pointer fixedImage = ImageType::New();
  fixedImage = btk::ImageHelper<ImageType>::CreateNewImageFromPhysicalSpaceOf(refImage);

  ImageType::Pointer movingImage = btk::ImageHelper<ImageType>::ReadImage(inputMaskFile);

  TransformType::Pointer identityTransform = TransformType::New();
  identityTransform->SetIdentity();

  // Extend input mask to the dimension of the fixed image with a resampler
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler -> SetInput( movingImage );
  resampler -> SetTransform( identityTransform );
  resampler -> UseReferenceImageOn();
  resampler -> SetReferenceImage( fixedImage );
  resampler -> Update();


  // Write output Mask into image
  CastFilterType::Pointer  caster =  CastFilterType::New();
  caster-> SetInput( resampler->GetOutput() );
  caster->Update();

  btk::ImageHelper< ImageMaskType >::WriteImage(caster->GetOutput(),outputMaskFile);

  return EXIT_SUCCESS;
}

