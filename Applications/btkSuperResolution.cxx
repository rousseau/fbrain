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

/* Standard includes */
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkEuler3DTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkCastImageFilter.h"

/*Btk includes*/
#include "btkSliceBySliceTransform.h"
//#include "btkSuperResolutionImageFilter.h"
#include "btkSuperResolutionRigidImageFilter.h"

#include "btkNLMTool.h"


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

  TCLAP::CmdLine cmd("Apply super-resolution algorithm using one or multiple input images.", ' ', "Unversioned");

  TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
  TCLAP::MultiArg<std::string> maskArg("m","mask","low-resolution image mask file",false,"string",cmd);
  TCLAP::MultiArg<std::string> transArg("t","transform","transform file",false,"string",cmd);
  TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Reconstructed image for initialization. "
      "Typically the output of btkImageReconstruction is used." ,true,"","string",cmd);
  TCLAP::ValueArg<std::string> outArg  ("o","output","Super resolution output image",true,"","string",cmd);
  TCLAP::ValueArg<int> iterArg  ("","iter","Number of iterations (default = 20)",false, 20,"int",cmd);
  TCLAP::ValueArg<float> lambdaArg  ("","lambda","Regularization factor (default = 0.1)",false, 0.1,"float",cmd);
  TCLAP::SwitchArg  boxcarSwitchArg("","boxcar","A boxcar-shaped PSF is assumed as imaging model"
      " (by default a Gaussian-shaped PSF is employed.).",cmd,false);
  TCLAP::ValueArg<int> loopArg  ("","loop","Number of loops (SR/denoising) (default = 5)",false, 5,"int",cmd);
    

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

  typedef btk::SuperResolutionRigidImageFilter< ImageType, ImageType >  ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();

  // Filter setup

  unsigned int numberOfImages = input.size();

  ImageType::IndexType  roiIndex;
  ImageType::SizeType   roiSize;

  for (unsigned int i=0; i<numberOfImages; i++)
  {
    // add image
    std::cout<<"Reading image : "<<input[i].c_str()<<std::endl;
    ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input[i].c_str() );
    imageReader -> Update();
    resampler -> AddInput( imageReader -> GetOutput() );

    // add region
    if ( mask.size() > 0 )
    {
      std::cout<<"Reading mask image : "<<mask[i].c_str()<<std::endl;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader -> SetFileName( mask[i].c_str() );
      maskReader -> Update();

      MaskType::Pointer mask = MaskType::New();
      mask -> SetImage( maskReader -> GetOutput() );

      resampler -> AddMask( mask );

      RegionType roi = mask -> GetAxisAlignedBoundingBoxRegion();
      roiIndex = roi.GetIndex();
      roiSize  = roi.GetSize();

    } else
        {
          std::cout<<"Creating a mask image (entire input image)"<<std::endl;
          roiSize  = imageReader -> GetOutput() -> GetLargestPossibleRegion().GetSize();
          roiIndex = imageReader -> GetOutput() -> GetLargestPossibleRegion().GetIndex();
        }

    RegionType imageRegion;
    imageRegion.SetIndex(roiIndex);
    imageRegion.SetSize (roiSize);
    resampler -> AddRegion( imageRegion );


    if (transform.size() > 0)
    {
    std::cout<<"Reading transform:"<<transform[i]<<std::endl;
    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader -> SetFileName( transform[i] );
    transformReader -> Update();

    TransformListType transforms = transformReader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator titr = transforms->begin();
    TransformType::Pointer trans = dynamic_cast< TransformType * >( titr->GetPointer() );

    for(unsigned int j=0; j< trans -> GetNumberOfSlices(); j++)
      resampler -> SetTransform(i, j, trans -> GetSliceTransform(j) ) ;

    }

  }

  // Set reference image
    std::cout<<"Reading the reference image : "<<refImage<<std::endl;
  ImageReaderType::Pointer refReader = ImageReaderType::New();
  refReader -> SetFileName( refImage );
  refReader -> Update();

  std::cout<<"Performing super resolution"<<std::endl;
  resampler -> UseReferenceImageOn();
  resampler -> SetReferenceImage( refReader -> GetOutput() );
  resampler -> SetIterations(iter);
  resampler -> SetLambda( lambda );
  if ( boxcarSwitchArg.isSet() )
    resampler -> SetPSF( ResamplerType::BOXCAR );
  resampler -> Update();
	  
  int numberOfLoops = loopArg.getValue();
    
    
  for (int i=0; i<numberOfLoops; i++){
    std::cout<<"Loop : "<<i+1<<std::endl;
       
    btk::NLMTool<float> myTool;
    myTool.SetInput(resampler -> GetOutput());
    myTool.SetPaddingValue(0);
    myTool.SetDefaultParameters();
    myTool.ComputeOutput();
    
    ImagePointer outputImage = ImageType::New();
    myTool.GetOutput(outputImage);

    resampler -> SetReferenceImage( outputImage );
    resampler -> Update();
  }
  //NLM denoising desired at the last step if number of loops > 0
  if(numberOfLoops>0){
      
    btk::NLMTool<float> myTool;
    myTool.SetInput(resampler -> GetOutput());
    myTool.SetPaddingValue(0);
    myTool.SetDefaultParameters();
    myTool.ComputeOutput();
    
    ImagePointer outputImage = ImageType::New();
    myTool.GetOutput(outputImage);
        
    resampler -> SetReferenceImage( outputImage );    
  }

  // Write image

  WriterType::Pointer writer =  WriterType::New();
  writer -> SetFileName( outImage );
  writer -> SetInput( resampler -> GetOutput() );

  if ( strcmp(outImage,"") != 0)
  {
    std::cout << "Writing " << outImage << " ... ";
    writer->Update();
    std::cout << "done." << std::endl;
  }

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

