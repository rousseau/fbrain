/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 29/11/2011
  Author(s): François Rousseau (rousseau@unistra.fr)

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

//#include "vcl_iostream.h"
//#include "vnl/vnl_double_2.h"
//#include "vnl/vnl_cost_function.h"
//#include "vnl/algo/vnl_amoeba.h"


//#include "itkEuler3DTransform.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImageMaskSpatialObject.h"

#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkMatrixOffsetTransformBase.h"
#include "btkSliceBySliceTransform.h"

//#include "btkSliceBySliceTransform.h"

//#include "btkSuperResolutionImageFilter.h"

#include <tclap/CmdLine.h>
#include <stdio.h>

#include "btkSuperResolutionManager.h"


#include "itkBSplineInterpolateImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"


/*
class vnl_rosenbrock : public vnl_cost_function
{
public:
  vnl_rosenbrock(): vnl_cost_function(2) {}
  double y;
  double f(const vnl_vector<double>& x)
  {
    double u = 10*(x[1] - x[0]*x[0]);
    double v = 1 - x[0];
    return u*u + v*v;
  }  
};
*/


int main( int argc, char *argv[] )
{
  
  try {
    
    TCLAP::CmdLine cmd("Apply iterated back projection to high resolution image using one low resolution image", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file.", true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","Low-resolution image mask file.", true,"string",cmd);
    TCLAP::ValueArg<std::string> refArg  ("r","reconstructed","Current reconstructed image (input HR image). ",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outArg  ("o","output","Output image obtained using iterated by back-projection.", true,"","string",cmd);
    TCLAP::SwitchArg  boxcarSwitchArg("","boxcar","A boxcar-shaped PSF is assumed as imaging model (by default a Gaussian-shaped PSF is employed.).",false);
    TCLAP::MultiArg<std::string> transArg("t","transform","transform file",false,"string",cmd);
    TCLAP::ValueArg<int> loopArg  ("l","loop","Number of loops (iterations of IBP algorithm).", false,1,"int",cmd);
    TCLAP::ValueArg<int> nlmArg  ("n","nlm","Type of filtering during IBP process (0: no filtering, 1: error map filtering, 2: current HR image filtering).", false,0,"int",cmd);
    TCLAP::ValueArg<int> simArg  ("s","sim","Simulation of LR images based on the input HR image (0: no simulation, 1: simulation).", false,0,"int",cmd);
    
    // Parse the argv array.
    cmd.parse( argc, argv );
    
    // Get the value parsed by each arg. 
    std::vector<std::string> input_file       = inputArg.getValue();
    std::vector<std::string> mask_file        = maskArg.getValue();
    std::vector<std::string> transform_file   = transArg.getValue();

    std::string reference_file   = refArg.getValue();
    std::string output_file      = outArg.getValue();
    int loops                    = loopArg.getValue();
    int nlm                      = nlmArg.getValue();
    int simulation               = simArg.getValue();
    
    // typedefs
    const   unsigned int    Dimension = 3;
    typedef btk::SliceBySliceTransform< double, Dimension > TransformType;
    itk::TransformFactory<TransformType>::RegisterTransform();
    typedef itk::MatrixOffsetTransformBase<double, 3, 3> ANTSTransformType;
    itk::TransformFactory<ANTSTransformType>::RegisterTransform();
    
    
    typedef float PixelType;
    
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
    
    //typedef btk::SuperResolutionImageFilter< ImageType, ImageType >  ResamplerType;
    //ResamplerType::Pointer resampler = ResamplerType::New();
    
    //typedef btk::SliceBySliceTransform< double, 3 >     TransformType;
    //itk::TransformFactory<TransformType>::RegisterTransform();

    //-------------------------------------------------------
    SuperResolutionManager btkSRM;
    
    
    btkSRM.data.ReadHRImage(reference_file);
    btkSRM.data.ReadLRImages(input_file);
    btkSRM.data.ReadMaskLRImages(mask_file);
    btkSRM.data.ReadAffineTransform(transform_file);
    
    btkSRM.Initialize();
    
    if(simulation==1){
      btkSRM.SimulateLRImages();
      btkSRM.data.WriteSimulatedLRImages(input_file);
    }
    
    btkSRM.IteratedBackProjection(loops,nlm);
    btkSRM.data.WriteOutputHRImage(output_file);
    
    
    
    //-------------------------------------------------------
    /*
    typedef itk::BSplineInterpolateImageFunction<ImageType, double, double>    itkInterpolator;
    itkInterpolator::Pointer interpolator = itkInterpolator::New();
    interpolator->SetSplineOrder(1);
    
    typedef itk::ImageRegionIterator< ImageType > itkIterator;
    typedef itk::ImageRegionIteratorWithIndex< ImageType > itkIteratorWithIndex;
    
    itkIteratorWithIndex LRImageIt( btkSRM.data.m_inputLRImages[0], btkSRM.data.m_inputLRImages[0]->GetLargestPossibleRegion());
    
    ImageType::IndexType lrIndex;	
    ImageType::PointType lrPoint;	
    ImageType::PointType hrPoint;
    typedef itk::ContinuousIndex<double, ImageType::ImageDimension> ContinuousIndexType;

    ContinuousIndexType hrIndex;
    */
    /*
    for(LRImageIt.GoToBegin(); !LRImageIt.IsAtEnd(); ++LRImageIt){
      lrIndex = LRImageIt.GetIndex();
      
      btkSRM.data.m_inputLRImages[0]->TransformIndexToPhysicalPoint(lrIndex,lrPoint);
      std::cout<<"lr index:"<<lrIndex[0]<<","<<lrIndex[1]<<","<<lrIndex[2]<<" -- ";
      std::cout<<"lr point:"<<lrPoint[0]<<","<<lrPoint[1]<<","<<lrPoint[2]<<" \n ";
      
      hrPoint = btkSRM.data.m_affineTransform[0]->TransformPoint(lrPoint);
      std::cout<<"hr point:"<<hrPoint[0]<<","<<hrPoint[1]<<","<<hrPoint[2]<<" -- ";
      btkSRM.data.m_inputHRImage->TransformPhysicalPointToContinuousIndex(hrPoint, hrIndex);
      std::cout<<"hr index:"<<hrIndex[0]<<","<<hrIndex[1]<<","<<hrIndex[2]<<" \n ";
      
    }
    */
    
    
    /*
    // Filter setup
    
    ImageType::IndexType  roiIndex;
    ImageType::SizeType   roiSize;
    
    std::cout<<"Reading the input (low-resolution) image:"<<input_file<<std::endl;
    ImageReaderType::Pointer imageReader = ImageReaderType::New();
    imageReader -> SetFileName( input_file );
    imageReader -> Update();
    resampler   -> AddInput( imageReader -> GetOutput() );
    
    // add region
    if ( mask_file != "" )
    {
      std::cout<<"Reading the mask image:"<<mask_file<<std::endl;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader -> SetFileName( mask_file );
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
    
    
    if ( transform_file != "" )
    {
      std::cout<<"Reading the transform file : "<<transform_file<<std::endl;
      TransformReaderType::Pointer transformReader = TransformReaderType::New();
      transformReader -> SetFileName( transform_file );
      transformReader -> Update();
      
      TransformListType transforms = transformReader->GetTransformList();
      TransformReaderType::TransformListType::const_iterator titr = transforms->begin();
      TransformType::Pointer trans = dynamic_cast< TransformType * >( titr->GetPointer() );
      
      for(unsigned int j=0; j< trans -> GetNumberOfSlices(); j++)
        resampler -> SetTransform(0, j, trans -> GetSliceTransform(j) ) ;
    }
    
    
    std::cout<<"Reading the initial estimate of the high resolution image : "<<reference_file<<std::endl;
    ImageReaderType::Pointer refReader = ImageReaderType::New();
    refReader -> SetFileName( reference_file );
    refReader -> Update();
    resampler -> UseReferenceImageOn();
    resampler -> SetReferenceImage( refReader -> GetOutput() );
    ImagePointer HRima = refReader -> GetOutput();
    
    int n = 10;
    std::vector<ImagePointer> v(n);  
    for(int i=0; i<n; i++){
      ImagePointer ima = ImageType::New();
      ima->SetRegions(HRima->GetLargestPossibleRegion());
      ima->SetSpacing( HRima->GetSpacing() );
      ima->SetOrigin( HRima->GetOrigin() );
      ima->SetDirection( HRima->GetDirection() );
      ima->Allocate();
      ima->FillBuffer(0); 
      
      v[i] = HRima;  
      std::cout<<i<<" ";
      int toto=0;
      
      printf("toto: "); 
      scanf("%d", &toto);
      //char c;
      //scanf("%s",c); 
    }
    
    
    if ( boxcarSwitchArg.isSet() )
      resampler -> SetPSF( ResamplerType::BOXCAR );
    
    std::cout<<"Compute the matrix H (y = Hx)"<<std::endl;
    //resampler -> CreateH();
    
    // Write simulated image
    
    WriterType::Pointer writer =  WriterType::New();
    writer -> SetFileName( output_file );
    //writer -> SetInput( resampler -> GetSimulatedImage(0) );
    */
    /*
     if ( output_file != "")
     {
     std::cout << "Writing " << output_file << " ... ";
     writer -> Update();
     std::cout << "done." << std::endl;
     }
     */
    
    
    // Set up a Rosenbrock compute object
    /*
    vnl_rosenbrock f;
    
    // Set up the initial guess
    
    vnl_vector<double> x = vnl_double_2(-1.9,2.0).as_vector();
    
    // Make a Levenberg Marquardt minimizer, attach f to it, and
    
    // run the minimization
    f.y = 0;
    vnl_amoeba::minimize(f, x);
    
    // Summarize the results
    vcl_cout << "Rosenbrock min at " << x << '\n';
    */
    
    
    
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}

