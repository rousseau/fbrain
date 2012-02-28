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
    TCLAP::MultiArg<std::string> inputArg ("i","input","Low-resolution image file.", true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg  ("m","mask","Low-resolution image mask file.", true,"string",cmd);
    TCLAP::ValueArg<std::string> refArg   ("r","reconstructed","Current reconstructed image (input HR image). ",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outArg   ("o","output","Output image obtained using iterated by back-projection.", true,"","string",cmd);
    TCLAP::SwitchArg  boxcarSwitchArg     ("","boxcar","A boxcar-shaped PSF is assumed as imaging model (by default a Gaussian-shaped PSF is employed.).",false);
    TCLAP::MultiArg<std::string> transArg ("t","transform","transform file",false,"string",cmd);
    TCLAP::ValueArg<int> loopArg          ("l","loop","Maximum number of loops (iterations of IBP algorithm). Default=50.", false,50,"int",cmd);
    TCLAP::ValueArg<int> nlmArg           ("n","nlm","Type of filtering during IBP process (0: no filtering (default), 1: error map filtering, 2: current HR image filtering).", false,0,"int",cmd);
    TCLAP::ValueArg<float> betaArg        ("b","beta","Smoothing parameter for NLM filtering (default = 1).", false,1,"float",cmd);
    TCLAP::ValueArg<int> simArg           ("s","sim","Simulation of LR images based on the input HR image and the input LR images (0: no simulation, 1: simulation).", false,0,"int",cmd);
    TCLAP::ValueArg<int> ibpOrderArg      ("","ibpOrder","Order for the B-spline interpolation during image backprojections (0: nearest neighbor, 1: trilinear etc.)", false,5,"int",cmd);
    
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
    float beta                   = betaArg.getValue();
    int ibpOrder                 = ibpOrderArg.getValue();
    
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
    
    //-------------------------------------------------------
    SuperResolutionManager btkSRM;
    
    
    btkSRM.data.ReadHRImage(reference_file);
    btkSRM.data.ReadLRImages(input_file);
    btkSRM.data.ReadMaskLRImages(mask_file);
    btkSRM.data.ReadAffineTransform(transform_file);
    
    btkSRM.tool.SetPSFInterpolationOrderIBP(ibpOrder);
    
    btkSRM.Initialize();
    
    if(simulation==1){
      btkSRM.SimulateLRImages();
      btkSRM.data.WriteSimulatedLRImages(input_file);
    }
    
    btkSRM.IteratedBackProjection(loops,nlm,beta);
    btkSRM.data.WriteOutputHRImage(output_file);
        
    
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;
}

