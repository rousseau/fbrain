/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 27/01/2012
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


/* Standard includes */
#include "string"
#include "iomanip"
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkContinuousIndex.h"


/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPandoraBoxImageFilters.h"

int main(int argc, char *argv[])
{
  try{

    TCLAP::CmdLine cmd("Resample an image using: 1) a reference image (better option) or 2) specific size or spacing.", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputImageArg("i","input_image","input image file (will be casted in float)",true,"","string");
    cmd.add( inputImageArg );
    TCLAP::ValueArg<std::string> outputImageArg("o","output_image","output image file",true,"","string");
    cmd.add( outputImageArg );
    TCLAP::ValueArg<std::string> referenceImageArg("r","reference_image","reference image file",false,"","string");
    cmd.add( referenceImageArg );

    TCLAP::MultiArg<float>   spacingArg("","spacing","resolution (spacing)) in mm (1 by default)",false,"float");
    cmd.add( spacingArg );
    TCLAP::MultiArg<float>   sizeArg   ("","size","size in voxels",false,"float");
    cmd.add( sizeArg );
    TCLAP::ValueArg<int>     orderArg("","order","order of the bspline interpolator (default:1)",false,1,"int");
    cmd.add( orderArg );
             
    // Parse the args.
    cmd.parse( argc, argv );

    // Get the value parsed by each arg ----------------------------------------------------
    std::string input_file       = inputImageArg.getValue();
    std::string output_file      = outputImageArg.getValue();
    std::string reference_file   = referenceImageArg.getValue();
            
    std::vector<float> size      = sizeArg.getValue();
    std::vector<float> spacing   = spacingArg.getValue();
    int order                    = orderArg.getValue();
    
    //ITK declaration ----------------------------------------------------------------------
    typedef float PixelType;
    const   unsigned int        Dimension = 3;
    
    typedef itk::Image< PixelType, Dimension >    itkImage;
    typedef itk::ImageFileReader< itkImage >      itkReader;
    typedef itk::ImageFileWriter< itkImage >      itkWriter;
    typedef itkImage::Pointer itkImagePointer;
    typedef itk::ResampleImageFilter<itkImage, itkImage>                    itkResampleFilter;
    typedef itk::IdentityTransform<double, 3>                               itkIdentityTransform;
    typedef itk::BSplineInterpolateImageFunction<itkImage, double, double>  itkBSplineInterpolator;
  
    typedef itk::ContinuousIndex<double,3>     itkContinuousIndex;
  

    std::cout<<"Reading the input image:"<<input_file<<"\n";




    itkReader::Pointer input_reader = itkReader::New();
    input_reader->SetFileName( input_file );
    input_reader->Update();
    
    itkResampleFilter::Pointer resample = itkResampleFilter::New();
    
    //parameters for interpolation (identity transform and bspline interpolator)
    itkIdentityTransform::Pointer transform = itkIdentityTransform::New();
    transform->SetIdentity();
    int interpolationOrder = order;
    itkBSplineInterpolator::Pointer bsInterpolator = itkBSplineInterpolator::New();
    bsInterpolator->SetSplineOrder(interpolationOrder);
  
    resample->SetTransform(transform);
    resample->SetInterpolator(bsInterpolator);
    
    if(reference_file != ""){
      std::cout<<"Using the following reference image:"<<reference_file<<"\n";
      itkReader::Pointer reference_reader = itkReader::New();
      reference_reader->SetFileName( reference_file );
      reference_reader->Update();
            
      resample->UseReferenceImageOn();
      resample->SetReferenceImage( reference_reader->GetOutput() );      
    }
    else{
      if( (size.size() < Dimension ) && (spacing.size() < Dimension ) ){
        std::cout<<"Please provide at least size or scaling factors for each dimension.\n";
        return 0;
      }
      
      itkImage::SizeType itkSize;
      itkImage::SpacingType itkSpacing;
      
      if( size.size() == Dimension )
        for(uint i=0; i<Dimension; i++)
          itkSize[i] = size[i];
          
      if( spacing.size() == Dimension )
        for(uint i=0; i<Dimension; i++)
          itkSpacing[i] = spacing[i];
          
      //Now, check if size or spacing is empty to automatically provide good size/spacing.
      itkImage::SizeType itkInputSize;
      itkImage::SpacingType itkInputSpacing;
      
      itkInputSize    = input_reader->GetOutput()->GetLargestPossibleRegion().GetSize();
      itkInputSpacing = input_reader->GetOutput()->GetSpacing();
      
      std::cout<<"Size of the input image: "<<itkInputSize[0]<<" "<<itkInputSize[1]<<" "<<itkInputSize[2]<<" \n";
      std::cout<<"Spacing of the input image: "<<itkInputSpacing[0]<<" "<<itkInputSpacing[1]<<" "<<itkInputSpacing[2]<<" \n";
                 
      if( size.size() < Dimension ){
        for(uint i=0; i<Dimension; i++)
          itkSize[i] = itkInputSize[i] * itkInputSpacing[i] / itkSpacing[i] ;        
      }    
      else
        for(uint i=0; i<Dimension; i++)
          itkSpacing[i] = itkInputSize[i] * itkInputSpacing[i] / itkSize[i] ;        
      
      std::cout<<"Size of the output image: "<<itkSize[0]<<" "<<itkSize[1]<<" "<<itkSize[2]<<" \n";
      std::cout<<"Spacing of the output image: "<<itkSpacing[0]<<" "<<itkSpacing[1]<<" "<<itkSpacing[2]<<" \n";
      
      resample->SetSize(itkSize);
      resample->SetOutputSpacing(itkSpacing);
    }

    resample->SetDefaultPixelValue(0.0);
   
    resample->SetInput( input_reader->GetOutput() );
    
    //Set the new origin correctly:
    itkImage::IndexType inputIndex;
    inputIndex[0] = 0;
    inputIndex[1] = 0;
    inputIndex[2] = 0;
    itkImage::PointType inputPoint;
    //This point should corresponds to the origin of the input image
    input_reader->GetOutput()->TransformIndexToPhysicalPoint(inputIndex,inputPoint);
    std::cout<<"Origin of the input image:"<<inputPoint<<"\n";
    
    itkContinuousIndex outputIndex;
    outputIndex[0] = -0.5 * input_reader->GetOutput()->GetSpacing()[0] + 0.5 * resample->GetOutputSpacing()[0]  ;
    outputIndex[1] = -0.5 * input_reader->GetOutput()->GetSpacing()[1] + 0.5 * resample->GetOutputSpacing()[1] ;
    outputIndex[2] = -0.5 * input_reader->GetOutput()->GetSpacing()[2] + 0.5 * resample->GetOutputSpacing()[2] ;
    
    itkImage::PointType outputPoint;
    input_reader->GetOutput()->TransformContinuousIndexToPhysicalPoint(outputIndex,outputPoint);
    std::cout<<"Origin of the output image:"<<outputPoint<<"\n";

    resample->SetOutputOrigin( outputPoint );
    resample->SetOutputDirection( input_reader->GetOutput()->GetDirection() );


    resample->Update();


    //Write the result image
    itkWriter::Pointer writer = itkWriter::New();  
    writer->SetFileName( output_file );
    writer->SetInput( resample->GetOutput() );
    writer->Update();
    
    } catch (TCLAP::ArgException &e)  // catch any exceptions
    { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return EXIT_SUCCESS;

}
