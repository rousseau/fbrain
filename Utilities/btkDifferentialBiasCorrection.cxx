/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

17 february 2012
rousseau@unistra.fr

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
*/

/*
This program implements a differential bias correction method :
K. Leung, G. Ridgway, S. Ourselin, N. Fox, and ADNI
Consistent-multi-time-point brain atrophy estimation from the boundary shift integral
Neuroimage 59 (2012) 3995-4005
*/



#include <tclap/CmdLine.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"
#include "itkLog10ImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkMultiplyImageFilter.h"
//#include "itkPowImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkDivideImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include <vector>

int main(int argc, char** argv)
{

  try {

    TCLAP::CmdLine cmd("Differential bias correction method for n images (Leung et al. Neuroimage 2012) ", ' ', "1.0", true);

    TCLAP::MultiArg<std::string> inputArg  ("i","input","input image files.", true,"string",cmd);
    TCLAP::MultiArg<std::string> outputArg ("o","output","output image files.", true,"string",cmd);
    TCLAP::MultiArg<std::string> biasArg   ("b","bias","estimated differential bias.", false,"string",cmd);
    TCLAP::ValueArg<int>         radiusArg ("r","radius","radius of the median filter. Default=5.", false,5,"int",cmd);
   
 
    // Parse the args.
    cmd.parse( argc, argv );


    // Get the value parsed by each arg.
    std::vector<std::string> input_file       = inputArg.getValue();
    std::vector<std::string> output_file      = outputArg.getValue();
    std::vector<std::string> bias_file        = biasArg.getValue();
    int radius                                = radiusArg.getValue();
    

    //ITK declaration
    typedef float PixelType;
    typedef itk::Image< PixelType, 3>          itkImage;
    typedef itkImage::Pointer                  itkPointer;  
    typedef itk::ImageFileReader< itkImage >   itkReader;
    typedef itk::ImageFileWriter< itkImage >   itkWriter;
    typedef itk::ImageDuplicator< itkImage >   itkDuplicator;
    typedef itk::Log10ImageFilter<itkImage, itkImage>    itkLog10Filter;
    typedef itk::MedianImageFilter<itkImage, itkImage>   itkMedianFilter;
    typedef itk::SubtractImageFilter<itkImage>           itkSubtractFilter;
    typedef itk::ExpImageFilter<itkImage, itkImage>      itkExpFilter;
    typedef itk::MultiplyImageFilter<itkImage, itkImage> itkMultiplyFilter;
    //typedef itk::PowImageFilter<itkImage, itkImage>      itkPowFilter;
    typedef itk::ImageRegionIterator< itkImage >                 itkIterator;
    typedef itk::DivideImageFilter<itkImage, itkImage, itkImage> itkDivideFilter;
    typedef itk::StatisticsImageFilter<itkImage>                 itkStatisticsImageFilter;


    itkMedianFilter::InputSizeType itkRadius;
    itkRadius.Fill(radius);
    
    
    std::vector<itkPointer>     inputImages;
    std::vector<itkPointer>     log10Images;
    std::vector<itkPointer>     outputImages;
    std::vector<itkPointer>     biasImages;
    
    unsigned int numberOfImages = input_file.size();
    inputImages.resize(numberOfImages);
    log10Images.resize(numberOfImages);
    outputImages.resize(numberOfImages);
    biasImages.resize(numberOfImages);
    
    //Read input images  
    for(unsigned int i=0;i<numberOfImages;i++){
      std::cout<<"Reading Input Image : "<<input_file[i]<<"\n";
      itkReader::Pointer reader = itkReader::New();
      reader->SetFileName( input_file[i]  );
      reader->Update();
      inputImages[i] = reader->GetOutput();
    }
    
    std::cout<<"Compute log10 images\n";
    for(unsigned int i=0;i<numberOfImages;i++){
      itkLog10Filter::Pointer log10 = itkLog10Filter::New();
      log10->SetInput(inputImages[i]);
      log10->Update();
      log10Images[i] = log10->GetOutput();
    }
      
    //Each ratio image is computed. (This is not optimal at all. Best option should be to compute the upper matrix of ratio images.)
    for(unsigned int i=0;i<numberOfImages;i++){
      std::cout<<"Computing correction for image "<<i+1<<"\n";
      std::vector<itkPointer>     ratioImages;
      ratioImages.resize(numberOfImages-1);
      // Rij = exp(median(log(Ii)-log(Ij)))
            
      int current = 0;
      for(unsigned int j=0;j<numberOfImages;j++)
        if(j!=i){
      
          std::cout<<"Compute ratio between "<<i+1<<" and "<<j+1<<"\n";     
          itkSubtractFilter::Pointer diff = itkSubtractFilter::New();
          diff->SetInput1(log10Images[i]);
          diff->SetInput2(log10Images[j]);
          diff->Update();
          
          itkMedianFilter::Pointer median = itkMedianFilter::New();
          median->SetRadius(itkRadius);
          median->SetInput( diff->GetOutput() );
          median->Update();
          
          itkExpFilter::Pointer exp = itkExpFilter::New();
          exp->SetInput(median->GetOutput());
          exp->Update();
          ratioImages[current] = exp->GetOutput();
    
          current++;
        }
      
      //duplicate the input image into the bias image to keep all header information
      itkDuplicator::Pointer duplicator = itkDuplicator::New();
      duplicator->SetInputImage( inputImages[i] );
      duplicator->Update();
      biasImages[i] = duplicator->GetOutput();
      biasImages[i]->FillBuffer(1);
      
      std::cout<<"Compute the differential bias\n";
      for(unsigned int k=0;k<numberOfImages-1;k++){
        itkMultiplyFilter::Pointer multiply = itkMultiplyFilter::New();
        multiply->SetInput1(biasImages[i]);
        multiply->SetInput2(ratioImages[k]);
        multiply->Update();
        biasImages[i] = multiply->GetOutput();
  
      }
    
      //Only available in ITK4.0
      /*
      itkPowFilter::Pointer pow = itkPowFilter::New();
      pow->SetInput1(biasImages[i]);
      pow->SetConstant2( 1/(numberOfImages-1) );
      pow->Update();
      biasImages[i] = pow->GetOutput();
      */
      PixelType p = 1.0/(numberOfImages-1.0);
      
      itkIterator itBias(biasImages[i],biasImages[i]->GetLargestPossibleRegion());
      for(itBias.GoToBegin(); !itBias.IsAtEnd(); ++itBias)
        itBias.Set( pow(itBias.Get(),p) );
      
        
      std::cout<<"Compute the output image\n";
      itkDivideFilter::Pointer divide = itkDivideFilter::New();
      divide->SetInput1(inputImages[i]);
      divide->SetInput2(biasImages[i]);
      divide->Update();
      outputImages[i] = divide->GetOutput();        
    }  

    //Write output images
    for(unsigned int i=0;i<numberOfImages;i++){
      std::cout<<"Writing output Image : "<<output_file[i]<<"\n";
      itkWriter::Pointer writer = itkWriter::New();
      writer->SetFileName( output_file[i]  );
      writer->SetInput( outputImages[i] );
      writer->Update();
    }

    //Write bias images
    if(bias_file.size() > 0){
      for(unsigned int i=0;i<numberOfImages;i++){
        std::cout<<"Writing bias Image : "<<bias_file[i]<<"\n";
        itkWriter::Pointer writer = itkWriter::New();
        writer->SetFileName( bias_file[i]  );
        writer->SetInput( biasImages[i] );
        writer->Update();
      }
    }
    return 1;

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
