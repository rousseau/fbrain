/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 24/01/2013
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
#include <tclap/CmdLine.h>

/* Itk includes */
#include "itkImage.h"
#include "itkImageRegionIterator.h"



/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPatch.h"
#include "btkPatchTool.h"


int main (int argc, char* argv[])
{
  try {  

  TCLAP::CmdLine cmd("Normalize the grayscale values of one image using a reference image by histogram matching.", ' ', "1.0", true);
  TCLAP::ValueArg<std::string>    inputImageArg ("i","input_file","input anatomical image file",true,"","string", cmd);
  TCLAP::ValueArg<std::string>    maskImageArg  ("m","mask_file","mask image file (short)",true,"","string", cmd);
  TCLAP::ValueArg<std::string>    refImageArg   ("r","reference_file","reference anatomical image file",true,"","string", cmd);
  TCLAP::ValueArg<std::string>    outputImageArg("o","output_file","output anatomical image file",true,"","string", cmd);
  TCLAP::ValueArg< int >          blockArg      ("","block","0: pointwise, 1: blockwise, 2: fast blockwise (default is 1)",false,1,"int",cmd);  
  TCLAP::ValueArg< int >          hwnArg        ("","hwn","patch half size (default is 1)",false,1,"int",cmd);  
  TCLAP::ValueArg< unsigned int > binArg        ("b","bin","number of histogram bins (default is 1024)",false,1024,"unsigned int",cmd);
  TCLAP::ValueArg< unsigned int > matchArg      ("p","points","number of match points (default is 7)",false,7,"unsigned int",cmd);
  TCLAP::ValueArg< bool >         thresholdArg  ("t","threshold","threshold at mean intensity (default is true, meaning excluding voxels of the background)",false,true,"bool",cmd);


  // Parse the args.
  cmd.parse( argc, argv );
  
  std::string input_file                 = inputImageArg.getValue();
  std::string mask_file                  = maskImageArg.getValue();
  std::string ref_file                   = refImageArg.getValue();
  std::string output_file                = outputImageArg.getValue();
  
  int          block                     = blockArg.getValue();  
  int          hwn                       = hwnArg.getValue();
  unsigned int numberOfBins              = binArg.getValue();
  unsigned int numberOfMatchPoints       = matchArg.getValue();
  bool         threshold                 = thresholdArg.getValue();
  
  //ITK declaration
  const   unsigned int                               Dimension = 3;
  typedef itk::Image< short, Dimension >             ShortImageType;
  typedef ShortImageType::Pointer                    ShortImagePointer;
  typedef itk::ImageRegionIterator< ShortImageType > itkShortIterator;
  typedef itk::Image< float, Dimension >             FloatImageType;
  typedef FloatImageType::Pointer                    FloatImagePointer;
  typedef itk::ImageRegionIterator< FloatImageType > itkFloatIterator;
 


  std::cout<<"Read input images \n";
  ShortImagePointer inputImage = btk::ImageHelper<ShortImageType>::ReadImage(input_file);
  ShortImagePointer refImage   = btk::ImageHelper<ShortImageType>::ReadImage(ref_file);  
  ShortImagePointer maskImage  = btk::ImageHelper<ShortImageType>::ReadOrCreateImage(mask_file, inputImage, 1);
  
  //Create empty ouput image
  ShortImagePointer outputImage = btk::ImageHelper<ShortImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

  std::cout<<"Performing local histogram matching\n";
  itkShortIterator maskImageIt( maskImage, maskImage->GetLargestPossibleRegion());
  itkShortIterator inputImageIt( inputImage, inputImage->GetLargestPossibleRegion());    
  
  FloatImagePointer weightImage = btk::ImageHelper<ShortImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

  //compute characteristics of the input image
  ShortImageType::RegionType  region  = inputImage->GetLargestPossibleRegion();
  ShortImageType::SizeType    size    = region.GetSize();
  ShortImageType::SpacingType spacing = inputImage->GetSpacing();
  
  int x,y,z;
  if(block == 0)
  {
    std::cout<<"pointwise HM"<<std::endl;

    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)size[2]; z++)
    {
      for(y=0; y < (int)size[1]; y++)
      {
				for(x=0; x < (int)size[0]; x++)
				{
	  			ShortImageType::IndexType index;
	  			index[0] = x;
	  			index[1] = y;
	  			index[2] = z;
	  
	  			if( maskImage->GetPixel(index) > 0 )
          {
 	    			btk::PatchTool<short,float> myPatchTool;
  
	    			btk::Patch<short> inputPatch;
	    			inputPatch.Initialize(inputImage,hwn);
	    			inputPatch.ComputePatch(index, inputImage);
  
				    btk::Patch<short> refPatch;
	   			 	refPatch.Initialize(refImage,hwn);
	    			refPatch.ComputePatch(index, refImage);

	  			  btk::Patch<float> outputPatch;
	  			  myPatchTool.CreatePatchFromAnotherPatch(inputPatch, outputPatch);
	  			  
	    			//outputPatch.Initialize(outputImage,hwn);
	    			//outputPatch.ComputePatch(index, outputImage);

	    			myPatchTool.PatchIntensityNormalizationUsingMeanAndVariance(inputPatch, refPatch, outputPatch);
	      
	    			outputImage->SetPixel( index, (short)outputPatch.GetCentralValue() );
	  			}
				}
  		}
    }
      
  }    
  if(block == 1)
  {
    std::cout<<"blockwise HM"<<std::endl;    
    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)size[2]; z++)
    {
      for(y=0; y < (int)size[1]; y++)
      {
				for(x=0; x < (int)size[0]; x++)
				{
	  			ShortImageType::IndexType index;
	  			index[0] = x;
				  index[1] = y;
				  index[2] = z;
	  
				  if( maskImage->GetPixel(index) > 0 )
          {
  	    		btk::PatchTool<short,float> myPatchTool;
              
	    			btk::Patch<short> inputPatch;
	    			inputPatch.Initialize(inputImage,hwn);
	    			inputPatch.ComputePatch(index, inputImage);
  
				    btk::Patch<short> refPatch;
	   			 	refPatch.Initialize(refImage,hwn);
	    			refPatch.ComputePatch(index, refImage);

	  			  btk::Patch<float> outputPatch;
	  			  myPatchTool.CreatePatchFromAnotherPatch(inputPatch, outputPatch);
	 
	    			myPatchTool.PatchIntensityNormalizationUsingMeanAndVariance(inputPatch, refPatch, outputPatch);

            #pragma omp critical	      
				    myPatchTool.AddPatchToImage(index, outputPatch, outputImage, weightImage, 1.0);
	  			}
				}
      }
    }
  }
  if(block == 2)
  {
    std::cout<<"fast blockwise HM"<<std::endl;

    btk::Patch<short> tmpPatch;
    tmpPatch.Initialize(inputImage,hwn);
    
    ShortImageType::SizeType halfPatchSize = tmpPatch.GetHalfPatchSize();

    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)size[2]; z++)
    {
      if( z%(halfPatchSize[2]+1) == 0)
      {
				for(y=0; y < (int)size[1]; y++)
        {
	  			if( y%(halfPatchSize[1]+1) == 0)
	  			{
	    			for(x=0; x < (int)size[0]; x++)
	    			{
	      			if( x%(halfPatchSize[0]+1) == 0)
	      			{

								ShortImageType::IndexType index;
								index[0] = x;
								index[1] = y;
								index[2] = z;
	  
								if( maskImage->GetPixel(index) > 0 )
								{
              
         	    		btk::PatchTool<short,float> myPatchTool;

				    			btk::Patch<short> inputPatch;
	    						inputPatch.Initialize(inputImage,hwn);
	    						inputPatch.ComputePatch(index, inputImage);
  
							    btk::Patch<short> refPatch;
	   						 	refPatch.Initialize(refImage,hwn);
	    						refPatch.ComputePatch(index, refImage);

				  			  btk::Patch<float> outputPatch;
	  						  myPatchTool.CreatePatchFromAnotherPatch(inputPatch, outputPatch);

	    						myPatchTool.PatchIntensityNormalizationUsingMeanAndVariance(inputPatch, refPatch, outputPatch);

                  #pragma omp critical	      
								  myPatchTool.AddPatchToImage(index, outputPatch, outputImage, weightImage, 1.0);
								}
				      }
				    }
				  }
				}
      }
    }
  }
  //Normalization of the output image using the weight image
  //this steps can be improved by considering float images for all the computations (currently outputImage is in short)
  //note that blockwise techniques can introduce also some smoothing in the final estimate because of the aggregation strategy (mean of estimates)
  if(block >= 1)
  {
    itkFloatIterator weightIt( weightImage, weightImage->GetLargestPossibleRegion() );
    itkShortIterator outputIt( outputImage, outputImage->GetLargestPossibleRegion() );
    for ( outputIt.GoToBegin(), weightIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt, ++weightIt)
    {
      if( weightIt.Get() > 0 )
      {
        outputIt.Set( outputIt.Get() / weightIt.Get() );
      }
    }
  }

  
  
  
  
  btk::ImageHelper<ShortImageType>::WriteImage(outputImage, output_file);

  
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
  return 1;  
}
