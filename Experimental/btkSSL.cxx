/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

12 / 02 / 2013
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


/* Standard includes */
#include <tclap/CmdLine.h>
#include "vector"
#include "sstream"

/* Itk includes */
#include "itkImage.h"
#include "itkDiscreteGaussianImageFilter.h"


/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPatch.h"
#include "btkPatchTool.h"
#include "btkNoise.h"

int main(int argc, char** argv)
{
  try {  

  TCLAP::CmdLine cmd("Semi-supervised Learning", ' ', "1.0", true);

  TCLAP::ValueArg<std::string> inputImageArg ("i","image_file","input image file (short)",true,"","string", cmd);
  TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string", cmd);  
  TCLAP::MultiArg<std::string> probaImageArg ("p","proba_file","probabilistic images of the textbook (float) (multiple inputs required) ",true,"string", cmd);  
  TCLAP::MultiArg<std::string> outputprobaImageArg ("w","output_proba_file","output probabilistic images (multiple inputs required) ",false,"string", cmd);  
  TCLAP::ValueArg<std::string> inputMaskArg  ("m","mask_file","filename of the mask image (short)",false,"","string", cmd);
  TCLAP::MultiArg< int >       hwnArg        ("","hwn","patch half size",true,"int",cmd);  
  TCLAP::MultiArg< int >       hwvsArg       ("","hwvs","half size of the volume search area, i.e. the spatial bandwidth",true,"int",cmd);
  TCLAP::ValueArg< float >     alphaArg      ("a","alpha","alpha value is the tradeoff between prior and graph-based consistency (default is 0.5)",false,0.5,"float",cmd);
  TCLAP::ValueArg< float >     betaArg       ("b","beta","beta value is the smoothing parameter of NLM (default is 1)",false,1,"float",cmd);
  TCLAP::ValueArg< float >     gammaArg      ("g","gamma","gamma value is the tradeoff between static and adaptive version (default is 1)",false,1,"float",cmd);
  TCLAP::ValueArg< float >     epsilonArg    ("e","epsilon","epsilon value is the minimum % of changes between 2 iterations (default is 0.001)",false,0.001,"float",cmd);
  TCLAP::ValueArg< int >       itermaxArg    ("","itermax","maximum of iterations (default is 10)",false,10,"int",cmd);
  TCLAP::ValueArg< int >       patchSelectArg("","patchselect","Patch selection method: 0: chi2, 1..., 2: mean and variance (default is 0)",false,0,"int",cmd);
  TCLAP::ValueArg< float >     fastArg       ("","fast","threshold for speeding up the labeling (default is 0.95)",false,0.95,"float",cmd);
  TCLAP::ValueArg< float >     varArg        ("","var","variance of the gaussian filter (default is 2.5)",false,2.5,"float",cmd);
  

  // Parse the args.
  cmd.parse( argc, argv );

    
  // Get the value parsed by each arg. 
  std::string input_file       = inputImageArg.getValue();
  std::string output_file      = outputImageArg.getValue();  
  std::vector<std::string> proba_file = probaImageArg.getValue();
  std::vector<std::string> output_proba_file = outputprobaImageArg.getValue();
  std::string mask_file        = inputMaskArg.getValue();
  
  
  std::vector<int> hwn             = hwnArg.getValue();
  std::vector<int> hwvs            = hwvsArg.getValue();
  float            alpha           = alphaArg.getValue();
  float            beta            = betaArg.getValue();
  float            gamma           = gammaArg.getValue();
  float            epsilon         = epsilonArg.getValue();
  int              iterMax         = itermaxArg.getValue();
  int              patchSelect     = patchSelectArg.getValue();
  float            probaThreshold  = fastArg.getValue();
  float            var             = varArg.getValue();
  
  if( hwn.size() != hwvs.size() )
  {
    std::cout<<"Warning : vector sizes are different (hwn and hwvs)\n";
    return 0;
  }
  std::cout<<"Threshold on probability maps : "<<probaThreshold<<std::endl;
  
  //ITK declaration
  const   unsigned int                               Dimension = 3;
  
  typedef itk::Image< short, Dimension >             ShortImageType;
  typedef ShortImageType::Pointer                    ShortImagePointer;
  typedef itk::ImageRegionIterator< ShortImageType > itkShortIterator;

  
  typedef itk::Image< float, Dimension >             FloatImageType;
  typedef FloatImageType::Pointer                    FloatImagePointer;
  typedef itk::ImageRegionIterator< FloatImageType > itkFloatIterator;

  typedef itk::DiscreteGaussianImageFilter< FloatImageType,FloatImageType > itkGaussianFilter;

  std::cout<<"Read input images \n";
  ShortImagePointer inputImage = btk::ImageHelper<ShortImageType>::ReadImage(input_file);
  ShortImagePointer maskImage  = btk::ImageHelper<ShortImageType>::ReadOrCreateImage(mask_file, inputImage, 1);
  
  std::vector< FloatImagePointer > probaImages = btk::ImageHelper<FloatImageType>::ReadImage(proba_file);
  std::vector< FloatImagePointer > probaImagesStepK = btk::ImageHelper<FloatImageType>::ReadImage(proba_file);
  std::vector< FloatImagePointer > probaImagesStepKPlusOne = btk::ImageHelper<FloatImageType>::ReadImage(proba_file);

  unsigned int numberOfClasses = probaImages.size();
  
  //compute characteristics of the input image
  ShortImageType::RegionType  region  = inputImage->GetLargestPossibleRegion();
  ShortImageType::SizeType    size    = region.GetSize();
  ShortImageType::SpacingType spacing = inputImage->GetSpacing();

  //Create empty ouput image
  ShortImagePointer outputImage = btk::ImageHelper<ShortImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

  int block = 0;
  float stoppingCriterion = epsilon; //percentage of changes between 2 iterations
  int numberOfChanges = size[0]*size[1]*size[2];
  int numberOfPoints = size[0]*size[1]*size[2]; //we should consider only the points of the mask
  int maxLoop = iterMax;
  int currentLoop = 1;
  
  
  //Add a step for normalizing weights if necessary
  
  //Estimation of noise level
  btk::Noise<short> myNoiseTool;
  myNoiseTool.SetSigma2Image( btk::ImageHelper<ShortImageType, FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer()) );
  myNoiseTool.ComputeGlobalSigma(inputImage, maskImage);
  
  int x,y,z;
  
  std::cout<<"Initializing the output image\n";  
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
	  	  double maxProba = 0;
		  int maxLabel = 0;
	  	  for(unsigned int l = 0; l < numberOfClasses; l++)
	  	    if( maxProba < probaImagesStepK[l]->GetPixel( index ) )
	  	    {
	  	      maxProba = probaImagesStepK[l]->GetPixel( index );
	  	      maxLabel = l+1;
	  	    }

		  outputImage->SetPixel(index, maxLabel);
        }
      }
    }
  }

  std::cout<<"Normalizing input probabilities\n";  
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
	  	  double sumProba = 0;
	  	  for(unsigned int l = 0; l < numberOfClasses; l++)
	  	    sumProba += probaImages[l]->GetPixel( index );

          if(sumProba > 0)
	  	    for(unsigned int l = 0; l < numberOfClasses; l++)
	  	    {
		      probaImages[l]->SetPixel(index, probaImages[l]->GetPixel(index) / sumProba);    
		      probaImagesStepK[l]->SetPixel(index, probaImagesStepK[l]->GetPixel(index) / sumProba);    
		      probaImagesStepKPlusOne[l]->SetPixel(index, probaImagesStepKPlusOne[l]->GetPixel(index) / sumProba);    
		    }
		}
      }
    }
  }   
  
  while( (currentLoop < maxLoop) && ( numberOfChanges > numberOfPoints * stoppingCriterion) )
  {

  if(block == 0)
  {
    std::cout<<"pointwise SSL"<<std::endl;
    
    numberOfChanges = 0;
    unsigned int nbp = 0;
    unsigned int nbpnull = 0;
	
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
	  	  
	  	  double maxProba = 0;
	  	  double sumProba = 0;
	  	  for(unsigned int l = 0; l < numberOfClasses; l++)
	  	  {
	  	    if( maxProba < probaImagesStepK[l]->GetPixel( index ) )
	  	    {
	  	      maxProba = probaImagesStepK[l]->GetPixel( index );	  	       
	  	    }
	  	    sumProba += probaImagesStepK[l]->GetPixel( index );
	  	  }
	  	  if(sumProba > 0)
	  	    maxProba /= sumProba;
	      	      
	  	  if( (maskImage->GetPixel(index) > 0) && (maxProba < probaThreshold) )
	  	  //if( maskImage->GetPixel(index) > 0 )
          {
            short currentLabel = outputImage->GetPixel( index ) ;
            
          	std::vector< double > newWeights(numberOfClasses);
          	
  			std::vector< btk::Patch<short> > neighbourPatches;
			std::vector< float > neighbourWeights;
          	
            for(unsigned int h = 0; h < hwn.size(); h++)
            {
              //Create the current patch----------------------------------------------------
	    	  btk::Patch<short> inputPatch;
	    	  inputPatch.Initialize(inputImage,hwn[h]);
	    	  inputPatch.ComputePatch(index, inputImage);

  	    	  //Get the global/local smoothing parameter -----------------------------------
  			  ShortImageType::SizeType hps = inputPatch.GetHalfPatchSize();
          	  float sigma2 = myNoiseTool.GetSigma2Image()->GetPixel(index);
  			  float smoothing = 2 * beta * sigma2 * (2*hps[0]+1) * (2*hps[1]+1) * (2*hps[2]+1);
  			
  			  //std::cout<<"Noise : "<<beta<<" "<<sigma2<<" "<<smoothing<<std::endl;
  
              btk::PatchTool<short,short> myPatchTool;
              myPatchTool.SetSpatialBandwidth(hwvs[h], inputPatch);
              myPatchTool.SetImageSize(size);
              myPatchTool.SetPatchSelectionMethod(patchSelect);
              myPatchTool.SetParamPatchSelectionMethod(0.95);
              myPatchTool.SetSelectionNeighbourPatchMethod(inputPatch);

              //Create the set of patch candidates------------------------------------------
  			  //This is clearly not optimised but the code is clearer ----------------------
  			  std::vector< btk::Patch<short> > neighbourPatchesTmp;
			  myPatchTool.GetNeighbourPatches(inputPatch, neighbourPatchesTmp, inputImage);

			  if(myPatchTool.GetPatchSelectionMethod()==2)
			  {
				inputPatch.ComputeMeanAndStdDevValues();
				for(unsigned int n = 0; n < neighbourPatchesTmp.size(); n++)
					neighbourPatchesTmp[n].ComputeMeanAndStdDevValues();
			  }

              //Compute weights of the set of patches---------------------------------------	
			  std::vector< float > neighbourWeightsTmp;		
			  myPatchTool.ComputeNeighbourWeights(inputPatch, neighbourPatchesTmp, neighbourWeightsTmp, smoothing);
			  
			  
			  //Concatenation of neighbour patches and weights
			  neighbourPatches.insert( neighbourPatches.end(), neighbourPatchesTmp.begin(), neighbourPatchesTmp.end() );
			  neighbourWeights.insert( neighbourWeights.end(), neighbourWeightsTmp.begin(), neighbourWeightsTmp.end() );
  
            }          	
  			
			//Need to take into account the central weight--------------------------------
			
			//Update labels --------------------------------------------------------------
  			
  			double maxWeight = 0;
  			int   maxLabel  = 0;
  			
  			double sumOfProbaImagesStepKPlusOne = 0;
  			
  			for(unsigned int l = 0; l < numberOfClasses; l++)
  			{
  				double sumOfWeights = 0;
  				double sumOfConsistencyWeights = 0;
  				newWeights[l] = 0;
  				  				
  				double localConsistency = 0;
  				  				
  				for(unsigned int n = 0; n < neighbourPatches.size(); n++)
  				{
  				    //float currentWeight = (1-alpha)*neighbourWeights[n];
  				    //sumOfWeights += currentWeight;
  				    
  				    //newWeights[l] += currentWeight * probaImagesStepK[l]->GetPixel( neighbourPatches[n].GetCentralPointInImage() );
  				    
  				    localConsistency += neighbourWeights[n] * probaImagesStepK[l]->GetPixel( neighbourPatches[n].GetCentralPointInImage() );
  				    sumOfConsistencyWeights += neighbourWeights[n];
  				    
  				    //Weighting using patch size
  				    //double weightedWeight = neighbourWeights[n] * neighbourPatches[n].GetFullPatchSize()[0] * neighbourPatches[n].GetFullPatchSize()[1];
  				    //localConsistency +=  weightedWeight * probaImagesStepK[l]->GetPixel( neighbourPatches[n].GetCentralPointInImage() );
  				    //sumOfConsistencyWeights += weightedWeight;  				    
  				    
  				}	
  				
  				if( sumOfConsistencyWeights > 0)
  				  localConsistency = localConsistency/sumOfConsistencyWeights;
  				
  				sumOfWeights = alpha + (1-alpha)*localConsistency;
  				//sumOfWeights += alpha;
  				
			    //Static version
  				//newWeights[l] += alpha * probaImages[l]->GetPixel( inputPatch.GetCentralPointInImage() ); 
  				//Adaptive version
  				//newWeights[l] += alpha * probaImagesStepK[l]->GetPixel( inputPatch.GetCentralPointInImage() ); 
  				
  				double newPrior = gamma * probaImages[l]->GetPixel( index ) + (1.0-gamma) * probaImagesStepK[l]->GetPixel( index );
  				
  				//newWeights[l] += alpha * newPrior;
  				newWeights[l] = alpha * newPrior + (1-alpha)*localConsistency;
  				
  				//probaImagesStepK(PlusOne) is the estimated proba at the current step. It is used for (possibly) updating the priors.
  			    probaImagesStepKPlusOne[l]->SetPixel(index, newWeights[l]);			  			      
  			    sumOfProbaImagesStepKPlusOne += newWeights[l];  
  			      
  			    if(newWeights[l] > maxWeight ) 
  			    {
  			      maxWeight = newWeights[l];
  			      maxLabel  = l+1;
  			    }
  			    
  			    //if( (x==10) && (y==7) && (z==12) )
  			      //std::cout<<"label : "<<l+1<<" w: "<<newWeights[l]<<" prior : "<<newPrior<<" consistency : "<<localConsistency<<" "<<sumOfConsistencyWeights<<", number of neighbours : "<<neighbourWeights.size()<<std::endl;
  			    	
  			}
  			
            //if( (x==10) && (y==7) && (z==12) )
  			  //std::cout<<"labelisation courante : "<<maxLabel<<" "<<maxWeight<<", intensity : "<<inputImage->GetPixel(index)<<std::endl;
  			  			
  			//outputImage is used as the current estimate of the final (crisp) labeling
  			#pragma omp critical
  			if(outputImage->GetPixel(index) != maxLabel)
  				numberOfChanges++;
  			  			
  			outputImage->SetPixel(index, maxLabel);
  
	  	  }
		}
  	  }
    }

  }//end of if block    

  //Copy probaImagesStepKPlusOne to probaImagesStepK
  for(unsigned int l = 0; l < numberOfClasses; l++)
  {
    probaImagesStepK[l] = btk::ImageHelper<FloatImageType>::DeepCopy( probaImagesStepKPlusOne[l].GetPointer() );
  }

  //Then apply Gaussian filter
  if(var > 0)
  {
    std::cout<<"Smoothing estimated label weights\n";
    for(unsigned int l = 0; l < numberOfClasses; l++)
    {
      itkGaussianFilter::Pointer filter = itkGaussianFilter::New();
      filter->SetInput(probaImagesStepK[l]);

      itkGaussianFilter::ArrayType variance;
      variance[0] = var;
      variance[1] = var;
      variance[2] = var;
      filter->SetVariance(variance);

      filter->SetUseImageSpacingOff();
      filter->Update();

      probaImagesStepK[l] = filter->GetOutput();
    }
  }
  std::cout<<"Normalizing probabilities ...\n";  
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
	  	  double sumProba = 0;
	  	  for(unsigned int l = 0; l < numberOfClasses; l++)
	  	    sumProba += probaImagesStepK[l]->GetPixel( index );

          if(sumProba > 0)
	  	    for(unsigned int l = 0; l < numberOfClasses; l++)
		      probaImagesStepK[l]->SetPixel(index, probaImagesStepK[l]->GetPixel(index) / sumProba);    
		}     
      }
    }
  }  
  
 //Compute the percentage of changed labels (convergence criterion)

  
    std::cout<<"iteration : "<<currentLoop<<std::endl;
    std::cout<<"number of changes : "<<numberOfChanges<<std::endl;
    std::cout<<"stopping criterion: "<<numberOfPoints * stoppingCriterion<<std::endl;
    
	currentLoop ++;
	
  }//end of while







  
  
  btk::ImageHelper<ShortImageType>::WriteImage(outputImage, output_file);
  
  if(output_proba_file.size() > 0)
    btk::ImageHelper<FloatImageType>::WriteImage(probaImagesStepK, output_proba_file);
  	
  
  return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
