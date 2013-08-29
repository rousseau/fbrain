/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

14 / 02 / 2013
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
  TCLAP::ValueArg<std::string> inputMaskArg  ("m","mask_file","filename of the mask image (short)",false,"","string", cmd);
  TCLAP::ValueArg< int >       hwnArg        ("","hwn","patch half size (default is 1)",false,1,"int",cmd);  
  TCLAP::ValueArg< int >       hwvsArg       ("","hwvs","half size of the volume search area, i.e. the spatial bandwidth (default is 5)",false,5,"int",cmd);
  TCLAP::ValueArg< float >     betaArg       ("b","beta","beta value is the smoothing parameter of NLM (default is 1)",false,1,"float",cmd);
  TCLAP::ValueArg< int >       normArg       ("n","norm","normalization with respect to the number of neighbours (default is yes (1))",false,1,"int",cmd);
  

  // Parse the args.
  cmd.parse( argc, argv );

    
  // Get the value parsed by each arg. 
  std::string input_file       = inputImageArg.getValue();
  std::string output_file      = outputImageArg.getValue();  
  std::string mask_file        = inputMaskArg.getValue();
  
  
  int          hwn             = hwnArg.getValue();
  int          hwvs            = hwvsArg.getValue();
  float        beta            = betaArg.getValue();
  int          norm            = normArg.getValue();
  
  
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
  ShortImagePointer maskImage  = btk::ImageHelper<ShortImageType>::ReadOrCreateImage(mask_file, inputImage, 1);
    
  //compute characteristics of the input image
  ShortImageType::RegionType  region  = inputImage->GetLargestPossibleRegion();
  ShortImageType::SizeType    size    = region.GetSize();
  ShortImageType::SpacingType spacing = inputImage->GetSpacing();

  //Create empty ouput image
  FloatImagePointer outputImage = btk::ImageHelper<ShortImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());
  
  //Estimation of noise level
  btk::Noise<short> myNoiseTool;
  myNoiseTool.SetSigma2Image( btk::ImageHelper<ShortImageType, FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer()) );
  myNoiseTool.ComputeGlobalSigma(inputImage, maskImage);
  
  
  //Gerer des images de "reference" en vecteur
  //Etre capable de concatener des vecteurs de voisins ou calcul sur un vecteur d'images
  //regarder l'evolution du poids d'un point en fonction du nombre d'images utilisees
  //Attention, cela dépend aussi de la méthode de selection des patches
  //-> tester les 3 strategies (rien, chi2, mean/var)
  //exp sur brainweb testera l'aspect geometrique des patches
  
  int x,y,z;
	
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
  			//Create the current patch----------------------------------------------------
	    	btk::Patch<short> inputPatch;
	    	inputPatch.Initialize(inputImage,hwn);
	    	inputPatch.ComputePatch(index, inputImage);
  
  			//Get the global/local smoothing parameter -----------------------------------
  			ShortImageType::SizeType hps = inputPatch.GetHalfPatchSize();
          	float sigma2 = myNoiseTool.GetSigma2Image()->GetPixel(index);
  			float smoothing = 2 * beta * sigma2 * (2*hps[0]+1) * (2*hps[1]+1) * (2*hps[2]+1);
  			  
            btk::PatchTool<short,short> myPatchTool;
            myPatchTool.SetSpatialBandwidth(hwvs, inputPatch);
            myPatchTool.SetImageSize(size);
            myPatchTool.SetPatchSelectionMethod(0);
            myPatchTool.SetParamPatchSelectionMethod(0.95);
            myPatchTool.SetSelectionNeighbourPatchMethod(inputPatch);
			  
  			//Create the set of patch candidates------------------------------------------
  			//This is clearly not optimised but the code is clearer ----------------------
  			std::vector< btk::Patch<short> > neighbourPatches;
			myPatchTool.GetNeighbourPatches(inputPatch, neighbourPatches, inputImage);
			
			//Compute weights of the set of patches---------------------------------------			
			std::vector< float > neighbourWeights;
			myPatchTool.ComputeNeighbourWeights(inputPatch, neighbourPatches, neighbourWeights, smoothing);

			//Need to take into account the central weight--------------------------------
			  			  			
			
			//Compute sum of weights
			float sum = 0;
			for(unsigned int i=0; i < neighbourWeights.size(); i++)
				sum += neighbourWeights[i];
				
			if(norm==1)	
		      sum /= neighbourWeights.size();
			
  			outputImage->SetPixel(index, sum);
  			
  
	  	  }
		}
  	  }
    }

  
  btk::ImageHelper<FloatImageType>::WriteImage(outputImage, output_file);

  
  return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
