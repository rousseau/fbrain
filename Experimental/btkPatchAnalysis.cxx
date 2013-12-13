/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

08 october 2013
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
#include "numeric"      // std::accumulate

/* Itk includes */
#include "itkImage.h"
#include "itkNormalVariateGenerator.h"
#include "itkImageRegionIterator.h"

#include "itkDivideImageFilter.h"

#include "vnl/vnl_random.h"

/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPatch.h"
#include "btkPatchTool.h"
#include "btkNoise.h"

#include "btkPatch2.h"
#include "btkPatchTool2.h"

int main(int argc, char** argv)
{
  try {
  
  	TCLAP::CmdLine cmd("Patch-based image analysis", ' ', "0.1", true);

  	TCLAP::ValueArg<std::string> inputImageArg ("i","image_file","input image file (short)",true,"","string", cmd);
  	TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (float)",true,"","string", cmd);
  	TCLAP::ValueArg<std::string> inputMaskArg  ("m","mask_file","filename of the mask image (short)",false,"","string", cmd);
	TCLAP::ValueArg< int >       xArg          ("x","x","x coordinate (default is 10)",false,10,"int",cmd);
	TCLAP::ValueArg< int >       yArg          ("y","y","y coordinate (default is 10)",false,10,"int",cmd);
	TCLAP::ValueArg< int >       zArg          ("z","z","z coordinate (default is 10)",false,10,"int",cmd);
  	TCLAP::ValueArg< int >       hwnArg        ("","hwn","patch half size",false,1,"int",cmd);
  	TCLAP::ValueArg< int >       hwvsArg       ("","hwvs","half size of the volume search area, i.e. the spatial bandwidth",false,5,"int",cmd);
  	TCLAP::ValueArg< float >     betaArg       ("b","beta","beta value is the smoothing parameter of NLM (default is 1)",false,1,"float",cmd);
    TCLAP::ValueArg< int >       patchSelectArg("p","patchselect","Patch selection method: 0: chi2, 1..., 2: mean and variance (default is 0)",false,0,"int",cmd);

  	// Parse the args.
  	cmd.parse( argc, argv );

  	// Get the value parsed by each arg. 
  	std::string input_file       = inputImageArg.getValue();
  	std::string output_file      = outputImageArg.getValue(); 
  	std::string mask_file        = inputMaskArg.getValue();
  	  
  	int x = xArg.getValue();
  	int y = yArg.getValue();
  	int z = zArg.getValue();
  	
  	int   hwn             = hwnArg.getValue();
  	int   hwvs            = hwvsArg.getValue();
	float beta            = betaArg.getValue();
  	int   patchSelect     = patchSelectArg.getValue();


  	//ITK declaration
  	const   unsigned int                               Dimension = 3;
  
 	typedef itk::Image< short, Dimension >             ShortImageType;
  	typedef ShortImageType::Pointer                    ShortImagePointer;

    typedef itk::Image< float, Dimension >             FloatImageType;
  	typedef FloatImageType::Pointer                    FloatImagePointer;
    typedef itk::ImageRegionIterator< FloatImageType > FloatIterator;

  	std::cout<<"Read input image \n";
    ShortImagePointer inputImage = btk::ImageHelper<ShortImageType>::ReadImage(input_file);
    ShortImagePointer maskImage  = btk::ImageHelper<ShortImageType>::ReadOrCreateImage(mask_file, inputImage, 1);

  	std::cout<<"Create empty output image\n";
    FloatImagePointer outputImage = btk::ImageHelper<ShortImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());
    FloatImagePointer weightImage = btk::ImageHelper<ShortImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

  	//compute characteristics of the input image
	ShortImageType::RegionType  region  = inputImage->GetLargestPossibleRegion();
	ShortImageType::SizeType    size    = region.GetSize();
  	ShortImageType::SpacingType spacing = inputImage->GetSpacing();

    //--------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------

    //Estimation of noise level
    btk::Noise<short> myNoiseTool;
    myNoiseTool.SetSigma2Image( btk::ImageHelper<ShortImageType, FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer()) );
    myNoiseTool.ComputeGlobalSigma(inputImage, maskImage);

	ShortImageType::IndexType index;
	index[0] = x;
  	index[1] = y;
	index[2] = z;

    ShortImageType::SizeType    patchSize;
    patchSize[0] = 3;
    patchSize[1] = 3;
    patchSize[2] = 3;

    btk::PatchTool2<short> myPatchTool2;
    myPatchTool2.ComputePatchSize(inputImage, hwn);
    myPatchTool2.ComputeSpatialBandwidth(inputImage, hwvs);

    double smoothing = 2 * beta * myNoiseTool.GetGlobalSigma2() * myPatchTool2.GetFullPatchSize()[0] * myPatchTool2.GetFullPatchSize()[1] * myPatchTool2.GetFullPatchSize()[2];
    std::cout<<"Smoothing = "<<smoothing<<", estimated variance = "<<myNoiseTool.GetGlobalSigma2()<<std::endl;

    FloatImagePointer meanImage     = btk::ImageHelper<ShortImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());
    FloatImagePointer varianceImage = btk::ImageHelper<ShortImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

    myPatchTool2.ComputeMeanAndVarianceImage(inputImage, maskImage, meanImage, varianceImage);


    /*
	if( maskImage->GetPixel(index) > 0 )
    {
        btk::Patch2<short> inputPatch2 = btk::Patch2<short>(inputImage, index, patchSize);

        //Neighbours selection (using mean and variance ... or nothing)
        std::vector< btk::Patch2<short> > neighbourPatchesTmp2;
        //myPatchTool2.GetNeighbourPatches(index, inputImage, neighbourPatchesTmp2);
        myPatchTool2.GetNeighbourPatchesUsingMeanAndVariance(index, inputImage, meanImage->GetPixel(index), varianceImage->GetPixel(index), meanImage, varianceImage, neighbourPatchesTmp2);

        //Post selection using chi2 or nearest neighbours
        //TO BE IMPLEMENTED

        std::cout<<"Nb voisins : "<<neighbourPatchesTmp2.size()<<std::endl;

        //Compute weights
        std::vector< float > neighbourWeightsTmp2(neighbourPatchesTmp2.size(),0);
        for(unsigned int n = 0; n < neighbourPatchesTmp2.size(); n++)
          neighbourWeightsTmp2[n] = exp( - myPatchTool2.ComputeL2NormBetweenPatches(inputPatch2, neighbourPatchesTmp2[n]) / smoothing);

        //weight normalization
        double sumOfWeights = std::accumulate( neighbourWeightsTmp2.begin(), neighbourWeightsTmp2.end(), 0 );
        if( sumOfWeights > 0 )
            for(unsigned int n = 0; n < neighbourWeightsTmp2.size(); n++)
                neighbourWeightsTmp2[n] = neighbourWeightsTmp2 / sumOfWeights;

        for(unsigned int n = 0; n < neighbourPatchesTmp2.size(); n++)
            outputImage->SetPixel(neighbourPatchesTmp2[n].GetCentralPointInImage(), 100*neighbourWeightsTmp2[n]);

    }
    */
/*
    typename ShortImageType::RegionType::IndexType centralPoint;
    centralPoint[0] = myPatchTool2.GetHalfPatchSize()[0];
    centralPoint[1] = myPatchTool2.GetHalfPatchSize()[1];
    centralPoint[2] = myPatchTool2.GetHalfPatchSize()[2];
*/

    int maxNeighbours = myPatchTool2.GetFullPatchSize()[0] * myPatchTool2.GetFullPatchSize()[1] * myPatchTool2.GetFullPatchSize()[2];

    #pragma omp parallel for private(x,y,z) schedule(dynamic)
    for(z=0; z < (int)size[2]; z++)
    for(y=0; y < (int)size[1]; y++)
    for(x=0; x < (int)size[0]; x++)
    {
      typename ShortImageType::IndexType p;
      p[0] = x;
      p[1] = y;
      p[2] = z;
      if( maskImage->GetPixel(p) > 0 )
      {

        //Get neighbour coordinates of the current voxel
        std::vector< typename ShortImageType::IndexType > neighbours;
        neighbours.reserve(maxNeighbours); //reserve sufficient space for vector building to avoid unnecessary resizing
        myPatchTool2.GetNeighboursUsingMeanAndVariance(p, inputImage, meanImage->GetPixel(p), varianceImage->GetPixel(p), meanImage, varianceImage, neighbours);
        //myPatchTool2.GetNeighbours(p, inputImage, neighbours);
        //myPatchTool2.RemoveCentralPointInNeighborhood(p, neighbours);

        double sumOfWeights = 0.0;
        double maxWeight = 0.0;


        //blockwise, mean aggregation
        btk::Patch2<float> outputPatch = btk::Patch2<float>(outputImage, patchSize);
        myPatchTool2.ComputeWeightedMeanOfPatches(inputImage, p, neighbours, outputPatch, sumOfWeights, maxWeight, smoothing);

        FloatIterator outputIt( outputPatch.GetData(), outputPatch.GetData()->GetLargestPossibleRegion());

        if( sumOfWeights > 0)
          for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
          {
            outputIt.Set( outputIt.Get() / sumOfWeights );
          }

        sumOfWeights = 1.0;
        #pragma omp critical
        myPatchTool2.AddPatchToImage(p, outputPatch, outputImage, weightImage, sumOfWeights);
        //---------------------------------------------------------------------------------------------------------------

        /*
        //pointwise
        double outputValue = 0;
        myPatchTool2.ComputeWeightedMeanAtPatchCenter(inputImage, p, neighbours, outputValue, sumOfWeights, maxWeight, smoothing);
        if( sumOfWeights > 0 )
          outputImage->SetPixel(p, outputValue/sumOfWeights );
        */
      }
    }

    FloatIterator it1( outputImage, outputImage->GetLargestPossibleRegion());
    FloatIterator it2( weightImage, weightImage->GetLargestPossibleRegion());
    for ( it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2)
        if(it2.Get() > 0)
            it1.Set( it1.Get() / it2.Get() );

    //Sinon remettre la valeur input



    btk::ImageHelper<FloatImageType>::WriteImage(outputImage, output_file);
    btk::ImageHelper<FloatImageType>::WriteImage(weightImage, "weightImage.nii.gz");

    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
