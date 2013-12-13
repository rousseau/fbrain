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
#include "btkPandoraBoxImageFilters.h"

int main(int argc, char** argv)
{
  try {
  
    TCLAP::CmdLine cmd("SSL2", ' ', "0.1", true);

  	TCLAP::ValueArg<std::string> inputImageArg ("i","image_file","input image file (short)",true,"","string", cmd);
  	TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (float)",true,"","string", cmd);
  	TCLAP::ValueArg<std::string> inputMaskArg  ("m","mask_file","filename of the mask image (short)",false,"","string", cmd);
    TCLAP::MultiArg<std::string> probaImageArg ("p","proba_file","probabilistic images of the textbook (float) (multiple inputs required) ",true,"string", cmd);
    TCLAP::MultiArg<std::string> outputprobaImageArg ("w","output_proba_file","output probabilistic images (multiple inputs required) ",false,"string", cmd);

  	// Parse the args.
  	cmd.parse( argc, argv );

  	// Get the value parsed by each arg. 
    std::string input_file                     = inputImageArg.getValue();
    std::string output_file                    = outputImageArg.getValue();
    std::string mask_file                      = inputMaskArg.getValue();
    std::vector<std::string> proba_file        = probaImageArg.getValue();
    std::vector<std::string> output_proba_file = outputprobaImageArg.getValue();

  	
    int   hwn             = 1;
    int   hwvs            = 1;
    float beta            = 1;


  	//ITK declaration
  	const   unsigned int                               Dimension = 3;
  
 	typedef itk::Image< short, Dimension >             ShortImageType;
  	typedef ShortImageType::Pointer                    ShortImagePointer;

    typedef itk::Image< float, Dimension >             FloatImageType;
  	typedef FloatImageType::Pointer                    FloatImagePointer;
    typedef itk::ImageRegionIterator< FloatImageType > FloatIterator;

  	std::cout<<"Read input image \n";
    FloatImagePointer inputImage = btk::ImageHelper<FloatImageType>::ReadImage(input_file);
    FloatImagePointer maskImage  = btk::ImageHelper<FloatImageType>::ReadOrCreateImage(mask_file, inputImage, 1);

  	std::cout<<"Create empty output image\n";
    ShortImagePointer outputImage = btk::ImageHelper<FloatImageType,ShortImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

    std::vector< FloatImagePointer > probaImages = btk::ImageHelper<FloatImageType>::ReadImage(proba_file);
    std::vector< FloatImagePointer > probaImagesStepK = btk::ImageHelper<FloatImageType>::ReadImage(proba_file);
    std::vector< FloatImagePointer > probaImagesStepKPlusOne = btk::ImageHelper<FloatImageType>::ReadImage(proba_file);

  	//compute characteristics of the input image
	ShortImageType::RegionType  region  = inputImage->GetLargestPossibleRegion();
	ShortImageType::SizeType    size    = region.GetSize();
  	ShortImageType::SpacingType spacing = inputImage->GetSpacing();

    //--------------------------------------------------------------------------------------


    //Estimation of noise level
    btk::Noise<float> myNoiseTool;
    myNoiseTool.SetSigma2Image( btk::ImageHelper<FloatImageType, FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer()) );
    myNoiseTool.ComputeGlobalSigma(inputImage, maskImage);

    //--------------------------------------------------------------------------------------

    //Initialize patch tools
    btk::PatchTool2<float> myPatchTool2;
    myPatchTool2.ComputePatchSize(inputImage, hwn);
    myPatchTool2.ComputeSpatialBandwidth(inputImage, hwvs);

    double smoothing = 2 * beta * myNoiseTool.GetGlobalSigma2() * myPatchTool2.GetFullPatchSize()[0] * myPatchTool2.GetFullPatchSize()[1] * myPatchTool2.GetFullPatchSize()[2];
    std::cout<<"Smoothing = "<<smoothing<<", estimated variance = "<<myNoiseTool.GetGlobalSigma2()<<std::endl;

    FloatImagePointer meanImage     = btk::ImageHelper<FloatImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());
    FloatImagePointer varianceImage = btk::ImageHelper<FloatImageType,FloatImageType>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

    myPatchTool2.ComputeMeanAndVarianceImage(inputImage, maskImage, meanImage, varianceImage);

    int maxNeighbours = myPatchTool2.GetFullPatchSize()[0] * myPatchTool2.GetFullPatchSize()[1] * myPatchTool2.GetFullPatchSize()[2];


    //-------------------------------------------------------------------------------------------------------------------------------------

    unsigned int numberOfClasses = probaImages.size();

    std::cout<<"Normalize probability images\n";
    btk::PandoraBoxImageFilters::ProbabilityImageNormalization(probaImages,probaImages);
    btk::PandoraBoxImageFilters::ProbabilityImageNormalization(probaImagesStepK,probaImagesStepK);
    btk::PandoraBoxImageFilters::ProbabilityImageNormalization(probaImagesStepKPlusOne,probaImagesStepKPlusOne);

    std::cout<<"Initialize output image using the maximum probability at each voxel\n";
    btk::PandoraBoxImageFilters::GetLabelWithMaxProbabilityImage(probaImages, outputImage);


    int itermax = 100;
    for(int i=0; i < itermax; i++)
    {
      std::cout<<"Iteration : "<<i<<std::endl;
      int x,y,z;
      #pragma omp parallel for private(x,y,z) schedule(dynamic)
      for(z=0; z < (int)size[2]; z++)
      for(y=0; y < (int)size[1]; y++)
      for(x=0; x < (int)size[0]; x++)
      {
        FloatImageType::IndexType p;
        p[0] = x;
        p[1] = y;
        p[2] = z;
        if( maskImage->GetPixel(p) > 0 )
        {
            //Get neighbours
            std::vector< FloatImageType::IndexType > neighbours;
            neighbours.reserve(maxNeighbours);
            myPatchTool2.GetNeighboursUsingMeanAndVariance(p, inputImage, meanImage->GetPixel(p), varianceImage->GetPixel(p), meanImage, varianceImage, neighbours);
            //myPatchTool2.GetNeighbours(p, inputImage, neighbours);

            //Compute new neighbours with updated labels
            std::vector< FloatImageType::IndexType > neighboursWithUpdatedLabels;
            for(unsigned int n=0; n<neighbours.size(); n++)
              if( outputImage->GetPixel(neighbours[n]) > 0)
                neighboursWithUpdatedLabels.push_back(neighbours[n]);

            //Compute the corresponding weights
            std::vector<double> weights;
            double sumOfWeights = 0;
            //Reserve the (known) size of the vector (should be faster)
            weights.reserve(neighboursWithUpdatedLabels.size());
            myPatchTool2.ComputeNeighbourWeights(inputImage, p, neighboursWithUpdatedLabels, weights, sumOfWeights, smoothing);

            //Propagate labels to the current voxel
            for(unsigned int l = 0; l < numberOfClasses; l++)
            {
              double localConsistency = 0;

              for(unsigned int n=0; n<neighboursWithUpdatedLabels.size(); n++)
                localConsistency += weights[n]*probaImagesStepK[l]->GetPixel( neighboursWithUpdatedLabels[n] );

              if( sumOfWeights > 0 )
                localConsistency /= sumOfWeights;

              probaImagesStepKPlusOne[l]->SetPixel(p, localConsistency);
            }
        }
      }
      for(unsigned int l = 0; l < numberOfClasses; l++)
        probaImagesStepK[l] = btk::ImageHelper<FloatImageType>::DeepCopy( probaImagesStepKPlusOne[l].GetPointer() );

      btk::PandoraBoxImageFilters::GetLabelWithMaxProbabilityImage(probaImagesStepKPlusOne, outputImage);

      if(i%5==0)
      {
        std::string tmpfile;
        std::ostringstream oss;
        oss << i;
        tmpfile = "output_"+oss.str()+".nii.gz";
        btk::ImageHelper<ShortImageType>::WriteImage(outputImage, tmpfile);
        for(unsigned int l = 0; l < numberOfClasses; l++)
        {
          std::ostringstream oss2;
          oss2 << l;
          tmpfile = "proba_"+oss2.str()+"_"+oss.str()+".nii.gz";
          btk::ImageHelper<FloatImageType>::WriteImage(probaImagesStepKPlusOne[l], tmpfile);

        }

      }
    }


    //-------------------------------------------------------------------------------------------------------------------------------------
    btk::ImageHelper<ShortImageType>::WriteImage(outputImage, output_file);
    btk::ImageHelper<FloatImageType>::WriteImage(probaImagesStepKPlusOne[0], "proba0.nii.gz");


    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------

    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
