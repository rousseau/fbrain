/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

20 May 2011
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

#include <tclap/CmdLine.h>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkImage.h"

#include "../Code/Segmentation/btkLabelPropagationTool.h"

#include <vector>
#include <sstream>

int main(int argc, char** argv)
{
  try {  

  TCLAP::CmdLine cmd("Label propagation", ' ', "1.0", true);

  TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
  cmd.add( inputImageArg );
  TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image file (short)",true,"","string");
  cmd.add( outputImageArg );
  
  TCLAP::MultiArg<std::string> labelImageArg("l","label_file","label image of the textbook (short) (possible multiple inputs) ",true,"string");
  cmd.add( labelImageArg ); 
  TCLAP::MultiArg<std::string> anatomicalImageArg("a","anatomical_file","anatomical image of the textbook (short) (possible multiple inputs) ",true,"string");
  cmd.add( anatomicalImageArg );
  
  TCLAP::ValueArg<std::string> inputMaskArg("m","mask_file","filename of the mask image",false,"","string");
  cmd.add( inputMaskArg );
  TCLAP::ValueArg<std::string> outputWeightArg("w","weight_file","filename of the output weight map",false,"","string");
  cmd.add( outputWeightArg );
  TCLAP::MultiArg< int > outputFuzzyWeightsArg("","fuzzy_weights","set of labels for fuzzy weight output 3D images (possible multiple inputs -- only using the block mode)",false,"int");
  cmd.add( outputFuzzyWeightsArg );   
  TCLAP::ValueArg<std::string> inputHwnArg("","hwn_map","filename of the input hwn map",false,"","string");
  cmd.add( inputHwnArg );
  TCLAP::ValueArg<std::string> inputHwvsArg("","hwvs_map","filename of the hwvs map",false,"","string");
  cmd.add( inputHwvsArg );
  TCLAP::ValueArg< float > paddingArg("p","pad","padding value (used if no mask image is provided), default is -1",false,-1,"float");
  cmd.add( paddingArg );
  TCLAP::ValueArg< int > hwnArg("","hwn","patch half size (default is 1)",false,1,"int");
  cmd.add( hwnArg );
  TCLAP::ValueArg< int > hwvsArg("","hwvs","half size of the volume search area, i.e. the spatial bandwidth (default is 5)",false,5,"int");
  cmd.add( hwvsArg );
  TCLAP::ValueArg< float > betaArg("b","beta","beta: smoothing parameter (high beta produces smoother result, default is 1)",false,1,"float");
  cmd.add( betaArg );
  TCLAP::ValueArg< int > blockArg("","block","0: pointwise, 1: blockwise, 2: fast blockwise (default is 2)",false,2,"int");
  cmd.add( blockArg );
  TCLAP::ValueArg< int > centerArg("c","center","weight of the central patch (possible value: 0, 1, -1 (max)) (default is -1)",false,-1,"int");
  cmd.add( centerArg );
  TCLAP::ValueArg< int > optimizedArg("","opt","optimized mode (use mean and standard deviation of patches) (0: no, 1: yes) (default is 1)",false,1,"int");
  cmd.add( optimizedArg );
  TCLAP::ValueArg< float > lowerMeanThresholdArg("","lmt","lower mean threshold (0.95 by default) -- for optimized mode only",false,0.95,"float");
  cmd.add( lowerMeanThresholdArg );
  TCLAP::ValueArg< float > lowerVarianceThresholdArg("","lvt","lower variance threshold (0.5 by default) -- for optimized mode only",false,0.5,"float");
  cmd.add( lowerVarianceThresholdArg );
  TCLAP::ValueArg< int > aggregationArg("","aggregation","aggregation strategy (0: max weight, 1: fixed weight (=1))",false,1,"int");
  cmd.add( aggregationArg );
  TCLAP::ValueArg< int > modeArg("","mode","mode (0: pair-wise, 1: groupwise, -1: groupwise denoising, -2: HR simulation)",false,0,"int");
  cmd.add( modeArg );
  TCLAP::ValueArg< int > normalizationArg("","normalization","intensity normalization -- makes the patches invariant to intensity changes (0: no, 1: yes)",false,0,"int");
  cmd.add( normalizationArg );
  TCLAP::ValueArg< int > minLabelArg("","minLabel","minimum label (default : -1)",false,-1,"int");
  cmd.add( minLabelArg );
  
  // Parse the args.
  cmd.parse( argc, argv );

    
  // Get the value parsed by each arg. 
  std::string input_file       = inputImageArg.getValue();
  std::string output_file      = outputImageArg.getValue();
  
  std::vector<std::string> anatomical_file = anatomicalImageArg.getValue();
  std::vector<std::string> label_file      = labelImageArg.getValue();

  std::string mask_file        = inputMaskArg.getValue();
  std::string weight_file      = outputWeightArg.getValue();
  std::vector<int> fuzzy_weights = outputFuzzyWeightsArg.getValue();

  std::string hwn_file         = inputHwnArg.getValue();   
  std::string hwvs_file        = inputHwvsArg.getValue();      
  
  float padding                = paddingArg.getValue();
  int hwn                      = hwnArg.getValue();
  int hwvs                     = hwvsArg.getValue();
  float beta                   = betaArg.getValue();
  int block                    = blockArg.getValue();
  int center                   = centerArg.getValue();
  int optimized                = optimizedArg.getValue();
  float lowerMeanThreshold     = lowerMeanThresholdArg.getValue();
  float lowerVarianceThreshold = lowerVarianceThresholdArg.getValue();
  int aggregation              = aggregationArg.getValue();
  int mode                     = modeArg.getValue();
  int normalization            = normalizationArg.getValue();
  int minLabel                 = minLabelArg.getValue();
  
  std::cout<<" input file : "<<input_file<<"\n";
  std::cout << "Number of anatomical image files is: " << anatomical_file.size() << "\n";
  for(unsigned int i=0; i<anatomical_file.size();i++)
    std::cout<<"   anatomical image file "<<i+1<<" : "<<anatomical_file[i]<<"\n";
  std::cout << "Number of label image files is: " << label_file.size() << "\n";
  for(unsigned int i=0; i<label_file.size();i++)
    std::cout<<"   label image file "<<i+1<<" : "<<label_file[i]<<"\n";   
  
  if(mask_file=="")
    std::cout<<"No mask image file given. The padding value ("<<padding<<") will be used instead.\n";
  else
    std::cout<<"   mask image file "<<mask_file<<"\n";  
  
  //ITK declaration
  typedef short PixelType;
  const   unsigned int        Dimension = 3;
  typedef itk::Image< PixelType, Dimension >    ImageType; //same type for input and output
  typedef ImageType::Pointer ImagePointer;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;

  ImagePointer maskImage;
  LabelFusionTool<PixelType> myTool;

  myTool.ReadInput(input_file);
  myTool.ReadAnatomicalImages(anatomical_file);
  myTool.ReadLabelImages(label_file);

  if(mask_file=="")                    //creating a mask image using the padding value
	  myTool.SetPaddingValue(padding); 
  else{                                //reading the mask image
	  ReaderType::Pointer maskReader = ReaderType::New();
	  maskReader->SetFileName( mask_file );
	  maskReader->Update();
	  maskImage = maskReader->GetOutput();
	  myTool.SetMaskImage(maskImage);
  }
		  
	  
  if (hwn_file!="")
    myTool.ReadOneImage(hwn_file, myTool.m_hwnImage);
  if (hwvs_file!="")
    myTool.ReadOneImage(hwvs_file, myTool.m_hwvsImage);
	  
  myTool.SetPatchSize(hwn);
  myTool.SetSpatialBandwidth(hwvs);
  myTool.SetSmoothing(beta);
  myTool.SetCentralPointStrategy(center);
  myTool.SetBlockwiseStrategy(block);
  myTool.SetAggregationStrategy(aggregation);
  myTool.SetOptimizationStrategy(optimized);  
  myTool.SetLowerThresholds(lowerMeanThreshold, lowerVarianceThreshold);
  myTool.SetNormalizationStrategy(normalization);
  myTool.SetMinLabel(minLabel);

  if(mode==0)
    myTool.ComputeOutput();            //pair-wise label propagation
  if(mode==1)
    myTool.ComputeOutput_groupwise();  //groupwise label propagation
  if(mode==-1)
    myTool.ComputeDenoisedOutput();    //denoising using the other anatomical images
    if(mode==-2)
      myTool.ComputeHROutput();        //HR estimation using the HR anatomical images

  //Write the result 
  ImagePointer outputImage = ImageType::New();
  myTool.GetOutput(outputImage);
  WriterType::Pointer writer = WriterType::New();  
  writer->SetFileName( output_file );
  writer->SetInput( outputImage );
  writer->Update();  

	  
  if (weight_file != "") {
    std::cout<<"Writing the weight image.\n";
    ImagePointer weightImage = ImageType::New();
    myTool.GetWeightImage(weightImage);
    WriterType::Pointer writerFilter = WriterType::New();  
    writerFilter->SetFileName( weight_file );
    writerFilter->SetInput( weightImage );
    writerFilter->Update();  
  }
	   
	if( (fuzzy_weights.size() > 0) && (block==1) ){
		for(unsigned int i=0; i<fuzzy_weights.size(); i++){
			int current_label = fuzzy_weights[i];

			itk::Image< float, 3 >::Pointer fuzzyWeightImage;
			myTool.InitImage(fuzzyWeightImage);
  		itk::ImageFileWriter< itk::Image< float, 3 > >::Pointer writerFilter = itk::ImageFileWriter< itk::Image< float, 3 > >::New();

			myTool.GetFuzzyWeightImage(fuzzyWeightImage, current_label);

			std::string outputFilename = "fuzzyWeight_" ;
			std::string extension = ".nii.gz";
			std::ostringstream oss;
      oss << current_label ;
			outputFilename += oss.str() + extension;
	    writerFilter->SetFileName( outputFilename );
  	  writerFilter->SetInput( fuzzyWeightImage );
    	writerFilter->Update();  
		}
	}
  
  return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
