/*==========================================================================

© Université de Strasbourg - Centre National de la Recherche Scientifique

Date: 07/08/2012
Author(s): Benoit Caldairou (benoit.caldairou@unistra.fr)

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

#include "btkTopologicalKMeans.h"
#include "btkFCMClassifier.h"
#include "btkImageHelper.h"

#include "itkBlackTopHatImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkMaskImageFilter.h"

#include "itkVectorIndexSelectionCastImageFilter.h"

#include <tclap/CmdLine.h>

//typedef
const unsigned int Dimension = 3;
typedef int 	GreyVoxelType;
typedef float	NormaliseVoxelType;
typedef int 	LabelVoxelType;

typedef itk::Image <GreyVoxelType,Dimension>			GreyImageType;
typedef itk::Image <NormaliseVoxelType, Dimension>		NormaliseImageType;
typedef itk::VectorImage <GreyVoxelType,Dimension>		VectorGreyImageType;
typedef itk::VectorImage <NormaliseVoxelType, Dimension> 	VectorNormaliseImageType;
typedef itk::Image <LabelVoxelType,Dimension>			LabelImageType;

typedef itk::ImageRegionIterator <LabelImageType>		LabelImageIteratorType;

typedef btk::ImageHelper <GreyImageType> 			GreyHelperType;
typedef btk::ImageHelper <LabelImageType> 			LabelHelperType; 

typedef itk::ComposeImageFilter <GreyImageType> 		GreyImageToVectorImageFilterType;
typedef itk::ComposeImageFilter <NormaliseImageType> 		NormaliseImageToVectorImageFilterType;

typedef itk::NormalizeImageFilter<GreyImageType, NormaliseImageType> 		NormaliseFilterType;
typedef itk::MaskImageFilter<GreyImageType, LabelImageType, GreyImageType>	MaskFilterType;

typedef itk::FlatStructuringElement<Dimension>										FlatStructuringElementType;
typedef itk::BlackTopHatImageFilter<GreyImageType, GreyImageType, FlatStructuringElementType>				TopHatFilterType;
typedef itk::GrayscaleMorphologicalClosingImageFilter<LabelImageType, LabelImageType, FlatStructuringElementType>	ClosingFilterType;

typedef btk::TopologicalKMeans <VectorGreyImageType, LabelImageType> 		GreyTopologicalKMeansType;
typedef btk::TopologicalKMeans <VectorNormaliseImageType, LabelImageType> 	NormaliseTopologicalKMeansType;

typedef btk::FCMClassifier<GreyImageType, LabelImageType, VectorNormaliseImageType> FCMClassifierType;

void lcrClosing(LabelImageType::Pointer segmentation);

int main(int argc, char **argv)
{
	try
	{
		/* ------------------------------------------------------------- Sets TCLAP and get input images ----------------------------------------------------------------------------- */
		//TCLAP Command Line Parser
		TCLAP::CmdLine cmd("Tissue Classification to retrieve the cortex of a fetal brain MRI", ' ', "0.1");
		
		//TCLAP Arguments
		TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
		cmd.add( inputImageArg );
		TCLAP::ValueArg<std::string> initSegArg("s","intra_vol","segmentation of the intracranian volume (short)",true,"","string");
		cmd.add( initSegArg );
		TCLAP::MultiArg<std::string> manualSegArg("l","brainstem_cerebellum","manual segmentation of brainstem and cerebellum (short)",true,"",cmd);
		TCLAP::MultiArg<std::string> outputImageArg("o","output_files","Sets output files : one for brain segmentation, the other for cortex segmentation (short)",true,"",cmd);
		
		//Parse command line
		cmd.parse(argc, argv);
		
		//Check how many outputs were set
		if(outputImageArg.getValue().size() != 2 || manualSegArg.getValue().size() != 2)
		{
			std::cout<<"There should be 2 outputs, one file for brainstem segmentation and one file for cerebellum segmentation"<<std::endl;
			return 1;
		}
		
		//Get arguments
		std::string inputFile = inputImageArg.getValue();
		std::string initSegFile = initSegArg.getValue();
		std::string brainstemFile = manualSegArg.getValue()[0];
		std::string cerebellumFile = manualSegArg.getValue()[1];
		std::string brainSegmentationFile = outputImageArg.getValue()[0];
		std::string cortexSegmentationFile = outputImageArg.getValue()[1];
		
		//Get Input images
		std::cout<<"Loads Input Images..."<<std::endl;
		
		GreyImageType::Pointer greyImage;
		greyImage = GreyHelperType::ReadImage(inputFile);
		
		LabelImageType::Pointer initSegmentation;
		initSegmentation = LabelHelperType::ReadImage(initSegFile);
		
		LabelImageType::Pointer brainstemImage;
		brainstemImage = LabelHelperType::ReadImage(brainstemFile);
		
		LabelImageType::Pointer cerebellumImage;
		cerebellumImage = LabelHelperType::ReadImage(cerebellumFile);
		
		std::cout<<"Input Images Loaded"<<std::endl<<std::endl;
		
		/* ----------------------------------------------------------------------------------FCM Classification -------------------------------------------------------------------------- */
		std::cout<<"2 class fcm classification..."<<std::endl;
		FCMClassifierType::Pointer fcmClassifier = FCMClassifierType::New();
		fcmClassifier->SetGreyImage(greyImage);
		fcmClassifier->SetMaskImage(initSegmentation);
		fcmClassifier->SetClassNumber(2);
		fcmClassifier->Update();
		
		LabelImageType::Pointer fcmSegmentation = fcmClassifier->GetLabelSegmentation();
		
		std::cout<<"fcm classification done"<<std::endl<<std::endl;
		
		/* --------------------------------------------------------------------------------- Brain Segmentation ---------------------------------------------------------------------------*/
		//Create Vector Image for Topological K-Means Input
		std::cout<<"Creates Vector Image..."<<std::endl;
		
		GreyImageToVectorImageFilterType::Pointer greyImageToVectorImageFilter = GreyImageToVectorImageFilterType::New();
		greyImageToVectorImageFilter->SetInput(0, greyImage);
		greyImageToVectorImageFilter->Update();
		
		VectorGreyImageType::Pointer greyInput;
		greyInput = greyImageToVectorImageFilter->GetOutput();
		
		std::cout<<"Vector Image created"<<std::endl<<std::endl;
		
		//Run Topological K-Means
		std::cout<<"Segmentation of the brain..."<<std::endl;
		
		GreyTopologicalKMeansType::Pointer greyTopologicalKMeansFilter = GreyTopologicalKMeansType::New();
		greyTopologicalKMeansFilter->SetInputImage(greyInput);
		greyTopologicalKMeansFilter->SetInitialSegmentation(fcmSegmentation);
		greyTopologicalKMeansFilter->SetLCR();
		greyTopologicalKMeansFilter->Update();
		
		LabelImageType::Pointer brainSegmentation = greyTopologicalKMeansFilter->GetOutput();
		
		std::cout<<"Brain Segmented"<<std::endl<<std::endl;
		
		//Closing of the LCR to fill the holes
		std::cout<<"Fills holes in LCR..."<<std::endl;
		lcrClosing(brainSegmentation);
		std::cout<<"Holes filed"<<std::endl<<std::endl;
		
		/* ---------------------------------------------------------------------------------- Cortex Segmentation -------------------------------------------------------------------------------*/
		//Create Black Top Hat Image
		std::cout<<"Creates black top hat image..."<<std::endl;
		
		FlatStructuringElementType::RadiusType radius;
		for(unsigned int i=0; i<3; i++)
			radius.SetElement(i, 2/greyImage->GetSpacing()[i]);
		FlatStructuringElementType structElement = FlatStructuringElementType::Ball(radius);
		
		TopHatFilterType::Pointer topHatFilter = TopHatFilterType::New();
		topHatFilter->SetInput(greyImage);
		topHatFilter->SetKernel(structElement);
		topHatFilter->Update();
		
		GreyImageType::Pointer blackTopHatImage = topHatFilter->GetOutput();
		
		std::cout<<"Black top hat image created"<<std::endl<<std::endl;
		
		//Mask Images by intracranian volume to get a better normalisation
		std::cout<<"Normalisation of grey and top hat images..."<<std::endl;
		
		MaskFilterType::Pointer maskFilter1 = MaskFilterType::New();
		maskFilter1->SetInput(greyImage);              
		maskFilter1->SetMaskImage(brainSegmentation);
		maskFilter1->Update();
		
		MaskFilterType::Pointer maskFilter2 = MaskFilterType::New();
		maskFilter2->SetInput(blackTopHatImage);
		maskFilter2->SetMaskImage(brainSegmentation);
		maskFilter2->Update();
		
		//Normalise Grey and Top-Hat Images
		NormaliseFilterType::Pointer normaliseFilter1 = NormaliseFilterType::New();
		NormaliseFilterType::Pointer normaliseFilter2 = NormaliseFilterType::New();
		
		normaliseFilter1->SetInput(maskFilter1->GetOutput());
		normaliseFilter1->Update();

		normaliseFilter2->SetInput(maskFilter2->GetOutput());
		normaliseFilter2->Update();
		
		std::cout<<"Normalisation done"<<std::endl<<std::endl;
		
		//Builds Vector Image
		std::cout<<"Creates Vector Image..."<<std::endl;
		
		NormaliseImageToVectorImageFilterType::Pointer normaliseImageToVectorImageFilter = NormaliseImageToVectorImageFilterType::New();
		normaliseImageToVectorImageFilter->SetInput(0, normaliseFilter1->GetOutput());
		normaliseImageToVectorImageFilter->SetInput(1, normaliseFilter2->GetOutput());
		normaliseImageToVectorImageFilter->Update();
		
		VectorNormaliseImageType::Pointer normaliseInput = normaliseImageToVectorImageFilter->GetOutput();
		
		std::cout<<"Vector Image created"<<std::endl<<std::endl;
		
		//Run Topological K-Means for Cortex
		std::cout<<"Segmentation of the cortex..."<<std::endl;
		
		NormaliseTopologicalKMeansType::Pointer normaliseTopologicalKMeansFilter = NormaliseTopologicalKMeansType::New();
		normaliseTopologicalKMeansFilter->SetInputImage(normaliseInput);
		normaliseTopologicalKMeansFilter->SetInitialSegmentation(brainSegmentation);
		normaliseTopologicalKMeansFilter->SetBrainstemImage(brainstemImage);
		normaliseTopologicalKMeansFilter->SetCerebellumImage(cerebellumImage);
		normaliseTopologicalKMeansFilter->SetCortex();
		normaliseTopologicalKMeansFilter->Update();
		LabelImageType::Pointer cortexSegmentation = normaliseTopologicalKMeansFilter->GetOutput();
		
		std::cout<<"Cortex segmented"<<std::endl<<std::endl;
		
		//Closing of the LCR to fill the holes
		std::cout<<"Fills holes in LCR..."<<std::endl;
		lcrClosing(cortexSegmentation);
		std::cout<<"Holes filed"<<std::endl<<std::endl;
		
		/* ------------------------------------------------------------------------- Outputs Writing ---------------------------------------------------------------------------------------------*/
		std::cout<<"Saves output images..."<<std::endl;
		LabelHelperType::WriteImage(brainSegmentation, brainSegmentationFile);
		LabelHelperType::WriteImage(cortexSegmentation, cortexSegmentationFile);
		std::cout<<"Output images saved"<<std::endl<<std::endl;
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}
	
	return 0;
}

void lcrClosing(LabelImageType::Pointer segmentation)
{
	LabelImageType::Pointer lcr = GreyTopologicalKMeansType::GetOneLabel(segmentation, 1);
	
	FlatStructuringElementType::RadiusType radius;
	for(unsigned int i=0; i<Dimension; i++)
	{
		radius.SetElement(i, 2/segmentation->GetSpacing()[i]);
	}
	FlatStructuringElementType structElement = FlatStructuringElementType::Ball(radius);
	
	ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
	closingFilter->SetInput(lcr);
	closingFilter->SetKernel(structElement);
	closingFilter->Update();
	LabelImageType::Pointer lcrCorrected = closingFilter->GetOutput();
	
	LabelImageIteratorType segmentationIterator(segmentation, segmentation->GetLargestPossibleRegion());
	LabelImageIteratorType lcrCorrectedIterator(lcrCorrected, lcrCorrected->GetLargestPossibleRegion());
	
	for(segmentationIterator.GoToBegin(), lcrCorrectedIterator.GoToBegin(); !segmentationIterator.IsAtEnd(); ++segmentationIterator, ++lcrCorrectedIterator)
	{
		if(segmentationIterator.Get() != 0 && lcrCorrectedIterator.Get() == 1)
			segmentationIterator.Set(1);
	}
}
