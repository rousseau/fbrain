/*==========================================================================

© Université de Strasbourg - Centre National de la Recherche Scientifique

Date: 19/07/2012
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
#include "btkImageHelper.h"

#include "itkBlackTopHatImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkFlatStructuringElement.h"

#include <tclap/CmdLine.h>

int main(int argc, char **argv)
{
	try
	{
		//TCLAP Command Line Parser
		TCLAP::CmdLine cmd("Tissue Classification to retrieve the cortex of a fetal brain MRI", ' ', "0.1");
		
		//TCLAP Arguments
		TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
		cmd.add( inputImageArg );
		TCLAP::ValueArg<std::string> initSegArg("s","init_seg","2 Class FCM classification of the intracranian volume (short)",true,"","string");
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
		
		typedef itk::NormalizeImageFilter<GreyImageType, NormaliseImageType> 	NormaliseFilterType;
		
		typedef itk::FlatStructuringElement<Dimension>										FlatStructuringElementType;
		typedef itk::BlackTopHatImageFilter<GreyImageType, GreyImageType, FlatStructuringElementType>				TopHatFilterType;
		typedef itk::GrayscaleMorphologicalClosingImageFilter<LabelImageType, LabelImageType, FlatStructuringElementType>	ClosingFilterType;
		
		typedef btk::TopologicalKMeans <VectorGreyImageType, LabelImageType> 		GreyTopologicalKMeansType;
		typedef btk::TopologicalKMeans <VectorNormaliseImageType, LabelImageType> 	NormaliseTopologicalKMeansType;
		
		//Get Input images
		GreyImageType::Pointer greyImage;
		greyImage = GreyHelperType::ReadImage(inputFile);
		
		LabelImageType::Pointer initSegmentation;
		initSegmentation = LabelHelperType::ReadImage(initSegFile);
		
		LabelImageType::Pointer brainstemImage;
		brainstemImage = LabelHelperType::ReadImage(brainstemFile);
		
		LabelImageType::Pointer cerebellumImage;
		cerebellumImage = LabelHelperType::ReadImage(cerebellumFile);
		
		/* ------------------------------------------------- Brain Segmentation -------------------------------------------------*/
		//Create Vector Image for Topological K-Means Input
		GreyImageToVectorImageFilterType::Pointer greyImageToVectorImageFilter = GreyImageToVectorImageFilterType::New();
		greyImageToVectorImageFilter->SetInput(0, greyImage);
		greyImageToVectorImageFilter->Update();
		
		VectorGreyImageType::Pointer greyInput;
		greyInput = greyImageToVectorImageFilter->GetOutput();
		
		//Run Topological K-Means
		GreyTopologicalKMeansType::Pointer greyTopologicalKMeansFilter = GreyTopologicalKMeansType::New();
		greyTopologicalKMeansFilter->SetInputImage(greyInput);
		greyTopologicalKMeansFilter->SetInitialSegmentation(initSegmentation);
		greyTopologicalKMeansFilter->SetLCR();
		greyTopologicalKMeansFilter->Update();
		
		LabelImageType::Pointer brainSegmentation = greyTopologicalKMeansFilter->GetOutput();
		
		//Closing of the LCR to fill the holes
		LabelImageType::Pointer lcr = GreyTopologicalKMeansType::GetOneLabel(brainSegmentation, 1);
		
		FlatStructuringElementType::RadiusType radius;
		radius.Fill(3);
		FlatStructuringElementType structElement = FlatStructuringElementType::Ball(radius);
		
		ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
		closingFilter->SetInput(lcr);
		closingFilter->SetKernel(structElement);
		closingFilter->Update();
		LabelImageType::Pointer lcrCorrected = closingFilter->GetOutput();
		
		LabelImageIteratorType brainSegmentationIterator(brainSegmentation, brainSegmentation->GetLargestPossibleRegion());
		LabelImageIteratorType lcrCorrectedIterator(lcrCorrected, lcrCorrected->GetLargestPossibleRegion());
		
		for(brainSegmentationIterator.GoToBegin(), lcrCorrectedIterator.GoToBegin(); !brainSegmentationIterator.IsAtEnd(); ++brainSegmentationIterator, ++lcrCorrectedIterator)
		{
			if(brainSegmentationIterator.Get() != 0 && lcrCorrectedIterator.Get() == 1)
				brainSegmentationIterator.Set(1);
		}
		
		/* --------------------------------------------------- Cortex Segmentation --------------------------------------------------*/
		//Create Black Top Hat Image
		for(unsigned int i=0; i<3; i++)
			radius.SetElement(0, 2/greyImage->GetSpacing()[0]);
		structElement = FlatStructuringElementType::Ball(radius);
		
		TopHatFilterType::Pointer topHatFilter = TopHatFilterType::New();
		topHatFilter->SetInput(greyImage);
		topHatFilter->SetKernel(structElement);
		topHatFilter->Update();
		
		GreyImageType::Pointer blackTopHatImage = topHatFilter->GetOutput();
		
		//Normalise Grey and Top-Hat Images
		NormaliseFilterType::Pointer normaliseFilter = NormaliseFilterType::New();
		
		normaliseFilter->SetInput(greyImage);
		normaliseFilter->Update();
		NormaliseImageType::Pointer normaliseGreyImage = normaliseFilter->GetOutput();
		
		normaliseFilter->SetInput(blackTopHatImage);
		normaliseFilter->Update();
		NormaliseImageType::Pointer normaliseBlackTopHatImage = normaliseFilter->GetOutput();
		
		//Builds Vector Image
		NormaliseImageToVectorImageFilterType::Pointer normaliseImageToVectorImageFilter = NormaliseImageToVectorImageFilterType::New(); std::cout<<"Update sur le vecteur d'image fait"<<std::endl;
		normaliseImageToVectorImageFilter->SetInput(0, normaliseGreyImage);
		normaliseImageToVectorImageFilter->SetInput(1, normaliseBlackTopHatImage);
		normaliseImageToVectorImageFilter->Update();
		
		VectorNormaliseImageType::Pointer normaliseInput = normaliseImageToVectorImageFilter->GetOutput();
		
		
		
		//Write outputs
		LabelHelperType::WriteImage(brainSegmentation, brainSegmentationFile);
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}
	
	return 0;
}