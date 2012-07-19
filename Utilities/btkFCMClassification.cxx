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

#include "btkFCMClassifier.h"
#include "btkImageHelper.h"

#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include <tclap/CmdLine.h>
#include <string>

const unsigned int Dimension = 3;
typedef int 	GreyVoxelType;
typedef int 	LabelVoxelType;
typedef float 	FuzzyVoxelType;

typedef itk::Image <GreyVoxelType,Dimension>		GreyImageType;
typedef itk::Image <LabelVoxelType,Dimension>		LabelImageType;
typedef itk::VectorImage <FuzzyVoxelType,Dimension>	FuzzyImageType;
typedef itk::Image <FuzzyVoxelType,Dimension>		FuzzyScalarType;

typedef GreyImageType::Pointer 	GreyImagePointer;
typedef LabelImageType::Pointer LabelImagePointer;
typedef FuzzyImageType::Pointer FuzzyImagePointer;

typedef itk::ImageFileReader<GreyImageType> GreyReaderType;
typedef itk::ImageFileWriter<GreyImageType> GreyWriterType;

typedef btk::ImageHelper<GreyImageType> GreyHelperType;
typedef btk::ImageHelper<LabelImageType> LabelHelperType; 

typedef btk::FCMClassifier <GreyImageType, LabelImageType, FuzzyImageType> 			FCMClassifierType;
typedef itk::VectorIndexSelectionCastImageFilter<FuzzyImageType, FuzzyScalarType> 		AdaptorType;

int main(int argc, char **argv)
{
	try
	{
		//TCLAP Command Line Parser
		TCLAP::CmdLine cmd("FCM Classification Algorithm", ' ', "0.1");
		
		//TCLAP Arguments
		TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
		cmd.add( inputImageArg );
		TCLAP::ValueArg<std::string> maskImageArg("m","mask_file","mask image file (short)",true,"","string");
		cmd.add( maskImageArg );
		TCLAP::ValueArg<std::string> labelImageArg("l","label_file","label image file (short)",true,"","string");
		cmd.add( labelImageArg );
		TCLAP::ValueArg<std::string> fuzzyImageArg("f","fuzzy_files","fuzzy maps files pattern (for example : fuzzyMaps.nii will turn to fuzzyMaps1.nii, fuzzyMaps2.nii and fuzzyMaps3.nii for a 3 class clustering)(short)",true,"","string");
		cmd.add( fuzzyImageArg );
		TCLAP::ValueArg<int> classNumberArg("n","class_number","number of class to be looked for",true,0,"int");
		cmd.add( classNumberArg );
		
		// Parse the args.
		cmd.parse( argc, argv );
		
		//Get Arguments
		std::string input_file = inputImageArg.getValue();
		std::string label_file = labelImageArg.getValue();
		std::string mask_file = maskImageArg.getValue();
		std::string fuzzy_files = fuzzyImageArg.getValue();
		unsigned int classNumber = classNumberArg.getValue();
		
// 		for(unsigned int i=0; i<classNumber; i++)
// 		{
// 			std::string fuzzyFile = fuzzy_files;
// 			std::stringstream fileNumber; fileNumber<<i;
// 			std::cout<<fuzzyFile<<" "<<fileNumber.str()<<std::endl;
// 			fuzzyFile.insert(fuzzyFile.find_first_of('.'), fileNumber.str());
// 			std::cout<<fuzzyFile<<std::endl;
// 		}
		
		//Read Grey and Mask Images
		GreyImagePointer greyImage = GreyHelperType::ReadImage(input_file);
		LabelImagePointer maskImage = LabelHelperType::ReadImage(mask_file);
		
		//Set FCMClassifier parameters
		FCMClassifierType::Pointer fcmClassifier = FCMClassifierType::New();
		fcmClassifier->SetGreyImage(greyImage);
		fcmClassifier->SetMaskImage(maskImage);
		fcmClassifier->SetClassNumber(classNumber);
		fcmClassifier->Update();
		
		//Write Ouputs
		FuzzyImageType::Pointer fuzzyMaps = fcmClassifier->GetFuzzyMaps();
		AdaptorType::Pointer adaptor = AdaptorType::New();
		adaptor->SetInput(fuzzyMaps);
		for(unsigned int i=0; i<classNumber; i++)
		{
			adaptor->SetIndex(i);
			std::string fuzzyFile = fuzzy_files;
			std::stringstream fileNumber; fileNumber<<(i+1);
			fuzzyFile.insert(fuzzyFile.find_first_of('.'), fileNumber.str());
			btk::ImageHelper< FuzzyScalarType >::WriteImage(adaptor->GetOutput(),fuzzyFile);
		}
		
		return 0;
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{ 
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}
}
