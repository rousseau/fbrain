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

#include "btkImageHelper.h"

#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include <tclap/CmdLine.h>
#include <string>

const unsigned int Dimension = 3;
typedef int 	VoxelType;

typedef itk::Image <VoxelType,Dimension>		ImageType;
typedef itk::BinaryBallStructuringElement<VoxelType,Dimension>	StructuringElementType;
typedef itk::GrayscaleMorphologicalClosingImageFilter <ImageType, ImageType, StructuringElementType> ClosingFilterType;

typedef ImageType::Pointer 	ImagePointer;

typedef btk::ImageHelper<ImageType> HelperType;

int main(int argc, char **argv)
{
	try
	{
		//TCLAP Command Line Parser
		TCLAP::CmdLine cmd("Greyscale Morphological Closing by a ball structuring element", ' ', "0.1");
		
		//TCLAP Arguments
		TCLAP::ValueArg<std::string> inputImageArg("i","image_file","input image file (short)",true,"","string");
		cmd.add( inputImageArg );
		TCLAP::ValueArg<std::string> outputImageArg("o","output_file","output image (short)",true,"","string");
		cmd.add( outputImageArg );
		TCLAP::ValueArg<int> radiusArg("r","struct_radius","radius of the structuring element (which is a ball here)",true,0,"int");
		cmd.add( radiusArg );
		
		// Parse the args.
		cmd.parse( argc, argv );
		
		//Get Arguments
		std::string inputFile = inputImageArg.getValue();
		std::string outputFile = outputImageArg.getValue();
		unsigned int radius = radiusArg.getValue();
		
		//Read image
		ImagePointer image = HelperType::ReadImage(inputFile);
		
		//Set Structuring Element
		StructuringElementType structElement;
		structElement.SetRadius(radius);
		structElement.CreateStructuringElement();
		
		//Create and Run Closing Filter
		ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();
		closingFilter->SetInput(image);
		closingFilter->SetKernel(structElement);
		
		//Write Ouput
		HelperType::WriteImage(closingFilter->GetOutput(), outputFile);
	}
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{ 
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; 
	}
	
	return 0;
}
