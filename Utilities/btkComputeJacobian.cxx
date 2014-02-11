/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 21/01/2014
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

/* ITK */
#include "itkImage.h"
#include "itkDisplacementFieldTransform.h"
#include "itkDisplacementFieldJacobianDeterminantFilter.h"

/* BTK */

#include "btkImageHelper.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>


int main(int argc, char * argv[])
{
    const unsigned int Dimension = 3;

    typedef itk::Image<float, Dimension>                    itkFloatImage;
	typedef itk::Vector< float, Dimension >                itkDeformationPixelType;
	typedef itk::Image<itkDeformationPixelType, Dimension>  itkDisplacementFieldType;
    //typedef itk::DisplacementFieldTransform< float, Dimension > itkDisplacementFieldType;
    //typedef itk::DisplacementFieldTransform< float,3 >::DisplacementFieldType itkDeformationField;

    typedef itk::DisplacementFieldJacobianDeterminantFilter<itkDisplacementFieldType, float, itkFloatImage>   itkFilterType;
	


    //TCLAP
	try { 
	
    TCLAP::CmdLine cmd("Computes the jacobian of a displacement field", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputArg("i","input","input deformation field",true,"","string",cmd);
    TCLAP::ValueArg<std::string> outputArg("o","output","output jacobian image",true,"","string",cmd);

    // Parse the argv array.
    cmd.parse( argc, argv );
    
    std::string inputFile  = inputArg.getValue();
    std::string outputFile = outputArg.getValue();


	//typedef itk::ImageFileReader< DeformationFieldType > FieldReaderType;
    //FieldReaderType::Pointer fieldReader = FieldReaderType::New();

    //fieldReader->SetFileName( argv[3] );

    //itkDeformationField::Pointer inputDisplacementField = btk::ImageHelper< itkDeformationField >::ReadImage(inputFile);
    itkDisplacementFieldType::Pointer inputDisplacementField = btk::ImageHelper< itkDisplacementFieldType >::ReadImage(inputFile);
    
    itkFilterType::Pointer filter = itkFilterType::New();
  	filter->SetInput( inputDisplacementField );
	filter->Update();

  	itkFloatImage::Pointer outputImage = filter->GetOutput();	
    btk::ImageHelper<itkFloatImage>::WriteImage(outputImage, outputFile);    

		
    } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  

  return EXIT_SUCCESS;

}
