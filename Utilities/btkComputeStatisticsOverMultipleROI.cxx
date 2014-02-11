/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 16/01/2014
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
#include "itkImageRegionIterator.h"

/* BTK */

#include "btkImageHelper.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>


int main(int argc, char * argv[])
{
    const unsigned int Dimension = 3;

    typedef itk::Image<float, Dimension> itkFloatImage;
    typedef itk::Image<short, Dimension> itkShortImage;

    typedef itk::ImageRegionIterator< itkFloatImage > itkFloatIterator;
    typedef itk::ImageRegionIterator< itkShortImage > itkShortIterator;

    //TCLAP
	try { 
	
    TCLAP::CmdLine cmd("Computes the average and the variance in every ROI defined by a label image", ' ', "1.0", true);
    
    TCLAP::ValueArg<std::string> inputArg("i","input","input image file of interest",true,"","string",cmd);
    TCLAP::ValueArg<std::string> labelArg("l","label","label image defining the ROI",true,"","string",cmd);
    TCLAP::ValueArg<std::string> meanArg("m","mean","mean image for each ROI",true,"","string",cmd);
    TCLAP::ValueArg<std::string> varianceArg("v","variance","variance image for each ROI",true,"","string",cmd);

    // Parse the argv array.
    cmd.parse( argc, argv );
    
    std::string inputFile = inputArg.getValue();
    std::string labelFile = labelArg.getValue();
    std::string meanFile = meanArg.getValue();
    std::string varianceFile = varianceArg.getValue();

    itkFloatImage::Pointer inputImage = btk::ImageHelper< itkFloatImage >::ReadImage(inputFile);
    itkShortImage::Pointer labelImage = btk::ImageHelper< itkShortImage >::ReadImage(labelFile);

    itkFloatImage::Pointer meanImage = btk::ImageHelper<itkFloatImage>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());
    itkFloatImage::Pointer varianceImage = btk::ImageHelper<itkFloatImage>::CreateNewImageFromPhysicalSpaceOf(inputImage.GetPointer());

	std::map<short,double> mMap;
	std::map<short,double> m2Map;
	std::map<short,long unsigned int> counterMap;
	std::map<short,double> meanMap;
	std::map<short,double> varianceMap;
	
	itkFloatIterator inputIterator(inputImage, inputImage->GetRequestedRegion() );
	itkShortIterator labelIterator(labelImage, labelImage->GetRequestedRegion() );	
	itkFloatIterator meanIterator(meanImage, meanImage->GetRequestedRegion() );
	itkFloatIterator varianceIterator(varianceImage, varianceImage->GetRequestedRegion() );
	
  	for(inputIterator.GoToBegin(),labelIterator.GoToBegin(); !inputIterator.IsAtEnd(); ++inputIterator, ++labelIterator){
	  float inputValue = inputIterator.Get();
	  short labelValue = labelIterator.Get();
	  
	  counterMap[labelValue] += 1;	
	  mMap[labelValue] += inputValue;
	  m2Map[labelValue]+= (inputValue * inputValue);	  
    } 	
    
    std::map<short, double>::iterator mMapIt;
    std::map<short, double>::iterator m2MapIt;
    std::map<short, long unsigned int>::iterator counterMapIt;
    
    for(mMapIt = mMap.begin(), m2MapIt = m2Map.begin(), counterMapIt = counterMap.begin(); mMapIt != mMap.end(); ++mMapIt, ++m2MapIt, ++counterMapIt){
      short label = (*mMapIt).first;
      long unsigned int n = (*counterMapIt).second;
      
      double mean = (*mMapIt).second / n;      
      double variance = ( (*m2MapIt).second / n) - (mean * mean) ;

	  meanMap[label]     = mean;
	  varianceMap[label] = variance; 
    }

  	for(labelIterator.GoToBegin(),meanIterator.GoToBegin(),varianceIterator.GoToBegin(); !labelIterator.IsAtEnd(); ++labelIterator, ++meanIterator, ++varianceIterator){
	  short labelValue = labelIterator.Get();
	 
	  meanIterator.Set( meanMap[labelValue] );
	  varianceIterator.Set( varianceMap[labelValue] );
    } 	
    
    btk::ImageHelper<itkFloatImage>::WriteImage(meanImage, meanFile);
    btk::ImageHelper<itkFloatImage>::WriteImage(varianceImage, varianceFile);    

	
		
		
    } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  

  return EXIT_SUCCESS;

}
