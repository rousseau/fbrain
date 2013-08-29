/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 07/01/2013
  Author(s): François Rousseau

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
#include "itkVector.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkPoint.h"

/* BTK */
#include "btkImageHelper.h"

/* OTHERS */
#include "iostream"
#include <tclap/CmdLine.h>

int main(int argc, char** argv)
{
  try {  

    TCLAP::CmdLine cmd("Compute Soft Mask using orthogonal images", ' ', "1.0", true);
    TCLAP::MultiArg<std::string> inputImageArg("i","input_file","input image file (float)",true,"string",cmd);
    TCLAP::MultiArg<std::string> outputImageArg("o","output_file","output image file (float)",true,"string",cmd);
    TCLAP::ValueArg<int>         valueArg("v","value","type of soft masks: 0=mask, 1=mean image (default:0)",false,0,"int",cmd);
    TCLAP::ValueArg<float>       distanceArg("d","distance","distance (to the barycentre) for mask thresholding",false,10000,"float",cmd);
    TCLAP::ValueArg<float>       thresholdArg("t","threshold","threshold on the dot product between z vector",false,0.5,"float",cmd);

    cmd.parse( argc, argv );
    
    // Get the value parsed by each arg. 
    std::vector<std::string> input_file       = inputImageArg.getValue();
    std::vector<std::string> output_file      = outputImageArg.getValue();
    int                      mask_value       = valueArg.getValue();
    float                    maxDistance      = distanceArg.getValue();
    float                    threshold        = thresholdArg.getValue();
    //ITK declaration
    const   unsigned int Dimension = 3;
    typedef itk::Image< float, Dimension >    FloatImage;
    typedef FloatImage::Pointer               FloatImagePointer;
    typedef FloatImage::Pointer               FloatImagePointer;
    typedef FloatImage::IndexType             IndexType;
    typedef itk::ContinuousIndex<double,3>    ContinuousIndex;


    std::cout<<"Reading the input images ..."<<std::endl;
    
    std::vector< FloatImagePointer > inputImages;
    inputImages = btk::ImageHelper<FloatImage>::ReadImage(input_file);
    
    std::cout<<"Compute the physical z vectors"<<std::endl;
    std::vector< itk::Vector<double, 3> > physicalZVectors;
    
    for(unsigned int i=0; i!=inputImages.size(); i++ )
    {
      FloatImage::IndexType startingIndex;  
      FloatImage::PointType startingPoint;  
      FloatImage::IndexType endingIndex;  
      FloatImage::PointType endingPoint;
      
      startingIndex[0] = 0;
      startingIndex[1] = 0;
      startingIndex[2] = 0;
      inputImages[i] -> TransformIndexToPhysicalPoint(startingIndex, startingPoint);

      endingIndex[0]   = 0;
      endingIndex[1]   = 0;
      endingIndex[2]   = 1;
      inputImages[i] -> TransformIndexToPhysicalPoint(endingIndex, endingPoint);
      
      itk::Vector<double, 3> physicalVector;
      physicalVector[0] = startingPoint[0] - endingPoint[0];
      physicalVector[1] = startingPoint[1] - endingPoint[1];
      physicalVector[2] = startingPoint[2] - endingPoint[2];
      
      
      double norm = sqrt(physicalVector * physicalVector);
      physicalVector[0] /= norm;
      physicalVector[1] /= norm;
      physicalVector[2] /= norm;
      
      physicalZVectors.push_back(physicalVector);
       
    }
    
    std::cout<<"Computing the soft mask for each input image ..."<<std::endl;
    std::vector< FloatImagePointer > outputImages;
    outputImages = btk::ImageHelper<FloatImage, FloatImage>::CreateNewImageFromPhysicalSpaceOf(inputImages);

    typedef itk::NearestNeighborInterpolateImageFunction< FloatImage >  NNInterpolatorType;

    typedef itk::ImageRegionIterator< FloatImage > itkFloatIterator;


    std::cout<<"Define a set of interpolators to use the ITK function isInsideBuffer "<<std::endl;
    std::vector< NNInterpolatorType::Pointer > nn_interpolator;
    for(unsigned int i=0; i!=inputImages.size(); i++ )    
    {
      NNInterpolatorType::Pointer nn = NNInterpolatorType::New();
      nn -> SetInputImage( inputImages[i] );
      nn_interpolator.push_back(nn);
    }
      

    //Can be done in a more elegant way...
    
    for(unsigned int i=0; i!=inputImages.size(); i++ )    
    {
      std::cout<<"Find the orthogonal images wrt the current image ..."<<std::endl;
      std::vector<unsigned int> orthogonalImages;
      
      for(unsigned int j=0; j!=inputImages.size(); j++)
      {
        double dotProduct = physicalZVectors[i] * physicalZVectors[j];
        if( std::abs(dotProduct) < threshold )
        {
          orthogonalImages.push_back(j);
        }
      }
      
      //Go through the output image ... ----------------------------------------
      itkFloatIterator outputImageIt( outputImages[i], outputImages[i] -> GetLargestPossibleRegion() );
      for(outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++outputImageIt )
      {            
        float outputValue = 0;
        
        
        FloatImage::IndexType index = outputImageIt.GetIndex();
        FloatImage::PointType point;
        outputImages[i] -> TransformIndexToPhysicalPoint(index, point);

        for(unsigned int j=0; j!=orthogonalImages.size(); j++)
        {

          if ( nn_interpolator[ orthogonalImages[j] ] -> IsInsideBuffer( point ) )
          {
            if(mask_value==0)
            {
              outputValue += 1.0;
            }  
            else
            {
              ContinuousIndex contIndex;
              outputImages[ orthogonalImages[j] ] ->TransformPhysicalPointToContinuousIndex(point, contIndex);     
              outputValue += nn_interpolator[ orthogonalImages[j] ]->EvaluateAtContinuousIndex(contIndex);
            }              
          }

        }
        
        outputImageIt.Set( outputValue / orthogonalImages.size() );
      }      
    }    


    //Now, compute the barycentre for each image and apply a distance-based threshold

    for(unsigned int i=0; i!=inputImages.size(); i++ )    
    {
      FloatImage::PointType barycentre;
      barycentre[0] = 0;
      barycentre[1] = 0;
      barycentre[2] = 0;
      double currentIndex = 0;
      itkFloatIterator outputImageIt( outputImages[i], outputImages[i] -> GetLargestPossibleRegion() );
      for(outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++outputImageIt )
        if( outputImageIt.Get() > 0 )
        {
          FloatImage::IndexType index = outputImageIt.GetIndex();
          FloatImage::PointType point;
          outputImages[i] -> TransformIndexToPhysicalPoint(index, point);
          double weight = (1.0*currentIndex)/(1.0*(currentIndex+1.0));
          barycentre[0] = weight * barycentre[0] + (1.0-weight) * point[0];
          barycentre[1] = weight * barycentre[1] + (1.0-weight) * point[1];
          barycentre[2] = weight * barycentre[2] + (1.0-weight) * point[2];          
          currentIndex += 1;
          
        }
        
      std::cout<<"Barycentre of image "<<i<<":\n "<<barycentre<<std::endl;
        
      //Threshold the outputImage based on the distance to the barycentre
      for(outputImageIt.GoToBegin(); !outputImageIt.IsAtEnd(); ++outputImageIt )
        if( outputImageIt.Get() > 0 )
        {
          FloatImage::IndexType index = outputImageIt.GetIndex();
          FloatImage::PointType point;
          outputImages[i] -> TransformIndexToPhysicalPoint(index, point);

          if( point.SquaredEuclideanDistanceTo(barycentre) > maxDistance )
            outputImageIt.Set(0);
        }  
            
    }  

        
    std::cout<<"Writing the output soft masks ... "<<std::endl;
    btk::ImageHelper<FloatImage>::WriteImage(outputImages, output_file);
  
  
  
    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
}
