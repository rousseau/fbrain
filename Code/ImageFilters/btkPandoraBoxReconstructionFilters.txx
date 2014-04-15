/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 31/01/2014
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

#include "btkPandoraBoxReconstructionFilters.h"

namespace btk
{

void PandoraBoxReconstructionFilters::Convert3DImageToSliceStack(std::vector<itkFloatImagePointer> & outputStack, itkFloatImagePointer & inputImage)
{
  unsigned int numberOfSlices = inputImage->GetLargestPossibleRegion().GetSize()[2];
  outputStack.resize(numberOfSlices);

  //Need to create a new region for each slice (ITK requires a size and an index for defining a new region)
  itkFloatImage::SizeType sliceSize;
  sliceSize[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
  sliceSize[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
  sliceSize[2] = 1;

  itkFloatImage::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  itkFloatImage::RegionType sliceRegion;
  sliceRegion.SetSize( sliceSize );
  sliceRegion.SetIndex( start );

  //Loop over each slice
  for(unsigned int i=0; i < numberOfSlices; i++)
  {
    outputStack[i] = itkFloatImage::New();
    outputStack[i]->SetRegions( sliceRegion );
    outputStack[i]->SetSpacing( inputImage->GetSpacing() );
    outputStack[i]->SetDirection( inputImage->GetDirection() );

    //Compute the new origin location
    itkFloatImage::PointType newOriginPoint;  //physical point location of the new origin
    itkFloatImage::IndexType newOriginIndex;  //index of the new origin
    newOriginIndex[0] = 0;
    newOriginIndex[1] = 0;
    newOriginIndex[2] = i;

    inputImage->TransformIndexToPhysicalPoint(newOriginIndex,newOriginPoint);
    outputStack[i]->SetOrigin( newOriginPoint );

    //Allocate and Copy image values
    outputStack[i]->Allocate();

    itkFloatImage::IndexType inputIndex;
    inputIndex[2] = i;
    itkFloatImage::IndexType outputIndex;
    outputIndex[2] = 0;

    for(unsigned int x=0; x<inputImage->GetLargestPossibleRegion().GetSize()[0]; x++)
    for(unsigned int y=0; y<inputImage->GetLargestPossibleRegion().GetSize()[1]; y++)
    {
      inputIndex[0] = x;
      inputIndex[1] = y;
      outputIndex[0] = x;
      outputIndex[1] = y;
      outputStack[i]->SetPixel(outputIndex, inputImage->GetPixel(inputIndex) );
    }
  }
}

void PandoraBoxReconstructionFilters::ImageFusionByInjection(itkFloatImagePointer & outputImage, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkTransformType::Pointer> > & affineSBSTransforms)
{
  std::cout<<"Starting ImageFusionByInjection"<<std::endl;
  itkFloatImage::PointType outputPoint;      //physical point in HR output image
  itkFloatImage::IndexType outputIndex;      //index in HR output image
  itkFloatImage::PointType transformedPoint; //Physical point location after applying affine transform
  itkContinuousIndex       inputContIndex;   //continuous index in LR image
  itkContinuousIndex       interpolationContIndex;   //continuous index in LR image for interpolation (i.e. z = 0)
  interpolationContIndex[2] = 0;

  //Create a weight image from the output image
  itkFloatImage::Pointer weightImage = btk::ImageHelper< itkFloatImage > ::CreateNewImageFromPhysicalSpaceOf(outputImage,0.0);

  //Set iterator for output and weight images
  itkFloatIteratorWithIndex itOuputImage(outputImage,outputImage->GetLargestPossibleRegion());
  itkFloatIterator          itWeightImage(weightImage,weightImage->GetLargestPossibleRegion());

  //Define a threshold for z coordinate based on FWHM = 2sqrt(2ln2)sigma = 2.3548 sigma
  float cst = 2*sqrt(2*log(2.0));
  float sigmaz = inputStacks[0][0]->GetSpacing()[2] /cst;
  int scale_search_Z = 8;
  float sz2 = sigmaz * scale_search_Z;

  std::cout<<"sz2 : "<<sz2<<std::endl;

  //ITK Interpolator
  itk::LinearInterpolateImageFunction<itkFloatImage, double>::Pointer interpolator = itk::LinearInterpolateImageFunction<itkFloatImage, double>::New();

  for(itOuputImage.GoToBegin(), itWeightImage.GoToBegin(); !itOuputImage.IsAtEnd(); ++itOuputImage, ++itWeightImage)
  {
    //Coordinate in the output image (index)
    outputIndex = itOuputImage.GetIndex();

    //Coordinate in mm (physical point)
    outputImage->TransformIndexToPhysicalPoint(outputIndex,outputPoint);

    //Loop over the input stacks
    for(unsigned int s=0; s<inputStacks.size(); s++)
    {
      unsigned int sizeX = inputStacks[s][0]->GetLargestPossibleRegion().GetSize()[0];
      unsigned int sizeY = inputStacks[s][0]->GetLargestPossibleRegion().GetSize()[1];

      //Loop over images of the current stack
      for(unsigned int i=0; i<inputStacks[s].size(); i++)
      {
        //Coordinate in mm in the current image (physical point)
        transformedPoint = affineSBSTransforms[s][i]->TransformPoint(outputPoint);

        //Coordinate in the current image (continuous index)
        inputStacks[s][i]->TransformPhysicalPointToContinuousIndex(transformedPoint,inputContIndex);

        //Check whether point is inside the ROI (2D image size and slice to distance less than a given threshold)
        if( (inputContIndex[0] >= 0) && (inputContIndex[1] >= 0) && (inputContIndex[0] < sizeX) && (inputContIndex[1] < sizeY) && (fabs(inputContIndex[2]) <= sz2) )
        {
          //Set input for interpolator
          interpolator->SetInputImage( inputStacks[s][i] );

          interpolationContIndex[0] = inputContIndex[0];
          interpolationContIndex[1] = inputContIndex[1];
          float pixelValue = interpolator->EvaluateAtContinuousIndex( interpolationContIndex );

          //Compute weight and corresponding intensity
          float weight = exp(-inputContIndex[2]*inputContIndex[2]/(2*sigmaz*sigmaz))/(sqrt(2*M_PI)*sigmaz);

          float newValue = itOuputImage.Get() + weight * pixelValue;
          float newWeight = weight + itWeightImage.Get();

          //Set computed values
          itOuputImage.Set(newValue);
          itWeightImage.Set(newWeight);
        }
      }
    }
  }

  //Divide the output image by the weight image
  for(itOuputImage.GoToBegin(), itWeightImage.GoToBegin(); !itOuputImage.IsAtEnd(); ++itOuputImage, ++itWeightImage)
  {
    if( itWeightImage.Get() > 0)
    {
      float newValue = itOuputImage.Get() / itWeightImage.Get();
      itOuputImage.Set( newValue );
    }
  }
  btk::ImageHelper<itkFloatImage>::WriteImage(weightImage, "weight.nii.gz");
}

} // namespace btk
