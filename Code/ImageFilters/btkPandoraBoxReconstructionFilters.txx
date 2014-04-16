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

void PandoraBoxReconstructionFilters::ConvertSliceStackTo3DImage(itkFloatImagePointer & outputImage, std::vector<itkFloatImagePointer> & inputStack)
{
  //Image information are computed from the first slice of the stack
  outputImage = itkFloatImage::New();
  outputImage->SetSpacing( inputStack[0]->GetSpacing() );
  outputImage->SetDirection( inputStack[0]->GetDirection() );
  outputImage->SetOrigin( inputStack[0]->GetOrigin() );

  itkFloatImage::SizeType imageSize;
  imageSize[0] = inputStack[0]->GetLargestPossibleRegion().GetSize()[0];
  imageSize[1] = inputStack[0]->GetLargestPossibleRegion().GetSize()[1];
  imageSize[2] = inputStack.size();

  itkFloatImage::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  itkFloatImage::RegionType imageRegion;
  imageRegion.SetSize( imageSize );
  imageRegion.SetIndex( start );
  outputImage->SetRegions( imageRegion );

  //Allocate and Copy image values
  outputImage->Allocate();

  itkFloatImage::IndexType inputIndex;
  inputIndex[2] = 0;
  itkFloatImage::IndexType outputIndex;

  for(unsigned int i=0; i<inputStack.size(); i++)
  for(unsigned int x=0; x<inputStack[i]->GetLargestPossibleRegion().GetSize()[0]; x++)
  for(unsigned int y=0; y<inputStack[i]->GetLargestPossibleRegion().GetSize()[1]; y++)
  {
    inputIndex[0] = x;
    inputIndex[1] = y;
    outputIndex[0] = x;
    outputIndex[1] = y;
    outputIndex[2] = i;
    outputImage->SetPixel(outputIndex, inputStack[i]->GetPixel(inputIndex) );
  }
}

void PandoraBoxReconstructionFilters::Project3DImageToSliceStack(std::vector<itkFloatImagePointer> & outputStack, itkFloatImagePointer & inputImage, std::vector<itkFloatImagePointer> & inputStack, std::vector<itkTransformType::Pointer> & affineSBSTransforms)
{
  outputStack.resize( inputStack.size() );

  //Use currently a linear interpolation
  itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>::Pointer bsInterpolator = itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>::New();
  bsInterpolator->SetSplineOrder(1);
  bsInterpolator->SetInputImage(inputImage);

  //Loop over the input stacks
  for(unsigned int s=0; s<outputStack.size(); s++)
  {
    //Allocate and initialize the output slice
    outputStack[s] = btk::ImageHelper< itkFloatImage > ::CreateNewImageFromPhysicalSpaceOf(inputStack[s],0.0);

    itkFloatImage::IndexType sliceIndex;      //index in the current slice
    sliceIndex[2] = 0;
    itkFloatImage::PointType slicePoint;      //physical point of the current slice
    itkFloatImage::PointType transformedPoint; //Physical point location after applying affine transform
    itkContinuousIndex       inputContIndex;   //continuous index in the 3D image

    //Loop over the pixel of the current slice
    for(unsigned int x=0; x<outputStack[s]->GetLargestPossibleRegion().GetSize()[0]; x++)
    for(unsigned int y=0; y<outputStack[s]->GetLargestPossibleRegion().GetSize()[1]; y++)
    {
      //Coordinate (index) of the current pixel in the current slice
      sliceIndex[0] = x;
      sliceIndex[1] = y;

      //Coordinate in the physical world (mm)
      outputStack[s]->TransformIndexToPhysicalPoint(sliceIndex,slicePoint);

      //Apply affine transform (slice to 3D image)
      transformedPoint = affineSBSTransforms[s]->TransformPoint(slicePoint);

      //Coordinate in the 3D image (continuous index)
      inputImage->TransformPhysicalPointToContinuousIndex(transformedPoint,inputContIndex);

      double interpolatedValue = bsInterpolator->EvaluateAtContinuousIndex(inputContIndex);

      outputStack[s]->SetPixel(sliceIndex, interpolatedValue);
    }
  }
}

void PandoraBoxReconstructionFilters::ComputePSFImage(itkFloatImagePointer & PSFImage, itkFloatImage::SpacingType HRSpacing, itkFloatImage::SpacingType LRSpacing)
{
  //Use a 3D Gaussian in this function
  float cst = 2*sqrt(2*log(2.0));
  float sigmaX = 1.2 * LRSpacing[0] / cst;
  float sigmaY = 1.2 * LRSpacing[1] / cst;
  float sigmaZ = LRSpacing[2] / cst;

  //To get 99%, psfSize must be >= 6*sigma
  //To get 95%, psfSize must be >= 4*sigma
  float scaleFactor = 4;

  //Compute size and index for the new PSF image
  itkFloatImage::SizeType psfSize;
  psfSize[0] = (int)ceil(scaleFactor*sigmaX / HRSpacing[0]) + 2;
  psfSize[1] = (int)ceil(scaleFactor*sigmaY / HRSpacing[1]) + 2;
  psfSize[2] = (int)ceil(scaleFactor*sigmaZ / HRSpacing[2]) + 2;

  itkFloatImage::IndexType psfIndex;
  psfIndex[0] = 0;
  psfIndex[1] = 0;
  psfIndex[2] = 0;

  itkFloatImage::RegionType psfRegion;
  psfRegion.SetSize(psfSize);
  psfRegion.SetIndex(psfIndex);

  PSFImage = itkFloatImage::New();
  PSFImage->SetRegions(psfRegion);
  PSFImage->SetSpacing(HRSpacing);
  PSFImage->Allocate();
  PSFImage->FillBuffer(0.0);

  //Set the origin of the PSF in the center of it
  itkContinuousIndex psfIndexCenter;
  psfIndexCenter[0] = (psfSize[0]-1)/2.0;
  psfIndexCenter[1] = (psfSize[1]-1)/2.0;
  psfIndexCenter[2] = (psfSize[2]-1)/2.0;
  itkFloatImage::PointType newPSFOrigin;
  PSFImage->TransformContinuousIndexToPhysicalPoint(psfIndexCenter,newPSFOrigin);
  newPSFOrigin[0] *= -1;
  newPSFOrigin[1] *= -1;
  newPSFOrigin[2] *= -1;
  PSFImage->SetOrigin(newPSFOrigin);

  itkFloatImage::IndexType currentIndex;      //index of the current voxel
  itkFloatImage::PointType currentPoint;      //physical coordinates of the current voxel
  itkFloatIteratorWithIndex itPSFImage(PSFImage,PSFImage->GetLargestPossibleRegion());

  double psfValue = 0;
  double psfSum = 0;
  //Loop over PSF points and compute the integration of the PSF function (currently 3D Gaussian)
  for(itPSFImage.GoToBegin(); !itPSFImage.IsAtEnd(); ++itPSFImage)
  {
    currentIndex = itPSFImage.GetIndex();

    PSFImage->TransformIndexToPhysicalPoint(currentIndex,currentPoint);

    //should be done by integration instead of evaluating the value at a certain location
    psfValue = exp(-currentPoint[0]*currentPoint[0]/(2*sigmaX*sigmaX))/(sqrt(2*M_PI)*sigmaX);
    psfValue*= exp(-currentPoint[1]*currentPoint[1]/(2*sigmaY*sigmaY))/(sqrt(2*M_PI)*sigmaY);
    psfValue*= exp(-currentPoint[2]*currentPoint[2]/(2*sigmaZ*sigmaZ))/(sqrt(2*M_PI)*sigmaZ);

    itPSFImage.Set(psfValue);
    psfSum += psfValue;
  }

  //Normalization of PSF values
  for(itPSFImage.GoToBegin(); !itPSFImage.IsAtEnd(); ++itPSFImage)
  {
    double normalizedValue = itPSFImage.Get() / psfSum;
    itPSFImage.Set(normalizedValue);
  }
}

void PandoraBoxReconstructionFilters::ComputerObservationModelParameters(vnl_sparse_matrix<float> & H, vnl_vector<float> & Y, vnl_vector<float> & X, itkFloatImagePointer & HRImage, std::vector< std::vector<itkFloatImagePointer> > & maskStacks, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkTransformType::Pointer> > & inverseAffineSBSTransforms, itkFloatImagePointer & PSFImage)
{
  //Principle: for each voxel of the LR images, we compute the influence of each voxel of the PSF (centered at the current LR voxel) and add the corresponding influence value (PSF value * interpolation weight) in the matrix H

  // Set size of matrices
  unsigned int ncols = HRImage->GetLargestPossibleRegion().GetNumberOfPixels();

  std::vector< std::vector<unsigned int> > offset;
  //offset is used to fill correctly the vector Y with the input LR image values (offset for the linear index)
  offset.resize(inputStacks.size());
  unsigned int nrows = 0;
  for(unsigned int im = 0; im < inputStacks.size(); im++){
    offset[im].resize( inputStacks[im].size() );
    for(unsigned int s = 0; s < inputStacks[im].size(); s++)
    {
      offset[im][s] = nrows;
      nrows += inputStacks[im][s]->GetLargestPossibleRegion().GetNumberOfPixels();
    }
  }
  H.set_size(nrows, ncols);
  Y.set_size(nrows);
  Y.fill(0.0);
  X.set_size(ncols);
  X.fill(0.0);

  //linear index : an integer value corresponding to (x,y,z) triplet coordinates (ITK index)
  unsigned int lrLinearIndex = 0;
  unsigned int hrLinearIndex = 0;

  //Temporary variables
  itkFloatImage::IndexType lrIndex;  //index of the current voxel in the LR image
  itkFloatImage::PointType lrPoint;  //physical point location of lrIndex
  itkFloatImage::IndexType psfIndex; //index of the current voxel in the PSF
  itkFloatImage::PointType psfPoint; //physical point location of psfIndex
  itkFloatImage::PointType transformedPoint; //Physical point location after applying affine transform
  itkContinuousIndex  hrContIndex;  //continuous index in HR image of psfPoint
  itkFloatImage::IndexType hrIndex;  //index in HR image of interpolated psfPoint

  //We use linear interpolation for the estimation of point influence in matrix H
  typedef itk::BSplineInterpolationWeightFunction<double, 3, 1> itkBSplineFunction;
  itkBSplineFunction::Pointer bsplineFunction = itkBSplineFunction::New();
  itkBSplineFunction::WeightsType bsplineWeights;
  bsplineWeights.SetSize(8); // (bsplineOrder + 1)^3
  itkBSplineFunction::IndexType   bsplineStartIndex;
  itkBSplineFunction::IndexType   bsplineEndIndex;
  itkBSplineFunction::SizeType    bsplineSize;
  itkFloatImage::RegionType       bsplineRegion;

  //Get the size of the HR image
  itkFloatImage::SizeType  hrSize  = HRImage->GetLargestPossibleRegion().GetSize();

  //Create a moving PSF to compute H
  itkFloatImage::Pointer movingPSFImage = btk::ImageHelper< itkFloatImage > ::DeepCopy(PSFImage.GetPointer());

  //Instantiate an iterator over the current PSF
  itkFloatIteratorWithIndex itPSF(movingPSFImage,movingPSFImage->GetLargestPossibleRegion());

  //Compute temporary variables for setting the correct origin of the moving PSF for every voxel of the LR image
  itkFloatImage::SizeType psfSize = PSFImage->GetLargestPossibleRegion().GetSize();
  itkContinuousIndex psfIndexCenter;
  psfIndexCenter[0] = (psfSize[0]-1)/2.0;
  psfIndexCenter[1] = (psfSize[1]-1)/2.0;
  psfIndexCenter[2] = (psfSize[2]-1)/2.0;
  itkFloatImage::PointType psfPointCenter;
  itkFloatImage::PointType psfOrigin;
  itkFloatImage::PointType zeroOrigin;
  zeroOrigin[0] = 0;
  zeroOrigin[2] = 0;
  zeroOrigin[1] = 0;

  unsigned int count = -1;

  for(unsigned int i=0; i<inputStacks.size(); i++)
  {
    for(unsigned int s=0; s<inputStacks[i].size(); s++)
    {
      //Get the size of the current LR slice
      itkFloatImage::SizeType  lrSize  = inputStacks[i][s]->GetLargestPossibleRegion().GetSize();

      //Instantiate an iterator over the current LR image
      itkFloatIteratorWithIndex itLRImage(inputStacks[i][s],inputStacks[i][s]->GetLargestPossibleRegion());

      //Instantiate an iterator over the current LR mask image
      itkFloatIterator itLRMask(maskStacks[i][s],maskStacks[i][s]->GetLargestPossibleRegion());

      //Set the correct direction for the PSF of the current image
      movingPSFImage->SetDirection(inputStacks[i][s]->GetDirection());

      //Loop over the voxels of the current LR image
      for(itLRImage.GoToBegin(), itLRMask.GoToBegin(); !itLRImage.IsAtEnd(); ++itLRImage, ++itLRMask)
      {
        count++;

        if(itLRMask.Get() > 0)
        {
          //Coordinate in the current LR image
          lrIndex = itLRImage.GetIndex();

          //Compute the corresponding linear index of lrIndex
          lrLinearIndex = offset[i][s] + lrIndex[0] + lrIndex[1]*lrSize[0];

          //Fill Y
          Y[lrLinearIndex] = itLRImage.Get();

          //Reset PSF origin
          movingPSFImage->SetOrigin( zeroOrigin );

          //Change the origin of the PSF so that itLRImage location corresponds to the center of the PSF
          movingPSFImage->TransformContinuousIndexToPhysicalPoint(psfIndexCenter,psfPointCenter);
          inputStacks[i][s]->TransformIndexToPhysicalPoint(lrIndex,lrPoint);
          psfOrigin[0] = lrPoint[0] - psfPointCenter[0];
          psfOrigin[1] = lrPoint[1] - psfPointCenter[1];
          psfOrigin[2] = lrPoint[2] - psfPointCenter[2];
          movingPSFImage->SetOrigin(psfOrigin);

          //Loop over the m_PSF
          for(itPSF.GoToBegin(); !itPSF.IsAtEnd(); ++itPSF)
          {
            //Get coordinate in PSF
            psfIndex = itPSF.GetIndex();

            //Compute the physical point of psfIndex
            movingPSFImage->TransformIndexToPhysicalPoint(psfIndex,psfPoint);

            //Apply estimated affine transform to psfPoint (need to apply the inverse since the transform goes from the HR image to the LR image)
            transformedPoint = inverseAffineSBSTransforms[i][s]->TransformPoint(psfPoint);

            //Get back to the index in the HR image
            HRImage->TransformPhysicalPointToContinuousIndex(transformedPoint,hrContIndex);

            //Check if the continuous index hrContIndex is inside the HR image
            if (HRImage->GetLargestPossibleRegion().IsInside(hrContIndex))
            {
              //Get the interpolation weight using itkBSplineInterpolationWeightFunction
              bsplineFunction->Evaluate(hrContIndex,bsplineWeights,bsplineStartIndex);

              //Get the support size for interpolation
              bsplineSize = bsplineFunction->GetSupportSize();

              //Check if the bspline support region is inside the HR image
              bsplineEndIndex[0] = bsplineStartIndex[0] + bsplineSize[0];
              bsplineEndIndex[1] = bsplineStartIndex[1] + bsplineSize[1];
              bsplineEndIndex[2] = bsplineStartIndex[2] + bsplineSize[2];

              if( (HRImage->GetLargestPossibleRegion().IsInside(bsplineStartIndex)) && (HRImage->GetLargestPossibleRegion().IsInside(bsplineEndIndex)) )
              {
                //Set the support region
                bsplineRegion.SetSize(bsplineSize);
                bsplineRegion.SetIndex(bsplineStartIndex);

                //Instantiate an iterator on HR image over the bspline region
                itkFloatIteratorWithIndex itHRImage(HRImage,bsplineRegion);

                //linear index of bspline weights
                unsigned int weightLinearIndex = 0;

                //Loop over the support region
                for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage)
                {
                  //Get coordinate in HR image
                  hrIndex = itHRImage.GetIndex();

                  //Compute the corresponding linear index
                  hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];

                  //Add weight*PSFValue to the corresponding element in H
                  H(lrLinearIndex,hrLinearIndex) += itPSF.Get() * bsplineWeights[weightLinearIndex];
                  weightLinearIndex += 1;
                }
              }
            }
          }
        }
      }
    }
  }
  // Normalize H
  for(unsigned int i = 0; i < H.rows(); i++){

    double sum = H.sum_row(i);

    vnl_sparse_matrix<float>::row & r = H.get_row(i);
    vnl_sparse_matrix<float>::row::iterator col_iter;

    for (col_iter = r.begin(); col_iter != r.end(); ++col_iter)
      (*col_iter).second = (*col_iter).second / sum;
  }

  //Fill X
  //Instantiate an iterator on HR image
  itkFloatIteratorWithIndex itHRImage(HRImage,HRImage->GetLargestPossibleRegion());
  for(itHRImage.GoToBegin(); !itHRImage.IsAtEnd(); ++itHRImage){
    hrIndex = itHRImage.GetIndex();
    hrLinearIndex = hrIndex[0] + hrIndex[1]*hrSize[0] + hrIndex[2]*hrSize[0]*hrSize[1];
    X[hrLinearIndex] = itHRImage.Get();
  }
}

void PandoraBoxReconstructionFilters::ImageFusionByInjection(itkFloatImagePointer & outputImage, itkFloatImagePointer & maskImage, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkTransformType::Pointer> > & affineSBSTransforms)
{
  //Strictly speaking, this is not an injection process, but il's faster to do it this way
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
  itkFloatIterator          itMaskImage(maskImage,maskImage->GetLargestPossibleRegion());

  //Define a threshold for z coordinate based on FWHM = 2sqrt(2ln2)sigma = 2.3548 sigma
  float cst = 2*sqrt(2*log(2.0));
  float sigmaz = inputStacks[0][0]->GetSpacing()[2] /cst;
  int scale_search_Z = 2;
  float sz2 = sigmaz * scale_search_Z;

  //ITK Interpolator
  itk::LinearInterpolateImageFunction<itkFloatImage, double>::Pointer interpolator = itk::LinearInterpolateImageFunction<itkFloatImage, double>::New();

  for(itOuputImage.GoToBegin(), itWeightImage.GoToBegin(), itMaskImage.GoToBegin(); !itOuputImage.IsAtEnd(); ++itOuputImage, ++itWeightImage, ++itMaskImage)
  {
    if(itMaskImage.Get() > 0)
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
}

void PandoraBoxReconstructionFilters::SimulateObservations(vnl_sparse_matrix<float> & H, vnl_vector<float> & X, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkFloatImagePointer> > & outputStacks)
{

  std::cout<<"SimulateObservations"<<std::endl;
  std::vector< std::vector<unsigned int> > offset;
  //offset is used to fill correctly the vector Y with the input LR image values (offset for the linear index)
  offset.resize(inputStacks.size());
  unsigned int nrows = 0;
  for(unsigned int im = 0; im < inputStacks.size(); im++){
    offset[im].resize( inputStacks[im].size() );
    for(unsigned int s = 0; s < inputStacks[im].size(); s++)
    {
      offset[im][s] = nrows;
      nrows += inputStacks[im][s]->GetLargestPossibleRegion().GetNumberOfPixels();
    }
  }

  //Compute H * x
  vnl_vector<float> Hx;
  H.mult(X,Hx);

  //Now, we have to convert this vector into a set of stacks
  outputStacks.resize(inputStacks.size());
  for(unsigned int im = 0; im < inputStacks.size(); im++)
    outputStacks[im].resize( inputStacks[im].size() );

  //Temporary variables
  itkFloatImage::IndexType lrIndex;  //index of the current voxel in the LR image
  unsigned int lrLinearIndex = 0;

  for(unsigned int i=0; i<inputStacks.size(); i++)
  {
    for(unsigned int s=0; s<inputStacks[i].size(); s++)
    {
      //Initialize every slice of the current stack
      outputStacks[i][s] = btk::ImageHelper< itkFloatImage > ::CreateNewImageFromPhysicalSpaceOf(inputStacks[i][s],0.0);

      //Get the size of the current LR slice
      itkFloatImage::SizeType  lrSize  = outputStacks[i][s]->GetLargestPossibleRegion().GetSize();

      //Instantiate an iterator over the current LR image
      itkFloatIteratorWithIndex itLRImage(outputStacks[i][s],outputStacks[i][s]->GetLargestPossibleRegion());

      for(itLRImage.GoToBegin(); !itLRImage.IsAtEnd(); ++itLRImage)
      {
        //Coordinate in the current LR image
        lrIndex = itLRImage.GetIndex();

        //Compute the corresponding linear index of lrIndex
        lrLinearIndex = offset[i][s] + lrIndex[0] + lrIndex[1]*lrSize[0];

        //Fill Y
        itLRImage.Set( Hx[lrLinearIndex] );

      }
    }
  }
}

} // namespace btk
