/*==========================================================================

© Université de Strasbourg - Centre National de la Recherche Scientifique

Date: 23/07/2012
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

#ifndef BTK_TOPOLOGICAL_KMEANS_TXX
#define BTK_TOPOLOGICAL_KMEANS_TXX

#include "btkTopologicalKMeans.h"
#include "btkImageHelper.h"

#include "itkImageDuplicator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkVariableLengthVector.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkFlatStructuringElement.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkMaskImageFilter.h"
#include "itkSubtractImageFilter.h"

namespace btk
{
	/* --------------------------Constructor------------------------------------------------- */
	template< typename TInputImage, typename TLabelImage>
	TopologicalKMeans<TInputImage, TLabelImage>::TopologicalKMeans()
	{
		//Two Outputs : Fuzzy Maps (as a vectorImage) and LabelSegmentation
		this->SetNumberOfRequiredOutputs(1);
		//Two Intputs : Original Image and Mask Image
		this->SetNumberOfRequiredInputs(2);
	}
	
	/* --------------------------------------Topological K-Means running functions------------------------------------- */
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::GenerateData()
	{
		//Get Inputs
		typename TInputImage::Pointer inputImage = this->GetInputImage();
		typename TLabelImage::Pointer initSeg = this->GetInitialSegmentation();
		std::cout<<"Got Inputs"<<std::endl;
		
		// Setup output (first a copy of the initial Segmentation)
		typename TLabelImage::Pointer finalSeg = this->GetOutput();
		typename itk::ImageDuplicator<TLabelImage>::Pointer duplicator = itk::ImageDuplicator<TLabelImage>::New();
		duplicator->SetInputImage(initSeg);
		duplicator->Update();
		finalSeg = duplicator->GetOutput();
		std::cout<<"Output Set"<<std::endl;
		
		//Algorithm
		InitialiseCentroids(inputImage, finalSeg);
		RunSegmentation(inputImage, finalSeg);
		
		return;
	}
	
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::InitialiseCentroids(typename TInputImage::Pointer inputImage, typename TLabelImage::Pointer segImage)
	{
		// 3 because there is three spheres in this model
		m_Centroids.SetSize(inputImage->GetNumberOfComponentsPerPixel(),3);
		m_Centroids.Fill(0);
		
		itk::ImageRegionConstIterator<TInputImage> inputImageIterator(inputImage, inputImage->GetLargestPossibleRegion());
		itk::ImageRegionConstIterator<TLabelImage> segImageIterator(segImage, segImage->GetLargestPossibleRegion());
		
		itk::VariableLengthVector<unsigned int> sum_voxel;
		sum_voxel.SetSize(inputImage->GetNumberOfComponentsPerPixel());
		sum_voxel.Fill(0);
		
		for(inputImageIterator.GoToBegin(), segImageIterator.GoToBegin(); !inputImageIterator.IsAtEnd(); ++inputImageIterator, ++segImageIterator)
		{
			if(segImageIterator.Get() != 0)
			{
				for(unsigned int i=0; i<m_Centroids.Rows(); i++)
					m_Centroids(i,segImageIterator.Get()-1) += inputImageIterator.Get()[i];
				sum_voxel[segImageIterator.Get()-1]++;
			}
		}
		
		std::cout<<"Nombre de lignes : "<<m_Centroids.Rows()<<std::endl;
		std::cout<<"Nombre de collonnes : "<<m_Centroids.Cols()<<std::endl;
		
		for(unsigned int i=0; i<m_Centroids.Rows(); i++)
			for(unsigned int j=0; j<m_Centroids.Cols(); j++)
				m_Centroids(i,j) = m_Centroids(i,j)/sum_voxel[j];
			
		for(unsigned int i=0; i<m_Centroids.Rows(); i++)
		{
			for(unsigned int j=0; j<m_Centroids.Cols(); j++)
				std::cout<<m_Centroids(i,j)<<" "<<std::flush;
			std::cout<<std::endl;
		}
		
		return;
	}
	
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::RunSegmentation(typename TInputImage::Pointer inputImage, typename TLabelImage::Pointer segImage)
	{
		typename TLabelImage::Pointer borderImage;
		
		unsigned int numLabelChange = 1;
		while(numLabelChange != 0)
		{
			numLabelChange = 0;
			
			borderImage = GetBorderImage(segImage,1,1); //Dilation
			numLabelChange += ClassifyBorderVoxel(inputImage, segImage, borderImage, 1);
			
			borderImage = GetBorderImage(segImage,1,0); //Erosion
			numLabelChange += ClassifyBorderVoxel(inputImage, segImage, borderImage, 1);
			
			borderImage = GetBorderImage(segImage,3,1); //Dilation
			numLabelChange += ClassifyBorderVoxel(inputImage, segImage, borderImage, 3);
			
			borderImage = GetBorderImage(segImage,3,0); //Erosion
			numLabelChange += ClassifyBorderVoxel(inputImage, segImage, borderImage, 3);
			
			ImageHelper<TLabelImage>::WriteImage(segImage,"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Classification_Change.nii");
			
			break;
		}
		
		return;
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::GetBorderImage(typename TLabelImage::Pointer segImage, typename TLabelImage::PixelType label, bool erodeOrDilate)
	{
		//Gets one label
		typename TLabelImage::Pointer labelImage = NULL;
		labelImage = GetOneLabel(segImage, label);
// 		ImageHelper<TLabelImage>::WriteImage(labelImage,"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_Label1.nii");
		
		//Dilate or erode this label
		typename TLabelImage::Pointer morphoImage = NULL;
		if(erodeOrDilate) morphoImage = DilateOneLabel(labelImage);
		else morphoImage = ErodeOneLabel(labelImage);
// 		ImageHelper<TLabelImage>::WriteImage(dilateImage,"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_DilateLabel1.nii");
		
		//Get the border by substracting morphology label and original label
		typename TLabelImage::Pointer diffImage = NULL;
		if(erodeOrDilate) diffImage = SubtractImage(labelImage, morphoImage);
		else diffImage = SubtractImage(morphoImage, labelImage);
// 		ImageHelper<TLabelImage>::WriteImage(diffImage,"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_DiffLabel1.nii");
		
		//Mask this diffImage to eliminate voxel out of intracranial volume
		typename TLabelImage::Pointer borderImage = NULL;
		borderImage = MaskImage(diffImage, segImage);
// 		ImageHelper<TLabelImage>::WriteImage(borderImage,"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_DiffMaskLabel1.nii");
		
		//Check if these border voxels are eligible for a label change
		CheckBorderVoxel(borderImage, segImage, label);
// 		ImageHelper<TLabelImage>::WriteImage(borderImage,"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_DiffMaskCorrectionLabel1.nii");
		
		return borderImage;
	}
	
	template< typename TInputImage, typename TLabelImage>
	unsigned int TopologicalKMeans<TInputImage, TLabelImage>::ClassifyBorderVoxel(typename TInputImage::Pointer inputImage, typename TLabelImage::Pointer segImage, typename TLabelImage::Pointer borderImage, typename TLabelImage::PixelType label)
	{
		itk::ImageRegionConstIterator<TInputImage> inputImageIterator(inputImage, inputImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<TLabelImage> segImageIterator(segImage, segImage->GetLargestPossibleRegion());
		itk::ImageRegionConstIterator<TLabelImage> borderImageIterator(borderImage, borderImage->GetLargestPossibleRegion());
		
		unsigned int pixelChange = 0;
		
		for(inputImageIterator.GoToBegin(), segImageIterator.GoToBegin(), borderImageIterator.GoToBegin(); !inputImageIterator.IsAtEnd(); ++inputImageIterator, ++segImageIterator, ++borderImageIterator)
		{
			if(borderImageIterator.Get() != 0)
			{
				//Vector to compare distance from label 2 and current label centroids
				itk::Vector<float, 2> distanceVector;
				distanceVector.Fill(0);
				
				//Compute distances
				typename TInputImage::PixelType currentPixel = inputImageIterator.Get();
				for(unsigned int i=0; i<currentPixel.GetNumberOfElements(); i++)
				{
					//Distance from label 2 centroids
					distanceVector[0] += pow(currentPixel[i] - m_Centroids(i, 1), 2);
					//Distance from voxel to changed
					distanceVector[1] += pow(currentPixel[i] - m_Centroids(i, label-1), 2);
				}
				
				//Set Label
				//If it changes to label 2
				if(distanceVector[0] < distanceVector[1] && segImageIterator.Get() == label)
				{
					segImageIterator.Set(2);
					pixelChange++;
				}
				//If it changes to current label
				if(distanceVector[0] > distanceVector[1] && segImageIterator.Get() != label)
				{
					segImageIterator.Set(label);
					pixelChange++;
				}
			}
		}
		
		return pixelChange;
	}
	
	/* ----------------------------------------------------- Base functions to get, dilate, erode a label and select the right border voxel ----------------------------------- */
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::CheckBorderVoxel(typename TLabelImage::Pointer borderImage, typename TLabelImage::Pointer segImage, typename TLabelImage::PixelType label)
	{
		typedef itk::ShapedNeighborhoodIterator<TLabelImage> NeighborIteratorType;
		
		//Sets Radius of Neighborhood
		itk::Size<3> radius;
		radius.Fill(1);
		
		//Declares Iterators
		itk::ImageRegionIterator<TLabelImage> borderImageIterator(borderImage, borderImage->GetLargestPossibleRegion());
		NeighborIteratorType segImageNeighborIterator(radius, segImage, segImage->GetLargestPossibleRegion());
		
		//Sets a Cross Neighborhood
		segImageNeighborIterator.ClearActiveList();
		typename NeighborIteratorType::OffsetType top = {{0,0,1}};
		segImageNeighborIterator.ActivateOffset(top);
		typename NeighborIteratorType::OffsetType bottom = {{0,0,-1}};
		segImageNeighborIterator.ActivateOffset(bottom);
		typename NeighborIteratorType::OffsetType left = {{-1,0,0}};
		segImageNeighborIterator.ActivateOffset(left);
		typename NeighborIteratorType::OffsetType right = {{1,0,0}};
		segImageNeighborIterator.ActivateOffset(right);
		typename NeighborIteratorType::OffsetType front = {{0,1,0}};
		segImageNeighborIterator.ActivateOffset(front);
		typename NeighborIteratorType::OffsetType back = {{0,-1,0}};
		segImageNeighborIterator.ActivateOffset(back);
		
		//Sets the opposite label 
		typename TLabelImage::PixelType oppositeLabel = 0;
		if(label == 1) oppositeLabel = 3;
		else oppositeLabel = 1;
		
		//Check border voxel wether it meets conditions to be eligible for label change
		for(borderImageIterator.GoToBegin(), segImageNeighborIterator.GoToBegin(); !borderImageIterator.IsAtEnd(); ++borderImageIterator, ++segImageNeighborIterator)
		{
			if(borderImageIterator.Get() == 1)
			{
				//Gets the list of active neighbors
				typename NeighborIteratorType::IndexListType indexList = segImageNeighborIterator.GetActiveIndexList();
				typename NeighborIteratorType::IndexListType::size_type indexListSize = segImageNeighborIterator.GetActiveIndexListSize();
				typename NeighborIteratorType::IndexListType::iterator indexListIterator;
				
				for(indexListIterator = indexList.begin(); indexListIterator != indexList.end(); indexListIterator++)
				{
					bool isInBound;
					typename TLabelImage::PixelType neighborValue = segImageNeighborIterator.GetPixel(*indexListIterator, isInBound);
					if(!isInBound || neighborValue == 0 || neighborValue == oppositeLabel)
						borderImageIterator.Set(0);
				}
			}
		}
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::GetOneLabel(typename TLabelImage::Pointer segImage, typename TLabelImage::PixelType label)
	{
		//Set threshold Filter
		typedef itk::BinaryThresholdImageFilter<TLabelImage, TLabelImage> ThresholdFilterType;
		typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
		thresholdFilter->SetLowerThreshold(label);
		thresholdFilter->SetUpperThreshold(label);
		thresholdFilter->SetInsideValue(1);
		thresholdFilter->SetOutsideValue(0);
		thresholdFilter->SetInput(segImage);
		thresholdFilter->Update();
		
		return thresholdFilter->GetOutput();
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::DilateOneLabel(typename TLabelImage::Pointer labelImage)
	{
		//Sets structuring element
		typedef itk::FlatStructuringElement<3> StructuringElementType;
		typename StructuringElementType::RadiusType elementRadius;
		elementRadius.Fill(1);
		StructuringElementType structElement = StructuringElementType::Ball(elementRadius);
		std::cout<<"Structuring Element built"<<std::endl;
		
		//Sets Dilate Filter
		typedef itk::GrayscaleDilateImageFilter<TLabelImage, TLabelImage, StructuringElementType> DilateFilterType;
		typename DilateFilterType::Pointer dilateFilter = DilateFilterType::New();
		dilateFilter->SetInput(labelImage);
		dilateFilter->SetKernel(structElement);
		dilateFilter->Update();
// 		ImageHelper<TLabelImage>::WriteImage(dilateFilter->GetOutput(),"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_DilateLabel1_DilateFunction.nii");
		
		return dilateFilter->GetOutput();
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::ErodeOneLabel(typename TLabelImage::Pointer labelImage)
	{
		//Sets structuring element
		typedef itk::FlatStructuringElement<3> StructuringElementType;
		typename StructuringElementType::RadiusType elementRadius;
		elementRadius.Fill(1);
		StructuringElementType structElement = StructuringElementType::Ball(elementRadius);
		std::cout<<"Structuring Element built"<<std::endl;
		
		//Sets Dilate Filter
		typedef itk::GrayscaleErodeImageFilter<TLabelImage, TLabelImage, StructuringElementType> ErodeFilterType;
		typename ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
		erodeFilter->SetInput(labelImage);
		erodeFilter->SetKernel(structElement);
		erodeFilter->Update();
// 		ImageHelper<TLabelImage>::WriteImage(erodeFilter->GetOutput(),"/home/caldairou/Tools/FBrain/test_fbrain/LEI_El_Essai_DilateLabel1_DilateFunction.nii");
		
		return erodeFilter->GetOutput();
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::SubtractImage(typename TLabelImage::Pointer labelImage, typename TLabelImage::Pointer dilateImage)
	{
		//Set Subtract Filter
		typedef itk::SubtractImageFilter<TLabelImage, TLabelImage, TLabelImage> SubtractFilterType;
		typename SubtractFilterType::Pointer subtractFilter = SubtractFilterType::New();
		subtractFilter->SetInput1(dilateImage);
		subtractFilter->SetInput2(labelImage);
		subtractFilter->Update();
		
		return subtractFilter->GetOutput();
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::MaskImage(typename TLabelImage::Pointer diffImage, typename TLabelImage::Pointer segImage)
	{
		//Set Mask Filter
		typedef itk::MaskImageFilter<TLabelImage, TLabelImage, TLabelImage> MaskFilterType;
		typename MaskFilterType::Pointer maskFilter = MaskFilterType::New();
		maskFilter->SetInput(diffImage);
		maskFilter->SetMaskImage(segImage);
		maskFilter->Update();
		
		return maskFilter->GetOutput();
	}
	
	/* ----------------------------------------------Input Acces-------------------------------------------- */
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::SetInputImage(const TInputImage* image)
	{
		SetNthInput(0, const_cast<TInputImage*>(image));
		
		m_Centroids.SetSize(image->GetNumberOfComponentsPerPixel(),3); // 3 because there is three labels
		m_Centroids.Fill(0);
	}
	
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::SetInitialSegmentation(const TLabelImage* initSeg)
	{
		SetNthInput(1, const_cast<TLabelImage*>(initSeg));
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TInputImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::GetInputImage()
	{
		return static_cast< const TInputImage * >
		( this->itk::ProcessObject::GetInput(0) );
	}
	
	template< typename TInputImage, typename TLabelImage>
	typename TLabelImage::Pointer TopologicalKMeans<TInputImage, TLabelImage>::GetInitialSegmentation()
	{
		return static_cast< const TLabelImage * >
		( this->itk::ProcessObject::GetInput(1) );
	}
}

#endif // BTK_TOPOLOGICAL_KMEANS_TXX
