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
#include "itkImageDuplicator.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkVariableLengthVector.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkFlatStructuringElement.h"

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
		
// 		//Setup Iterator over input Image
// 		itk::ImageRegionConstIterator<TInputImage> greyImageIterator(inputImage, inputImage->GetLargestPossibleRegion());
// 		std::cout<<"greyImageIterator Set"<<std::endl;
// 		
// 		//Setup neighborhood radius
// 		typename TLabelImage::SizeType neighborRadius;
// 		neighborRadius[0] = 1;
// 		neighborRadius[1] = 1;
// 		neighborRadius[2] = 1;
// 		std::cout<<"neighbor Radius Set"<<std::endl;
// 		
// 		//Setup neighborhood Iterator
// 		itk::ShapedNeighborhoodIterator<TLabelImage> iterator(neighborRadius, finalSeg, finalSeg->GetLargestPossibleRegion());
// 		std::cout<<"neighborIterator Set"<<std::endl;
		
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
