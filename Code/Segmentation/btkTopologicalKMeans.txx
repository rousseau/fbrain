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
	
	/* --------------------------------------FCM running functions------------------------------------- */
	template< typename TInputImage, typename TLabelImage>
	void TopologicalKMeans<TInputImage, TLabelImage>::GenerateData()
	{
		//Get Inputs
		typename TInputImage::Pointer inputImage = this->GetInputImage();
		typename TLabelImage::Pointer initSeg = this->GetInitialSegmentation();
		
		// Setup output
		typename TLabelImage::Pointer finalSeg = this->GetOutput();
		typename itk::ImageDuplicator<TLabelImage>::Pointer duplicator = itk::ImageDuplicator<TLabelImage>::New();
		duplicator->SetInputImage(initSeg);
		finalSeg = duplicator->GetOutput();
// 		finalSeg->SetRegions(initSeg->GetLargestPossibleRegion());
// 		finalSeg->Allocate();
		
		itk::ImageRegionConstIterator<TInputImage> greyImageIterator(inputImage, inputImage->GetLargestPossibleRegion());
		itk::ImageRegionConstIterator<TLabelImage> initSegIterator(initSeg, initSeg->GetLargestPossibleRegion());
		itk::ImageRegionIterator<TLabelImage> finalSegIterator(finalSeg, finalSeg->GetLargestPossibleRegion());
		
		
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
