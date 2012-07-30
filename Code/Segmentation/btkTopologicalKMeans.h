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

#ifndef BTK_TOPOLOGICAL_KMEANS_H
#define BTK_TOPOLOGICAL_KMEANS_H

//ITK Includes
#include "itkImageToImageFilter.h"
#include "itkVariableSizeMatrix.h"

//WARNING
//Type TInputImage shall be a vector Image
//Type TLabelImage shall be a scalar Image

namespace btk
{
	/**
	* FCM class for fuzzy classification
	* @author Benoit Caldairou
	*/
	
	template <class TInputImage, class TLabelImage>
	class TopologicalKMeans : public itk::ImageToImageFilter< TInputImage, TLabelImage >
	{
		public :
			/** Standard class typedefs. */
			typedef TopologicalKMeans 	Self;
			typedef itk::ImageToImageFilter< TInputImage, TLabelImage > 	Superclass;
			typedef itk::SmartPointer < Self > 				Pointer;
			/** Centroids Type */
			typedef itk::VariableSizeMatrix < float > 	CentroidsVectorType; //If there is several images
			
			/** Method for creation through the object factory. */
			itkNewMacro(Self);
			
			/** Run-time type information (and related methods). */
			itkTypeMacro(TopologicalKMeans, itk::ImageToImageFilter);
			
			/** Constructor and destructor */
			TopologicalKMeans();
			~TopologicalKMeans(){}
			
			void SetInputImage(const TInputImage* input);
			void SetInitialSegmentation(const TLabelImage* initSeg);
			
		protected:
			
			/** Does the real work. */
			virtual void GenerateData();
			
			/**  Create the Output */
			itk::DataObject::Pointer MakeOutput(unsigned int idx);
			
			/** Get Functions */
			typename TInputImage::Pointer GetInputImage();
			typename TLabelImage::Pointer GetInitialSegmentation();
			
		private:
			TopologicalKMeans(const Self &); //purposely not implemented
			void operator=(const Self &);  //purposely not implemented
			
			void InitialiseCentroids(typename TInputImage::Pointer inputImage, typename TLabelImage::Pointer segImage);
			void RunSegmentation(typename TInputImage::Pointer inputImage, typename TLabelImage::Pointer segImage);
			unsigned int ClassifyBorderVoxel(typename TInputImage::Pointer inputImage, typename TLabelImage::Pointer segImage, typename TLabelImage::Pointer borderImage, typename TLabelImage::PixelType label);
			typename TLabelImage::Pointer GetBorderImage(typename TLabelImage::Pointer segImage, typename TLabelImage::PixelType label, bool erodeOrDilate);
			typename TLabelImage::Pointer GetOneLabel(typename TLabelImage::Pointer segImage, typename TLabelImage::PixelType label);
			typename TLabelImage::Pointer DilateOneLabel(typename TLabelImage::Pointer labelImage);
			typename TLabelImage::Pointer ErodeOneLabel(typename TLabelImage::Pointer labelImage);
			typename TLabelImage::Pointer SubtractImage(typename TLabelImage::Pointer labelImage, typename TLabelImage::Pointer dilateImage);
			typename TLabelImage::Pointer MaskImage(typename TLabelImage::Pointer diffImage, typename TLabelImage::Pointer segImage);
			void CheckBorderVoxel(typename TLabelImage::Pointer borderImage, typename TLabelImage::Pointer segImage, typename TLabelImage::PixelType label);
			
		private :
			CentroidsVectorType m_Centroids;
	};
	
}

#include "btkTopologicalKMeans.txx"

#endif // BTK_TOPOLOGICAL_KMEANS_H
