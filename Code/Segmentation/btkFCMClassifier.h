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

#ifndef BTK_FCM_CLASSIFIER_H
#define BTK_FCM_CLASSIFIER_H

//ITK Includes
#include "itkImageToImageFilter.h"
#include <itkVariableLengthVector.h>

//WARNING
//Type TGreyImage and TLabelImage are for scalar images
//Type TFuzzyImage shall be a Vector Image

namespace btk
{
	/**
	* FCM class for fuzzy classification
	* @author Benoit Caldairou
	*/
	
	template <class TGreyImage, class TLabelImage, class TFuzzyImage>
	class FCMClassifier : public itk::ImageToImageFilter< TGreyImage, TLabelImage >
	{
		public :
			/** Standard class typedefs. */
			typedef FCMClassifier 	Self;
			typedef itk::ImageToImageFilter< TGreyImage, TLabelImage > 	Superclass;
			typedef itk::SmartPointer < Self > 				Pointer;
			
			/** Method for creation through the object factory. */
			itkNewMacro(Self);
			
			/** Run-time type information (and related methods). */
			itkTypeMacro(FCMClassifier, itk::ImageToImageFilter);
			
			/** Centroids Type */
			typedef itk::VariableLengthVector < float > 	CentroidsVectorType;
			
			/** Constructor and destructor */
			FCMClassifier();
			~FCMClassifier(){}
			
			/** Acces Functions */
			TLabelImage* GetLabelSegmentation();
			TFuzzyImage* GetFuzzyMaps();
			
			void SetGreyImage(const TGreyImage* input);
			void SetMaskImage(const TLabelImage* mask);
			
			unsigned int GetClassNumber() {return m_ClassNumber;}
			void SetClassNumber(unsigned int classNumber) {m_ClassNumber = classNumber;}
			
		protected:
			
			/** Does the real work. */
			virtual void GenerateData();
			
			/**  Create the Output */
			itk::DataObject::Pointer MakeOutput(unsigned int idx);
			
			/** Get Functions */
			typename TGreyImage::Pointer GetGreyImage();
			typename TLabelImage::Pointer GetMaskImage();
			
		private:
			FCMClassifier(const Self &); //purposely not implemented
			void operator=(const Self &);  //purposely not implemented
			void InitialiseCentroids(typename TGreyImage::Pointer inputImage, typename TLabelImage::Pointer maskImage);
			void FCMOptimisation(typename TGreyImage::Pointer inputImage, typename TLabelImage::Pointer maskImage, typename TFuzzyImage::Pointer fuzzyMaps);
			void MakeLabelImage(typename TLabelImage::Pointer maskImage, typename TLabelImage::Pointer labelImage, typename TFuzzyImage::Pointer fuzzyMaps);
		
		private :
			
			unsigned int m_ClassNumber;
			CentroidsVectorType m_Centroids;
// 			typename TGreyImage::Pointer m_Input;
// 			typename TLabelImage::Pointer m_Mask;
// 			typename TLabelImage::Pointer m_LabelSegmentation;
// 			typename TFuzzyImage::Pointer m_FuzzySegmentation; //Number of fuzzy maps depends on the number of reasearched clusters
		
	};
	
}

#include "btkFCMClassifier.txx"

#endif // BTK_FCM_CLASSIFIER_H
