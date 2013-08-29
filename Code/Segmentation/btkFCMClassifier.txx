/*==========================================================================

© Université de Strasbourg - Centre National de la Recherche Scientifique

Date: 19/04/2012
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

#ifndef BTK_FCM_CLASSIFIER_TXX
#define BTK_FCM_CLASSIFIER_TXX

#include "btkFCMClassifier.h"
#include "btkImageHelper.h"

#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkMinimumMaximumImageCalculator.h"

#include <cmath>

namespace btk
{
	/* --------------------------Constructor------------------------------------------------- */
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::FCMClassifier()
	{
		//Two Outputs : Fuzzy Maps (as a vectorImage) and LabelSegmentation
		this->SetNumberOfRequiredOutputs(2);
		//Two Intputs : Original Image and Mask Image
		this->SetNumberOfRequiredInputs(2);
		
		this->SetNthOutput( 0, this->MakeOutput(0) );
		this->SetNthOutput( 1, this->MakeOutput(1) );
		
		m_ClassNumber=2; //Default Value to avoid a crash
	}
	
	/* --------------------------------------FCM running functions------------------------------------- */
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	void FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::GenerateData()
	{
		//Get Inputs
		typename TGreyImage::Pointer inputImage = this->GetGreyImage();
		typename TLabelImage::Pointer maskImage = this->GetMaskImage();
		
		// Setup output 1
		typename TLabelImage::Pointer labelSegmentation = this->GetLabelSegmentation();
		labelSegmentation->SetRegions(inputImage->GetLargestPossibleRegion());
		labelSegmentation->Allocate();
		
		// Setup output 2
		typename TFuzzyImage::Pointer fuzzyMaps = this->GetFuzzyMaps();
		fuzzyMaps->SetVectorLength(m_ClassNumber);
		fuzzyMaps->SetRegions(inputImage->GetLargestPossibleRegion());
		fuzzyMaps->Allocate();
		
		//Setup Centroids Vector
		m_Centroids.SetSize(m_ClassNumber);
		m_Centroids.Fill(0);
		
		//Initialization of Centroids
		InitialiseCentroids(inputImage, maskImage);
		
		//FCM Optimization
		FCMOptimisation(inputImage, maskImage, fuzzyMaps);
		
		//Build LabelSegmentation
		MakeLabelImage(maskImage, labelSegmentation, fuzzyMaps);
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	void FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::InitialiseCentroids(typename TGreyImage::Pointer inputImage, typename TLabelImage::Pointer maskImage)
	{
		itk::ImageRegionConstIterator<TGreyImage> greyImageIterator(inputImage,inputImage->GetLargestPossibleRegion());
		itk::ImageRegionConstIterator<TLabelImage> maskImageIterator(maskImage,maskImage->GetLargestPossibleRegion());
		
		//Get Minimum and Maximum Value inside the mask
		typedef itk::MinimumMaximumImageCalculator<TGreyImage> MinMaxImageCalculatorType;
		typename MinMaxImageCalculatorType::Pointer minMaxCalculator = MinMaxImageCalculatorType::New();
		minMaxCalculator->SetImage(inputImage);
		minMaxCalculator->Compute();
		
		float max = (float) minMaxCalculator->GetMinimum();
		float min = (float) minMaxCalculator->GetMaximum();
		
		for(greyImageIterator.GoToBegin(), maskImageIterator.GoToBegin(); !greyImageIterator.IsAtEnd(); ++greyImageIterator, ++maskImageIterator)
		{
			if(maskImageIterator.Get() > 0)
			{
				if(greyImageIterator.Get() > max)
				{
					max = (float) greyImageIterator.Get();
				}
				if(greyImageIterator.Get() < min)
				{
					min = (float) greyImageIterator.Get();
				}
			}
		}
		
		//Initialisation
		m_Centroids[0] = max;
		m_Centroids[m_ClassNumber] = min;
		for(unsigned int i=1; i<m_ClassNumber-1; i++)
		{
			m_Centroids[i] = m_Centroids[0] + i*(max-min)/(m_ClassNumber-1);
		}
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	void FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::MakeLabelImage(typename TLabelImage::Pointer maskImage, typename TLabelImage::Pointer labelImage, typename TFuzzyImage::Pointer fuzzyMaps)
	{
		itk::ImageRegionConstIterator<TLabelImage> maskImageIterator(maskImage, maskImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<TLabelImage> labelImageIterator(labelImage, labelImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<TFuzzyImage> fuzzyImageIterator(fuzzyMaps, fuzzyMaps->GetLargestPossibleRegion());
		
		//Maximum fuzzy value in fuzzyMaps settles the sharp label
		for(maskImageIterator.GoToBegin(), labelImageIterator.GoToBegin(), fuzzyImageIterator.GoToBegin(); !maskImageIterator.IsAtEnd(); ++maskImageIterator, ++labelImageIterator, ++fuzzyImageIterator)
		{
			if(maskImageIterator.Get() > 0)
			{
				float max = 0;
				for(unsigned int i=0; i<m_ClassNumber; i++)
				{
					if(fuzzyImageIterator.Get()[i] > max)
					{
						max = fuzzyImageIterator.Get()[i];
						labelImageIterator.Set(i+1);
					}
				}
			}
			else
			{
				labelImageIterator.Set(0);
			}
		}
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	void FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::FCMOptimisation(typename TGreyImage::Pointer inputImage, typename TLabelImage::Pointer maskImage, typename TFuzzyImage::Pointer fuzzyMaps)
	{
		typedef typename TFuzzyImage::PixelType FuzzyPixelType;
		itk::ImageRegionConstIterator<TGreyImage> greyImageIterator(inputImage, inputImage->GetLargestPossibleRegion());
		itk::ImageRegionConstIterator<TLabelImage> maskImageIterator(maskImage, maskImage->GetLargestPossibleRegion());
		itk::ImageRegionIterator<TFuzzyImage> fuzzyImageIterator(fuzzyMaps, fuzzyMaps->GetLargestPossibleRegion());
		
		float J=1;
		float J_Old=0;
		
		while(fabs(J-J_Old)/J > 0.01)
		{
			J_Old = J;
			
			//Compute Fuzzy Maps
			for(greyImageIterator.GoToBegin(), maskImageIterator.GoToBegin(), fuzzyImageIterator.GoToBegin(); !greyImageIterator.IsAtEnd(); ++greyImageIterator, ++maskImageIterator, ++fuzzyImageIterator)
			{
				if(maskImageIterator.Get() > 0)
				{
					FuzzyPixelType currentFuzzyPixel(m_ClassNumber);
					float total = 0;
					bool isAtOne = 0;
					
					for(unsigned int i=0; i<m_ClassNumber; i++)
					{
						//In case centroids and greyvalue are the same => impossible to divide by zero
						if((float)greyImageIterator.Get() - m_Centroids[i] != 0)
						{
							currentFuzzyPixel[i] = 1/pow((float)greyImageIterator.Get() - m_Centroids[i], 2.0);
							total += currentFuzzyPixel[i];
						}
						else
						{
							currentFuzzyPixel[i] = -1.0;
							isAtOne = 1;
							break;
						}
					}
					
					if(isAtOne != 1)
					{
						for(unsigned int i=0; i<m_ClassNumber; i++)
						{
							currentFuzzyPixel[i] = currentFuzzyPixel[i]/total;
						}
					}
					else
					{
						//If the current grey value is equal to one centroids
						for(unsigned int i=0; i<m_ClassNumber; i++)
						{
							if(currentFuzzyPixel[i] == -1.0)
								currentFuzzyPixel[i] = 1.0;
							else
								currentFuzzyPixel[i] = 0.0;
						}
					}
					
					fuzzyImageIterator.Set(currentFuzzyPixel);
				}
			}
			
			
			//Compute Centroids
			itk::VariableLengthVector<float> total(m_ClassNumber);
			total.Fill(0);
			m_Centroids.Fill(0);
			
			for(greyImageIterator.GoToBegin(), maskImageIterator.GoToBegin(), fuzzyImageIterator.GoToBegin(); !fuzzyImageIterator.IsAtEnd(); ++greyImageIterator, ++maskImageIterator, ++fuzzyImageIterator)
			{
				if(maskImageIterator.Get() > 0)
				{
					for(unsigned int i=0; i<m_ClassNumber; i++)
					{
						m_Centroids[i] += fuzzyImageIterator.Get()[i] * greyImageIterator.Get();
						total[i] += fuzzyImageIterator.Get()[i];
					}
				}
			}
			
			for(unsigned int i=0; i<m_ClassNumber; i++)
				m_Centroids[i] = m_Centroids[i]/total[i];
			
			
			//Compute J
			J=0;
			
			for(greyImageIterator.GoToBegin(), maskImageIterator.GoToBegin(), fuzzyImageIterator.GoToBegin(); !fuzzyImageIterator.IsAtEnd(); ++greyImageIterator, ++maskImageIterator, ++fuzzyImageIterator)
			{
				if(maskImageIterator.Get() > 0)
				{
					for(unsigned int i=0; i<m_ClassNumber; i++)
					{
						J += pow((float)greyImageIterator.Get() - m_Centroids[i], 2.0) * fuzzyImageIterator.Get()[i];
					}
				}
			}
		}
	}
	
	/* -----------------------------------------Output Acces ---------------------------------------*/
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	itk::DataObject::Pointer FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::MakeOutput(unsigned int idx)
	{
		itk::DataObject::Pointer output;
		
		switch ( idx )
		{
			case 0:
				output = ( TLabelImage::New() ).GetPointer();
				break;
			case 1:
				output = ( TFuzzyImage::New() ).GetPointer();
				break;
			default:
				std::cerr << "No output " << idx << std::endl;
				output = NULL;
				break;
		}
		return output.GetPointer();
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	TLabelImage* FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::GetLabelSegmentation()
	{
		return dynamic_cast< TLabelImage * >(this->itk::ProcessObject::GetOutput(0) );
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	TFuzzyImage* FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::GetFuzzyMaps()
	{
		return dynamic_cast< TFuzzyImage * >(this->itk::ProcessObject::GetOutput(1) );
	}
	
	/* ----------------------------------------------Input Acces-------------------------------------------- */
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	void FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::SetGreyImage(const TGreyImage* image)
	{
        this->SetNthInput(0, const_cast<TGreyImage*>(image));
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	void FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::SetMaskImage(const TLabelImage* mask)
	{
        this->SetNthInput(1, const_cast<TLabelImage*>(mask));
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	typename TGreyImage::Pointer FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::GetGreyImage()
	{
        return static_cast< TGreyImage * >
		( this->itk::ProcessObject::GetInput(0) );
	}
	
	template< typename TGreyImage, typename TLabelImage, typename TFuzzyImage>
	typename TLabelImage::Pointer FCMClassifier<TGreyImage, TLabelImage, TFuzzyImage>::GetMaskImage()
	{
        return static_cast< TLabelImage * >
		( this->itk::ProcessObject::GetInput(1) );
	}
}
#endif //BTK_FCM_CLASSIFIER_TXX
