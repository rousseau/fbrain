/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 24/10/2012
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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
#ifndef BTKLOWTOHIGHRESOLUTIONFILTER_TXX
#define BTKLOWTOHIGHRESOLUTIONFILTER_TXX

#include "btkLowToHighResolutionFilter.hxx"
namespace btk
{
template<class TImage>
LowToHighResolutionFilter<TImage>::LowToHighResolutionFilter():m_UseReference(false),m_ReferenceImageNum(0),
    m_NumberOfImages(0), m_Margin(0.0)
{
    m_InverseTransforms.resize(0);
    m_Interpolators.resize(0);
    m_Transforms.resize(0);
    m_Masks.resize(0);
    m_Images.resize(0);
}
//------------------------------------------------------------------------------------------------
template<class TImage>
void LowToHighResolutionFilter<TImage>::Initialize()throw(itk::ExceptionObject)
{

    if(m_Images.size() == 0 || m_Masks.size() == 0 || m_Transforms.size() == 0)
    {
        btkException("Images, Masks or Transforms not sets");
    }
    else
    {
        m_NumberOfImages = m_Images.size();
    }
    if(!m_UseReference)
    {
        m_ReferenceImage = m_Images[m_ReferenceImageNum];
        m_ReferenceRegion = m_ReferenceImage->GetLargestPossibleRegion().GetSize();
        m_ReferenceMask = m_Masks[m_ReferenceImageNum];
    }
    else
    {
        if(m_ReferenceImage.GetPointer()== NULL || m_ReferenceMask.GetPointer() == NULL)
        {
            btkException("Reference Image or maks not set");
        }
    }

    m_Interpolators.resize(m_NumberOfImages);
    m_InverseTransforms.resize(m_NumberOfImages);

    for(unsigned int i = 0; i< m_NumberOfImages; i++)
    {
        m_Interpolators[i] = Interpolator::New();
        m_Interpolators[i]-> SetInputImage(m_Images[i]);

        m_InverseTransforms[i] = TransformType::New();
        m_InverseTransforms[i] -> SetIdentity();
        m_InverseTransforms[i] -> SetCenter( m_Transforms[i] -> GetCenter() );
        m_InverseTransforms[i] -> SetParameters( m_Transforms[i] -> GetParameters() );

        m_Transforms[i] -> GetInverse( m_InverseTransforms[i] );

    }

    /* resampling matrix */

    typename ImageType::SpacingType fixedSpacing = m_ReferenceImage->GetSpacing();
    typename ImageType::SizeType fixedSize = m_ReferenceImage->GetLargestPossibleRegion().GetSize();
    typename ImageType::PointType fixedOrigin = m_ReferenceImage->GetOrigin();

    typename ImageType::SpacingType resampleSpacing;
    typename ImageType::SizeType resampleSize;

    // Combine masks to redefine the resampling region

    typename ImageType::IndexType indexMin, indexMax, index, targetIndex;
    typename ImageType::PointType point;

    indexMin = m_ReferenceImage-> GetLargestPossibleRegion().GetIndex();

    indexMax[0] = fixedSize[0]-1;
    indexMax[1] = fixedSize[1]-1;
    indexMax[2] = fixedSize[2]-1;

    for (unsigned int i=0; i<m_NumberOfImages; i++)
    {
        if (i != m_ReferenceImageNum || m_UseReference == true)
        {
            Iterator regionIt(m_Images[i],m_Images[i]->GetLargestPossibleRegion());
            for(regionIt.GoToBegin(); !regionIt.IsAtEnd(); ++regionIt )
            {

                index = regionIt.GetIndex();
                m_Images[i] -> TransformIndexToPhysicalPoint(index,point);
                m_ReferenceImage -> TransformPhysicalPointToIndex(point,targetIndex);

                for (unsigned int k=0; k<3; k++)
                {
                    if (targetIndex[k]<indexMin[k])
                        indexMin[k]=targetIndex[k];
                    else
                        if (targetIndex[k]>indexMax[k])
                            indexMax[k]=targetIndex[k];
                }
            }
        }
    }

    resampleSpacing[0] = fixedSpacing[0];
    resampleSpacing[1] = fixedSpacing[0];
    resampleSpacing[2] = fixedSpacing[0];

    resampleSize[0] = floor((indexMax[0]-indexMin[0]+1)*fixedSpacing[0]/resampleSpacing[0] + 0.5);
    resampleSize[1] = floor((indexMax[1]-indexMin[1]+1)*fixedSpacing[1]/resampleSpacing[1] + 0.5);
    resampleSize[2] = floor(((indexMax[2]-indexMin[2]+1)*fixedSpacing[2] + 2*m_Margin)/resampleSpacing[2] + 0.5);

    /* Create high resolution image */

    m_HighResolutionImage = ImageType::New();

    typename ImageType::IndexType start;
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;

    typename ImageType::RegionType region;
    region.SetIndex(start);
    region.SetSize(resampleSize);

    m_HighResolutionImage -> SetRegions( region );
    m_HighResolutionImage -> Allocate();
    m_HighResolutionImage -> SetOrigin( fixedOrigin );
    m_HighResolutionImage -> SetSpacing( resampleSpacing );
    m_HighResolutionImage -> SetDirection( m_ReferenceImage-> GetDirection() );
    m_HighResolutionImage -> FillBuffer( 0 );

    typename ImageType::IndexType newOriginIndex;
    typename ImageType::PointType newOrigin;

    newOriginIndex[0] = floor(indexMin[0]*fixedSpacing[0]/resampleSpacing[0] + 0.5);
    newOriginIndex[1] = floor(indexMin[1]*fixedSpacing[1]/resampleSpacing[1] + 0.5);
    newOriginIndex[2] = floor(indexMin[2]*fixedSpacing[2]/resampleSpacing[2] + 0.5) - floor( m_Margin / resampleSpacing[2] + 0.5);

    m_HighResolutionImage -> TransformIndexToPhysicalPoint(newOriginIndex, newOrigin);
    m_HighResolutionImage -> SetOrigin( newOrigin );




}
//------------------------------------------------------------------------------------------------
template<class TImage>
void LowToHighResolutionFilter<TImage>::Update()throw(itk::ExceptionObject)
{
    try
    {
        this->Initialize();
    }
    catch( itk::ExceptionObject &exp)
    {
        throw exp;
    }

    //Create combination of masks
    typename ImageType::PointType physicalPoint, transformedPoint;
    typename ImageType::IndexType index;

    m_MaskCombination = MaskType::New();
    m_MaskCombination -> SetRegions( m_HighResolutionImage -> GetLargestPossibleRegion() );
    m_MaskCombination -> Allocate();
    m_MaskCombination -> SetSpacing(   m_HighResolutionImage -> GetSpacing() );
    m_MaskCombination -> SetDirection( m_HighResolutionImage -> GetDirection() );
    m_MaskCombination -> SetOrigin(    m_HighResolutionImage -> GetOrigin() );
    m_MaskCombination -> FillBuffer( 0 );

    Iterator It( m_HighResolutionImage, m_HighResolutionImage -> GetLargestPossibleRegion() );
    MaskIterator maskIt( m_MaskCombination, m_MaskCombination -> GetLargestPossibleRegion() );

    typename MaskInterpolator::Pointer imageMaskInterpolator = MaskInterpolator::New();

    for(It.GoToBegin(),maskIt.GoToBegin(); !It.IsAtEnd(); ++It, ++maskIt )
    {
        index = It.GetIndex();
        m_MaskCombination -> TransformIndexToPhysicalPoint( index, physicalPoint );

        for (unsigned int i=0; i < m_NumberOfImages; i++)
        {
            transformedPoint = m_Transforms[i] -> TransformPoint(physicalPoint);
            imageMaskInterpolator -> SetInputImage( m_Masks[i] );

            if ( imageMaskInterpolator -> IsInsideBuffer(transformedPoint) )
            {
                if ( imageMaskInterpolator -> Evaluate(transformedPoint) > 0 )
                {
                    maskIt.Set(1);
                    continue;
                }
            }
        }
    }

    // Average low resulotion images

//    int value;
//    unsigned int counter;

//    Iterator imageIt( m_HighResolutionImage, m_HighResolutionImage->GetLargestPossibleRegion() );

//    for (imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
//    {
//        index = imageIt.GetIndex();

//        m_HighResolutionImage->TransformIndexToPhysicalPoint( index, physicalPoint);

//        value = 0;
//        counter = 0;

//        for (unsigned int i=0; i < m_NumberOfImages; i++)
//        {

//            transformedPoint = m_Transforms[i]->TransformPoint(physicalPoint);

//            if ( m_Interpolators[i]->IsInsideBuffer(transformedPoint) )
//            {
//                value+= m_Interpolators[i]->Evaluate(transformedPoint);
//                counter++;
//            }
//        }
//        if( counter>0 )
//        {
//            imageIt.Set(value/counter);
//        }

//        if(m_MaskCombination -> GetPixel(index) == 0)
//        {
//            imageIt.Set(0);
//        }

//    }

    m_Output = m_HighResolutionImage;


}
//------------------------------------------------------------------------------------------------
template<class TImage>
std::vector<typename LowToHighResolutionFilter<TImage>::ImageType::Pointer>&
LowToHighResolutionFilter<TImage>::ResampleImages() throw(itk::ExceptionObject)
{

    typedef itk::ResampleImageFilter<ImageType,ImageType> ResampleType;

    m_ResampledImages.resize(m_NumberOfImages);

    for(unsigned int i = 0; i< m_NumberOfImages; i++ )
    {
        typename ResampleType::Pointer resampler =  ResampleType::New();
        resampler -> SetTransform( m_Transforms[i] );
        resampler -> SetInput( m_Images[i] );
        resampler -> SetReferenceImage( m_HighResolutionImage );
        resampler -> SetUseReferenceImage( true );
        resampler -> SetDefaultPixelValue( 0 );
        try
        {
        resampler -> Update();
        }
        catch(itk::ExceptionObject & exp)
        {
            throw (exp);
        }

        m_ResampledImages[i] = resampler->GetOutput();
    }

    return m_ResampledImages;





}
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
}

#endif
