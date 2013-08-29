/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 25/06/2013
  Author(s):Frederic Champ(champ(at)unistra.fr)
  
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

#ifndef BTKRESOLUTIONSAMPLERFILTER_H
#define BTKRESOLUTIONSAMPLERFILTER_H



// ITK includes
#include "itkSmartPointer.h"
#include "itkResampleImageFilter.h"
#include "itkAffineTransform.h"

// Local includes
#include "btkMacro.h"


namespace btk
{

template <class TImage>
class ResolutionSamplerFilter : public itk::ResampleImageFilter< TImage,TImage >
{
    public:
        /** Standard class typedefs. */
        typedef ResolutionSamplerFilter                                 Self;
        typedef itk::ResampleImageFilter< TImage,TImage>                Superclass;
        typedef itk::SmartPointer< Self >                               Pointer;

        /** Image typedefs */
        typedef typename TImage::Pointer                                ImagePointer;
        typedef typename TImage::ConstPointer                           ImageConstPointer;
        typedef typename TImage::IndexType                              ImageIndexType;
        typedef typename TImage::PointType                              ImagePointType;
        typedef typename TImage::SizeType                               ImageSizeType;
        typedef typename TImage::SpacingType                            ImageSpacingType;
        typedef typename TImage::DirectionType                          ImageDirectionType;
        typedef typename TImage::RegionType                             ImageRegionType;

        /** Transform typedefs. */
        typedef typename itk::AffineTransform<double,3>                  TransformType;
        typedef typename TransformType::Pointer                          TransformPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        itkTypeMacro(ResolutionSamplerFilter, itk::ResampleImageFilter);

        /**
         * @brief Set/Get Input Image.
         * @param InputImage The image which has to be resampled.
         */
        btkSetMacro(InputImage,ImageConstPointer);
        btkGetMacro(InputImage,ImageConstPointer);

        /**
         * @brief Set/Get the new resolution in mm per voxel's sides.
         * @param Resolution The new resolution in mm.
         */
        btkSetMacro(Resolution,unsigned int);
        btkGetMacro(Resolution,unsigned int);

        /**
        * @brief Update process.
        */
       virtual void Update();

    protected:

        /**
          * @brief Constructor.
          */
         ResolutionSamplerFilter();

         /**
          * @brief Destructor
          */
         virtual ~ResolutionSamplerFilter();


    private:

        /** Resolution */
        unsigned int  m_Resolution;

        /** Input Image */
        ImageConstPointer m_InputImage;

        /** Output image */
        ImagePointer m_OutputImage;


};

} // end namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkResolutionSamplerFilter.txx"
#endif

#endif // BTKRESOLUTIONSAMPLERFILTER_H
