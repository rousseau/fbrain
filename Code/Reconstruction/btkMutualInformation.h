/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 03/06/2013
  Author(s):Frederic Champ (champ(at)unistra.fr)
  
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

#ifndef __btkMutualInformation_h
#define __btkMutualInformation_h

// ITK includes
#include "itkImageToImageMetric.h"
#include "itkJoinImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkStatisticsAlgorithm.h"
#include "itkRandomVariateGeneratorBase.h"
#include "itkImageIteratorWithIndex.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageDuplicator.h"

// Local Includes
#include "btkMacro.h"





namespace btk
{
using namespace itk;

/**
 * @class MutualInformation
 * @brief Computes the mutual information to compare two images without any transforms.
 *
 * MutualInformation uses the code of ITK_DIR/Examples/Statistics/ImageMutualInformation1.cxx
 *
 * @warning If you use both SetNumberOfBins and SetPercentageOfBins, you will set a percentage of the Number of bins set
 * , not a percentage of the total number of bins.
 *  for example for an image with 255 grey levels.
 *  If you use SetPercentageOfBins(80), you will have a joint histogram with 318 bins.
 *  If you use SetNumberOfBins(80), you will have a joint histogram with 80 bins.
 *  If you use both, you will have a joint histogram with 64 bins.
 *
 * @author Frederic Champ
 * @ingroup Registration
 */

template <class TImage>
class  MutualInformation : public ProcessObject
{
    public:


        typedef MutualInformation                                               Self;
        typedef SmartPointer<Self>                                              Pointer;
        typedef SmartPointer<const Self>                                        ConstPointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        typedef typename TImage::Pointer                                        ImagePointer;
        typedef typename TImage::ConstPointer                                   ImageConstPointer;
        typedef typename TImage::PointType                                      PointType;
        typedef typename TImage::IndexType                                      IndexType;
        typedef typename TImage::SizeType                                       SizeType;
        typedef typename TImage::RegionType                                     RegionType;
        typedef typename TImage::PixelType                                      PixelType;
        /** Type of the sequence iterator. */
        typedef ImageRegionConstIteratorWithIndex< TImage >                     ConstIteratorType;
        typedef ImageRegionIteratorWithIndex< TImage >                          IteratorType;

        /** Type of Joind Image filter. */
        typedef          JoinImageFilter<TImage,TImage>                         JoinFilterType;
        typedef typename JoinFilterType::Pointer                                JoinFilterPointer;
        typedef typename JoinFilterType::OutputImagePointer                     JoinImagePointer;
        typedef typename JoinFilterType::OutputImageType                        VectorImageType;

        /** Type of Join Histrogram Filter. */
        typedef typename Statistics::ImageToHistogramFilter< VectorImageType >  HistogramFilterType;
        typedef typename HistogramFilterType::Pointer                           HistogramFilterPointer;
        typedef typename HistogramFilterType::HistogramSizeType                 HistogramSizeType;
        typedef typename HistogramFilterType::HistogramMeasurementVectorType    HistogramMeasurementVectorType;


        /** Type of Join Histrogram. */
        typedef typename HistogramFilterType::HistogramType                     HistogramType;
        typedef typename HistogramType::ConstPointer                            HistogramConstPointer;
        typedef typename HistogramType::ConstIterator                           HistogramIterator;


        /**
         * @brief Set/Get verbose mod.
         * @param Boolean.
         */
        itkSetMacro(VerboseMod,bool);
        itkGetMacro(VerboseMod,bool);

        /**
         * @brief Set/Get form of Mutual Information.
         * @param Method (string): NonNormalized, Normalized.
         * Non Normalized : Entropy1 + Entropy2 - JointEntropy.
         * Normalized : Mutual Information divided by the mean entropy of the input images.
         */
        itkSetMacro(Method,std::string);
        itkGetMacro(Method,std::string);

        /**
         * @brief Set/Get number of bins used to compute the joint histogram.
         *  Default = 100 bins
         *  More number of bins is small more the method is fast
         *  but this may introduce a bias
         * @param Number of bins.
         */
        itkSetMacro(NumberOfBins,unsigned int);
        itkGetMacro(NumberOfBins,unsigned int);

        /**
         * @brief Set/Get the percentage of bins used to calculate the mutual information
         *  Default = 100 %
         *  More percentage of bins is small more the method is fast
         *  but this may introduce a bias
         * @param Percentag of bins (between 0 and 100).
         */
        itkSetMacro(PercentageOfBins,unsigned int);
        itkGetMacro(PercentageOfBins,unsigned int);

        /**
         * @brief Set/Get fixed image (reference image)
         * @param Image 3D pointer.
         */
        btkSetMacro(FixedImage, ImagePointer);
        btkGetMacro(FixedImage, ImagePointer);

        /**
         * @brief Set/Get moving image (query image)
         * @param Image 3D pointer.
         */
        btkSetMacro(MovingImage, ImagePointer);
        btkGetMacro(MovingImage, ImagePointer);


        /**
         * @brief Initialize method.
         */
        void Initialize();

        /**
         * @brief Get the similarity value between the two images set in parameters.
         */
        float GetValue()throw (ExceptionObject);


    protected:

        /**
         * @brief Constructor.
         */
        MutualInformation();

        /**
         * @brief Destructor.
         */
        virtual ~MutualInformation();

    private:

        /** Verbose Mod */
        bool                           m_VerboseMod;

        /** Number of bins used to compute join histogram
         *  Default = 100 bins
         */

        unsigned int                   m_NumberOfBins;

        /** Percentage of bins used to calculate the mutual information
         *  Default = 100 %
         */

        unsigned int                   m_PercentageOfBins;

        /** Joint histogram size */
        unsigned int                    m_histoSize;

        /** Minimum bin of joint histogram */
        PixelType                       m_binMin;

        /** Maximum bin of joint histogram */
        PixelType                       m_binMax;

        /** Mutual information method: NonNormalized, Normalized (cf. method definition).*/
        std::string                     m_Method;

        /** Moving image pointer */
        ImagePointer                   m_MovingImage;

        /** Fixed image pointer */
        ImagePointer                   m_FixedImage;

        /** Join Image pointer */
        JoinImagePointer               m_JoinImage;

        /** Minimum bins vector ( image dimension )  of joint histogram */
        HistogramMeasurementVectorType m_BinMinimum;

        /** Maximum bins vector ( image dimension ) of joint histogram */
        HistogramMeasurementVectorType m_BinMaximum;







};

}// end namespace btk
#ifndef ITK_MANUAL_INSTANTIATION
#include "btkMutualInformation.txx"
#endif
#endif // btkMutualInformation_h


