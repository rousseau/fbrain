/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/06/2013
  Author(s): Frederic Champ (champ(at)unistra.fr)
  
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

#ifndef __btkOutlierCorrectionFilter_h
#define __btkOutlierCorrectionFilter_h


// ITK includes
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkKdTree.h"
#include "itkKdTreeGenerator.h"
#include "itkListSample.h"
#include "itkGradientImageFilter.h"

// Local includes
#include "btkMutualInformation.h"
#include "btkDiffusionSequence.h"
#include "btkFileHelper.h"
#include "btkSphericalHarmonicsDiffusionDecompositionFilter.h"
#include "btkOrientationDiffusionFunctionModel.h"
#include "btkRBFInterpolateImageFunctionS2S.h"
#include "btkS2SSimilarityFilter.h"
#include "btkWeightedEstimationFilter.h"
#include "itkEllipsoidInteriorExteriorSpatialFunction.h"
#include "btkAffineSliceBySliceTransform.h"
#include "btkDiffusionDataset.h"





namespace btk
{
/**
 * @brief Detects the corrupted images and replaces the values of voxels with estimated values.
 *
 * Values can be estimated by SH: Spherical Harmonics, WSH: Weighted Spherical Harmonics, RBF : Radial Basis Function (old).
 *
 * @author Frederic Champ
 * @ingroup Reconstruction
 */

template <class TImage>
class OutlierCorrectionFilter : public ProcessObject
{
    public:
        /** Standard class typedefs. */
        typedef OutlierCorrectionFilter                                 Self;
        typedef ProcessObject                                           Superclass;
        typedef itk::SmartPointer< Self >                               Pointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        typedef typename TImage::Pointer                                ImagePointer;
        typedef typename TImage::IndexType                              ImageIndexType;
        typedef typename TImage::SizeType                               ImageSizeType;
        typedef typename TImage::RegionType                             ImageRegionType;
        typedef ImageRegionIteratorWithIndex< TImage >                  ImageIterator;

        typedef itk::Image< itk::CovariantVector< float, 4 >, 4 >       TGradient;

        typedef btk::DiffusionSequence                                  TSequence;
        typedef typename TSequence::PointType                           SequencePointType;
        typedef typename TSequence::IndexType                           SequenceIndexType;
        typedef typename TSequence::RegionType                          SequenceRegionType;
        typedef typename TSequence::SizeType                            SequenceSizeType;
        typedef typename TSequence::Pointer                             SequencePointer;
        typedef typename TSequence::ConstPointer                        SequenceConstPointer;
        typedef ImageRegionIteratorWithIndex< TSequence >               SequenceIterator;

        /** Diffusion dataset typedefs. */
        typedef DiffusionDataset                                        DatasetType;
        typedef DiffusionDataset::Pointer                               DatasetPointer;

        typedef itk::Point<float,3u>                                    ModelPoint;
        typedef itk::ContinuousIndex<float,3u>                          IndexContinuous;

        /** Type of outliers indexes ntainer */
        typedef std::vector<SequenceIndexType>                          IndexVectorType;

        /** Type of Image extractor */
        typedef itk::ExtractImageFilter< TSequence,TImage >             ImageExtractor;
        typedef typename ImageExtractor::Pointer                        ImageExtractorPointer;

        typedef itk::JoinSeriesImageFilter< TImage, TSequence >         JoinSequenceFilter;
        typedef typename JoinSequenceFilter::Pointer                    JoinSequencePointer;

        /** Type of image joiner */
        typedef itk::JoinImageFilter< TImage, TSequence >               ImageJoiner;

        /** vector of similarities calculated */
        typedef std::vector< float >                                    S2SVectorType;
        typedef std::vector< S2SVectorType >                            SimilarityS2SType;

        /** Type of transform for the metric */
        typedef itk::AffineTransform<double>                            TransformType;
        typedef TransformType::Pointer                                  TransformPointer;

        /** Type of Sequence duplicator */
        typedef itk::ImageDuplicator<TSequence>                         DuplicatorType;
        typedef typename DuplicatorType::Pointer                        DuplicatorPointer;

        /** Type of Spherical Harmonics Diffusion Decomposition */
        typedef SphericalHarmonicsDiffusionDecompositionFilter          DecompositionType;
        typedef typename DecompositionType::Pointer                     DecompositionPointer;

        /** Type of Spherical Harmonics Diffusion Decomposition */
        typedef WeightedEstimationBase  WeightedDecompositionType;
        typedef typename WeightedDecompositionType::Pointer             WeightedDecompositionPointer;

        /** Type of signal model */
        typedef OrientationDiffusionFunctionModel                       ModelType;
        typedef typename    ModelType::Pointer                          ModelPointer;

        /** Type of RBF interpolator */
        typedef RBFInterpolateImageFunctionS2S<TSequence,double>        RBFInterpolatorType;
        typedef RBFInterpolatorType::Pointer                            RBFInterpolatorPointer;

        /** Type of S2SSimilarity */
        typedef S2SSimilarityFilter<TImage>                             S2SSimilarityType;
        typedef typename S2SSimilarityType::Pointer                     S2SSimilarityPointer;

        /** Type of WeightedEstimationFilter */
        typedef typename WeightedEstimationFilter::Pointer              WeighteEstimationFilterPointer;

         /** Type of bool indexes */
        typedef std::vector< std::vector<bool> >                        boolIndexesVector;


        /**
         * @brief Set input diffusion sequence.
         * @param InputSequence Input diffusion sequence.
         */

       // btkSetMacro(InputSequence, SequenceConstPointer);

        void SetInputSequence(SequenceConstPointer InputSequence)
        {
            m_InputSequence = InputSequence;
            this->Initialize();
        }


        /**
         * @brief Set verbose mod
         * @param VerboseMod true or false;
         */
        btkSetMacro(VerboseMod, bool);


        /**
         * @brief Set the r_gra parameter for the RBF interpolation
         * @param Rgra
         */
        btkSetMacro(Rgra,double);

        /**
         * @brief Set tolerance factor
         * @param Factor Tolerance factor compared to the standard deviation
         */
        btkSetMacro(Factor, unsigned int);

        /**
         * @brief Set number of iterations for slice correction.
         * @param NumberOfIterations Number of iteration.
         */
        btkSetMacro(NumberOfIterations, unsigned int);

        /**
         * @brief Set method used to estimate signal model
         * @param Method SH : Spherical Harmonic, WSH : weighted spherical harmonic or RBF : RBF interpolation ( RBF treat only outlier slices);
         */
        btkSetMacro(Method, std::string);

        /**
         * @brief Get a vector with outliers indexes
         * @return A vector of 4D indexes
         */
        btkGetMacro(OutliersIndexes, IndexVectorType);

        /**
         * @brief Set/Get Ellipsoid size
         */
        btkSetMacro(Radius, double);
        btkGetMacro(Radius, double);


        /**
         * @brief Get output
         * @return The corrected sequence
         */
        SequencePointer GetOutput()
        {
            return m_CorrectedSequence;
        }

        /**
         * @brief Write similarities
         * @param FileName The path of the output file
         */
        void WriteSimilarities(std::string SimilaritiesFileName)
        {
            m_SimilarityFileName = SimilaritiesFileName;
        }


        boolIndexesVector GetBoolIndexes();

        /**
         * @brief Update the process.
         */
        virtual void Update() throw (ExceptionObject);

    protected:
        /**
         * @brief Constructor.
         */
        OutlierCorrectionFilter();

        /**
         * @brief Destructor.
         */
        virtual ~OutlierCorrectionFilter();

        /**
         * @brief Initialize process
         */
        void Initialize();

        /**
         * @brief Detect outliers
         * @return A vectors of outliers indexes;
         */
        void GetOutliers() throw (ExceptionObject);

        /**
         * @brief Correct outliers with values inferred from the spherical harmonics.
         * @param InputSequence Sequence to be corrected
         * @return Corrected DW images vector;
         */
        void CorrectOutliers(SequenceConstPointer InputSequence) ;

        SequencePointer ExtractOutlier(SequenceConstPointer InputSequence, std::vector<unsigned int> corruptedVolume, unsigned int Slice);


    private:

        /** Radius of neighbor search */
        double m_Radius;

        /**  Verbose Mod */
        bool                    m_VerboseMod;

        /** the r_gra parameter for RBF method */
        double                  m_Rgra;

        /** Delimiter in the outliers indexes file */
        std::string             m_Delimiter;

        /** The path of the outliers indexes file */
        std::string             m_SimilarityFileName;

        /** The path of the outliers indexes file */
        std::string             m_OutliersFileName;

        /** Method used : SH (Shperical harmonics) ,WSH (weighted spherical harmonics) or RBF (radial basis function) */
        std::string             m_Method;

        /** Size of 4D sequence */
        SequenceSizeType        m_Size4D;

        /** Tolrance factor compared to the standard deviation. */
        unsigned int            m_Factor;

        /** Number of iterations for slice correction. */
        unsigned int            m_NumberOfIterations;

        /** Outliers Indexes vector */
        IndexVectorType         m_OutliersIndexes;

        /** Indexes boolean vector (false if indexes are considered as outliers) */
        boolIndexesVector       m_BoolIndexes;

        /**  Input sequence. */
        SequenceConstPointer    m_InputSequence;

        /** Corrected Sequence */
        SequencePointer         m_CorrectedSequence;




};

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkOutlierCorrectionFilter.txx"
#endif

#endif // BTKOUTLIERCORRECTION_H

