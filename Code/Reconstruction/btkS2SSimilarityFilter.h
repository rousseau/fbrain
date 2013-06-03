/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date:
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

#ifndef __btkS2SSimilarityFilter_h__
#define __btkS2SSimilarityFilter_h__


// TCLAP includes
#include "tclap/CmdLine.h"

// ITK includes
#include "itkImage.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkJoinSeriesImageFilter.h"

// Local includes
#include "btkFileHelper.h"
#include "btkImageHelper.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkDiffusionSequenceFileHelper.h"
#include "btkMutualInformation.h"




namespace btk
{
/**
 * @brief Calculate slice to slice similarities
 * @author Frederic Champ
 * @ingroup Reconstruction
 */

template <class TImage>
class S2SSimilarityFilter: public ProcessObject
{
    public:
        /** Standard class typedefs. */
        typedef S2SSimilarityFilter                                             Self;
        typedef ProcessObject                                                   Superclass;
        typedef itk::SmartPointer< Self >                                       Pointer;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** 3D Image types definition */
        typedef typename TImage::Pointer                                        ImagePointer;
        typedef typename TImage::ConstPointer                                   ImageConstPointer;
        typedef typename TImage::IndexType                                      ImageIndexType;
        typedef typename TImage::SizeType                                       ImageSizeType;
        typedef typename TImage::RegionType                                     ImageRegionType;


        /** Sequence types definition */
        typedef btk::DiffusionSequence                                          TSequence;
        typedef typename btk::DiffusionSequence::ConstPointer                   SequenceConstPointer;
        typedef typename btk::DiffusionSequence::IndexType                      SequenceIndexType;
        typedef typename btk::DiffusionSequence::SizeType                       SequenceSizeType;
        typedef typename btk::DiffusionSequence::RegionType                     SequenceRegionType;

        /** Image extraction filter definition */
        typedef typename itk::ExtractImageFilter< TSequence,TImage >            SequenceExtractor;
        typedef typename SequenceExtractor::Pointer                             SequenceExtractorPointer;

        /** Slice extraction filter definition */
        typedef typename itk::ExtractImageFilter< TImage,TImage >               ImageExtractor;
        typedef typename ImageExtractor::Pointer                                ImageExtractorPointer;

         /** join image filter definition */
        typedef typename itk::JoinSeriesImageFilter< TImage,TSequence >         ImageSequenceJoiner;
        typedef typename itk::ResampleImageFilter< TImage,TImage >              ImageResampler;
        typedef typename itk::ImageDuplicator< TImage >                         DuplicatorType;

         /** Metric definition */
        typedef typename btk::MutualInformation<TImage>                         MetricType;
        typedef typename MetricType::Pointer                                    MetricPointer;

        /** Values container definition */
        typedef std::vector< float  >                                           S2SVectorType;
        typedef std::vector<std::vector< float  > >                             SimilarityS2SType;



        /**
         * @brief Set verbose mod
         * @param VerboseMod true or false;
         */
        btkSetMacro(VerboseMod, bool);

        /**
         * @brief Set Method used to calculate similarites.
         * @param Method: NonNormalized or Normalized.
         */
        btkSetMacro(Method, std::string);


        /**
         * @brief Set the delimiter between values in the .csv file.
         * @param Delimiter The dilimiter, default ";".
         */
        btkSetMacro(Delimiter, std::string);

        /**
         * @brief Set the slice number for which the similiraty will be calculated
         * @param SliceNumber
         */
        void SetSliceNumber(unsigned int SliceNumber);

        /**
         * @brief Set the volume number for which the similiraty will be calculated
         * @param VolumeNumber
         */
        void SetVolumeNumber(unsigned int VolumeNumber);


        /**
         * @brief Set the reference volume.
         * @param Reference The 3D reference volume.
         */
        btkSetMacro(Reference, ImageConstPointer);

        /**
         * @brief Set input diffusion sequence.
         * @param InputSequence Input diffusion sequence.
         */
        btkSetMacro(InputSequence, SequenceConstPointer);


        /**
         * @brief Write similarities in a .csv file.
         * @param FileName Name of the .csv file.
         */
        void WriteData(std::string FileName);

        /**
         * @brief the matrice of similarities. lines = slices , columns = DWI volumes.
         */
        btkGetMacro(Similarity,SimilarityS2SType);

        /**
         * @brief Get the vector of means of similarites for each slice.
         */
        btkGetMacro(Mean, S2SVectorType);

        /**
         * @brief Get the vector of standard deviations for each slice.
         */
        btkGetMacro(StandardDeviation,S2SVectorType );

        /**
         * @brief Update the process.
         */
        virtual void Update();

    protected:

        /**
         * @brief Constructor.
         */
        S2SSimilarityFilter();

        /**
         * @brief Destructor.
         */
        virtual ~S2SSimilarityFilter();

        /**
         * @brief Allocate structures
         */
        void Initialize();


    private:


        /**
         * @brief Verbose Mod
         */
        bool                    m_VerboseMod;


        /**
        * @brief the slice number.
        */
        unsigned int            m_SliceNumber;

        /**
        * @brief number of slices.
        */
        unsigned int            m_NumberOfSlices;

        /**
        * @brief the image number.
        */
        unsigned int            m_VolumeNumber;

        /**
        * @brief number of images.
        */
        unsigned int            m_NumberOfVolumes;

        /**
         * @brief Method used to calulate similarities
         */
        std::string             m_Method;

        /**
         * @brief delimiter between similarites values.
         */
        std::string             m_Delimiter;

        /**
        * @brief Mean of similarties for each slice.
        */
        S2SVectorType           m_Mean;

        /**
        * @brief Standard deviation for each slice.
        */
        S2SVectorType           m_StandardDeviation;

        /**
        * @brief Similarities container.
        */
        SimilarityS2SType       m_Similarity;

        /**
        * @brief Slice Size.
        */
        SequenceSizeType        m_Slice4DSize;

        /**
         * @brief Slice Size.
         */
        SequenceRegionType      m_SliceRegion;

        /**
         * @brief Reference volume.
         */
        ImageConstPointer       m_Reference;

        /**
         * @brief Input sequence.
         */
        SequenceConstPointer    m_InputSequence;




};

}  // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkS2SSimilarityFilter.txx"
#endif

#endif // btkS2SSimilarityFilter_h

