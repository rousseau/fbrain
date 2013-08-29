/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 28/08/2013
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

#ifndef BTKDWIREGISTRATIONFILTER_H
#define BTKDWIREGISTRATIONFILTER_H

// ITK includes
#include "itkProcessObject.h"
#include "itkImageRegistrationMethod.h"
#include "itkPowellOptimizer.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"

// Local includes
#include "btkDiffusionSlice.h"
#include "btkDiffusionDataset.h"
#include "btkDiffusionSequence.h"
#include "btkResolutionSamplerFilter.h"
#include "btkAffineTransform.h"
#include "btkWeightedSumOfImagesFilter.h"

namespace btk
{

/**
  * @brief Slice by slice registration for DWIs
  *
  * @author Frederic Champ
  * @ingroup Registration
  *
  */

class DwiReconstructionFilter: public itk::ProcessObject
{
    public:
        /** Standard class typedef. */
        typedef DwiReconstructionFilter                                       Self;
        typedef ProcessObject                                           Superclass;
        typedef itk::SmartPointer<Self>                                 Pointer;
        typedef itk::SmartPointer<const Self>                           ConstPointer;

        /** Image typedefs. */
        typedef itk::Image<short,3>                                     ImageType;
        typedef ImageType::Pointer                                      ImagePointer;
        typedef ImageType::ConstPointer                                 ImageConstPointer;

        /** Diffusion sequence typedefs. */
        typedef DiffusionSequence                                       SequenceType;
        typedef SequenceType::Pointer                                   SequencePointer;
        typedef SequenceType::ConstPointer                              SequenceConstPointer;

        /** Diffusion dataset typedefs. */
        typedef DiffusionDataset                                        DatasetType;
        typedef DiffusionDataset::Pointer                               DatasetPointer;


        /** Registration method typedefs. */
        typedef itk::ImageRegistrationMethod< ImageType,ImageType >     RegistrationType;
        typedef RegistrationType::Pointer                               RegistrationPointer;
        typedef RegistrationType::ParametersType                        RegistrationParameters;

        /** Powell optimizer typedefs. */
        typedef itk::PowellOptimizer                                    PowellOptimizerType;
        typedef PowellOptimizerType::Pointer                            PowellOptimizerPointer;

        /** gradient direction optimizer typedefs. */
        typedef itk::RegularStepGradientDescentOptimizer                GradOptimizerType;
        typedef GradOptimizerType::Pointer                              GradOptimizerPointer;

        /** Mattes metric typedefs. */
        typedef itk::MattesMutualInformationImageToImageMetric<
                                                         ImageType,
                                                         ImageType>     MattesMetricType;
        typedef MattesMetricType::Pointer                               MattesMetricPointer;

        /** Normalized correlation metric typedefs. */
        typedef itk::NormalizedCorrelationImageToImageMetric<
                                                          ImageType,
                                                          ImageType >   NormalizedMetricType;
        typedef NormalizedMetricType::Pointer                           NormalizedMetricPointer;

        /** Interpolator typedefs. */
        typedef itk::BSplineInterpolateImageFunction<
                                               ImageType,
                                               double,
                                               double >                 InterpolatorType;
        typedef InterpolatorType::Pointer                               InterPolatorPointer;

        /** Transform typedefs. */
        typedef AffineTransform< double,3 >                             TransformType;
        typedef TransformType::Pointer                                  TransformPointer;

        /** Extractor typedefs. */
        typedef itk::ExtractImageFilter< SequenceType,ImageType >       ExtractorType;
        typedef ExtractorType::Pointer                                  ExtractorPointer;


        /** Joiner Series typedefs. */
        typedef itk::JoinSeriesImageFilter< ImageType,SequenceType>     JoinerType;
        typedef JoinerType::Pointer                                     JoinerPointer;

        /** resampler typedefs. */
        typedef itk::ResampleImageFilter< ImageType,ImageType>          ResamplerType;

        /** resolution filter typedefs. */
        typedef ResolutionSamplerFilter<ImageType>                      ResolutionFilterType;
        typedef ResolutionFilterType::Pointer                           ResolutionPointer;

        /** Image vector typedefs. */
        typedef std::vector< ImagePointer >                             ImageVector;

        /** Image vector typedefs. */
        typedef WeightedSumOfImagesFilter<
                                          itk::Image<double,3>,
                                          double >                      WeightedSumType;
        typedef WeightedSumType::Pointer                                WeightedSumPointer;


        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        itkTypeMacro(DwiReconstructionFilter,ProcessObject);

        /**
         * @brief Set/Get number of iterations.
         * @param NumberOfIterations Number of iterations.
         */
        btkSetMacro(NumberOfIterations,unsigned int);
        btkGetMacro(NumberOfIterations,unsigned int);

        /**
         * @brief Set/Get reference image.
         * @param ReferenceImage Pointer of reference image.
         */
        btkSetMacro(ReferenceImage,ImageConstPointer);
        btkGetMacro(ReferenceImage,ImageConstPointer);

        /**
         * @brief Set/Get input sequence.
         * @param InputSequence Pointer of input sequence.
         */
        btkSetMacro(InputSequence,SequenceConstPointer);
        btkGetMacro(InputSequence,SequenceConstPointer);

        /**
         * @brief Set/Get output sequence.
         * @param OutputSequence Pointer of output sequence.
         */
        btkSetMacro(OutputSequence,SequencePointer);
        btkGetMacro(OutputSequence,SequencePointer);

        /**
         * @brief Set volume registration boolean.
         * @param VolumeRegistration.
         */
        btkSetMacro(VolumeRegistration,bool);


        /**
         * @brief Set slice registration boolean.
         * @param SliceRegistration.
         */
        btkSetMacro(SliceRegistration,bool);


        /**
         * @brief Get output dataset.
         * @return Dataset  Diffusion dataset.
         */
        DatasetPointer GetOutput()
        {
            return m_Dataset;
        }


        /**
         * @brief Update process
         */
        virtual void Update();


    protected:
        /**
         * @brief Constructor.
         */
        DwiReconstructionFilter();
        /**
         * @brief Destructor.
         */
        virtual ~DwiReconstructionFilter();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /** Intialize method */
        void Initialize();

    private:

        /** Volume registration boolean */
        bool              m_VolumeRegistration;

        /** Slice registration boolean */
        bool              m_SliceRegistration;

        /** Number of iterations */
        unsigned int      m_NumberOfIterations;

        /** Referecne image */
        ImageConstPointer m_ReferenceImage;

        /** Input sequence */
        SequenceConstPointer   m_InputSequence;

        /** Output sequence */
        SequencePointer   m_OutputSequence;

        /** Dataset */
        DatasetPointer    m_Dataset;

        /** Vecotr of gradient images */
        ImageVector       m_GradientImages;

        /** Input 4D region */
        SequenceType::RegionType m_Region4D;

        /** Input 4D size */
        SequenceType::SizeType   m_Size4D;


};
} // end namespace btk

#endif // BTKDWIREGISTRATIONFILTER_H
