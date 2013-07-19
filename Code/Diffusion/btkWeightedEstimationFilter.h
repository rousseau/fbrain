/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/06/2013
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

#ifndef BTKWEIGHTEDESTIMATIONFILTER_H
#define BTKWEIGHTEDESTIMATIONFILTER_H

// ITK includes
#include "itkImageToImageFilter.h"
#include "itkVariableSizeMatrix.h"
#include "itkGradientImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkEllipsoidInteriorExteriorSpatialFunction.h"

// Local includes
#include "btkDiffusionSequence.h"
#include "btkWeightedEstimationBase.h"


namespace btk
{
/**
 * @brief Weighted spherical harmonics decomposition and model estimation (diffusion sequence).
 * @author Frederic Champ
 * @ingroup Diffusion
 */

class WeightedEstimationFilter : public itk::ImageToImageFilter< DiffusionSequence, DiffusionSequence >
{
    public:
        typedef WeightedEstimationFilter                                        Self;
        typedef itk::ImageToImageFilter< DiffusionSequence, DiffusionSequence > Superclass;
        typedef itk::SmartPointer< Self >                                       Pointer;
        typedef itk::SmartPointer< const Self >                                 ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(WeightedEstimationFilter,itk::ImageToImageFilter);

        /** DiffusionSequence types definitions */
        typedef DiffusionSequence                                               TSequence;
        typedef TSequence::Pointer                                              SequencePointer;
        typedef TSequence::ConstPointer                                         SequenceConstPointer;
        typedef TSequence::IndexType                                            SequenceIndexType;

        /** Gradient image types definitions */
        typedef itk::GradientImageFilter< itk::Image<short,4>, float>           GradientFilterType;
        typedef GradientFilterType::Pointer                                     GradientFilterPointer;
        typedef itk::Image< itk::CovariantVector< float, 4 >, 4 >               TGradient;
        typedef TGradient::Pointer                                              GradientPointer;
        typedef itk::ImageRegionConstIteratorWithIndex<TGradient>               GradientIterator;
        typedef TGradient::PointType                                            GradPoint;

        /** Ellispoid function types definitions */
        typedef itk::EllipsoidInteriorExteriorSpatialFunction<3>                EllipsoidFunctionType;
        typedef EllipsoidFunctionType::Pointer                                  EllipsoidFunctionPointer;
        typedef EllipsoidFunctionType::InputType                                EllipsoidFunctionVector;

        /** Type of Spherical Harmonics Diffusion Decomposition */
        typedef WeightedEstimationBase                                          WeightedEstimationType;
        typedef WeightedEstimationType::Pointer                                 WeightedEstimationPointer;

        /** Containers definition */
        typedef btk::DiffusionSequence::GradientTable                           GradientTableType;
        typedef std::vector<SequenceIndexType>                                  IndexVectorType;
        typedef std::vector< std::vector<bool> >                                boolIndexesVector;
        typedef std::vector<float>                                              VectorType;


        /**
         * @brief Set/Get booleans vector. lines = number of slices, columuns = number of volumes
         * for each slice of each volumes you have a boolean: true if the slice is an outliers false otherwise.
         */
        btkSetMacro(BoolIndexes, boolIndexesVector);
        btkGetMacro(BoolIndexes, boolIndexesVector);

        /**
         * @brief Set/Get Sigma in case you use the WSH method (default = 1.0).
         */
        btkSetMacro(Sigma, float);
        btkGetMacro(Sigma, float);

        /**
         * @brief Set/Get ellipsoid size factor in case you use the WSH method (default = 1.0).
         */
        btkSetMacro(EllipsoidSize, float);
        btkGetMacro(EllipsoidSize, float);

        /** Set/Get the verbose mod */
        btkSetMacro(VerboseMod, bool);
        btkGetMacro(VerboseMod, bool);

        /**
         * @brief Set input diffusion sequence.
         * @param InputSequence Input diffusion sequence.
         */
        void SetInputSequence(SequenceConstPointer InputSequence);


    protected:
        /**
         * @brief Constructor.
         */
        WeightedEstimationFilter();

        /**
         * @brief Destructor.
         */
        virtual ~WeightedEstimationFilter();

        /**
         * @brief Pre-processing.
         */
        virtual void BeforeThreadedGenerateData();

        /**
         * @brief Allocate correctly the output image.
         */
        virtual void AllocateOutputs();

        /**
         * @brief The execute method for each thread.
         */
        virtual void ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId);

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:

        /** verbose mod */
        bool                m_VerboseMod;

        /** Sigma */
        float               m_Sigma;

        /** Ellipsoid Size */
        float               m_EllipsoidSize;

        /** Percentage of progression */
        float               m_Percent;

        /** Nb of loop */
        unsigned int        m_NbLoop;

        /** Size of the input sequence */
        TSequence::SizeType m_Size4D;

        /** Outliers Indexes */
        boolIndexesVector   m_BoolIndexes;

        /** Gradient Table */
        GradientTableType   m_GradientTable;

        /**  InputSequence */
        SequencePointer     m_InputSequence;

        /** Gradient Image */
        GradientPointer     m_GradientImage;






};
}// end namespace btk
#endif // BTKWEIGHTEDDECOMPOSITIONFILTER_H
