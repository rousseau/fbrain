/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/08/2013
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

#ifndef BTKDIFFUSIONDATASET_H
#define BTKDIFFUSIONDATASET_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkPoint.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkExtractImageFilter.h"

// Local includes
#include "btkDiffusionSlice.h"
#include "btkDiffusionSequence.h"
#include "btkS2SSimilarityFilter.h"

namespace btk
{

class DiffusionDataset : public itk::ProcessObject
{
    public:
        typedef DiffusionDataset                                        Self;
        typedef itk::ProcessObject                                      Superclass;
        typedef itk::SmartPointer< Self >                               Pointer;
        typedef itk::SmartPointer< const Self>                          ConstPointer;

         /** Images typedefs. */
        typedef itk::Image<short,3>                                     ImageType;
        typedef ImageType::Pointer                                      ImagePointer;

        typedef DiffusionSequence::SpacingType                          SpacingType;
        typedef DiffusionSequence::PointType                            Point4DType;
        typedef DiffusionSequence::DirectionType                        Direction4DType;
        typedef DiffusionSequence::RegionType                           Region4DType;
        typedef std::vector<DiffusionSlice::Pointer>                    DataVectorType;

         /** Extractor typedefs. */
        typedef itk::ExtractImageFilter<DiffusionSequence,ImageType>    ExtractorType;
        typedef ExtractorType::Pointer                                  ExtractorPointer;

        /** Slice by slice similarity filter typedefs. */
        typedef S2SSimilarityFilter<ImageType>                          SimilarityType;
        typedef SimilarityType::Pointer                                 SimilarityPointer;

        /** Gradient Table typedef */
        typedef DiffusionSequence::GradientTable                        GradientTableType;

        /** Transform vector typedef */
        typedef std::vector<DiffusionSlice::TransformPointer>           TransformVectorType;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        itkTypeMacro(DiffusionDataset, itk::ProcessObject);

        /**
         * @brief Set/Get the input diffusion sequence.
         * @param InputSequence The input diffusion sequence.
         */
        btkSetMacro(InputSequence, DiffusionSequence::ConstPointer);
        btkGetMacro(InputSequence, DiffusionSequence::ConstPointer);


        /**
         * @brief Set Initial slice transforms
         * @param initialTransforms Initial slice transforms
         */
        void SetInitialTransforms(TransformVectorType initialTransforms)
        {
            m_InitialTransforms = initialTransforms;

        }

        /**
         * @brief Remove slice which avec an OutliersStatus = true
         */
        void RemoveOutliers();

        /**
         * @brief Get original sequence origin
         */
        btkGetMacro(Origin,Point4DType);

        /**
         * @brief Get original sequence spacing.
         */
        btkGetMacro(Spacing,SpacingType);

        /**
         * @brief Get original sequence direction.
         */
        btkGetMacro(Direction,Direction4DType);

        /**
         * @brief Get original sequence region.
         */
        btkGetMacro(SequenceRegion,Region4DType);

        /**
         * @brief Get original sequence gradient table.
         */
        btkGetMacro(GradientTable,GradientTableType);

        /**
         * @brief Get original B0.
         */
        btkSetMacro(B0,ImagePointer);
        btkGetMacro(B0,ImagePointer);

        /**
         * @brief Get data size.
         */
        unsigned int GetSize()
        {
            return m_Data.size();
        }


        /**
         * @brief Return begin iterator of data vector.
         * @return begin iterator
         */
        DataVectorType::iterator begin()
        {
            return m_Data.begin();
        }

        /**
         * @brief Return end iterator of data vector.
         * @return end iterator
         */
        DataVectorType::iterator end()
        {
            return m_Data.end();
        }


        /**
         * @brief Initialize parameters.
         */
        void Initialize();


    protected:

        /**
         * @brief Constructor.
         */
        DiffusionDataset();

        /**
         * @brief Destructor.
         */
        virtual ~DiffusionDataset();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:

        /** Original sequence origin */
        Point4DType                     m_Origin;

         /** Original sequence spacing */
        SpacingType                     m_Spacing;

        /** Original sequence direction */
        Direction4DType                 m_Direction;

          /** Original sequence region */
        Region4DType                    m_SequenceRegion;

         /** vector of transforms */
        TransformVectorType             m_InitialTransforms;

        /** Inpute diffusion sequence */
        DiffusionSequence::ConstPointer m_InputSequence;

        /** Original B0 */
        ImagePointer                    m_B0;

         /** Original gradient table */
        GradientTableType               m_GradientTable;

        /** Vector of slices */
        DataVectorType                  m_Data;

        /** size of data = number of slices */
        unsigned int                    m_Size;




};
} // end namespace btk

#endif // BTKDIFFUSIONDATASET_H
