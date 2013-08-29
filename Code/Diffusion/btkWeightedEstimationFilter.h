/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/08/2013
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
#include "itkNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkAffineTransform.h"

// Local includes
#include "btkDiffusionSequence.h"
#include "btkWeightedEstimationBase.h"
#include "btkDiffusionDataset.h"



namespace btk
{

/**
 * @brief Weighted spherical harmonics decomposition and model estimation (diffusion sequence).
 *
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

        /** Diffusion dataset typedefs. */
        typedef DiffusionDataset                                                DatasetType;
        typedef DatasetType::Pointer                                            DatasetPointer;
        typedef DatasetType::Iterator                                           DataIterator;

         /** Point and index typedef. */
        typedef itk::Point<double,3>                                            Point3DType;
        typedef itk::Index<3>                                                   Index3DType;

        /** Transform typedefs. */
        typedef itk::AffineTransform< double,3 >                                TransformType;

        /** Neighborhood iterator typedefs. */
        typedef itk::NeighborhoodIterator<itk::Image<short,3> >     NeighborhoodIterator;

        /** Type of Spherical Harmonics Diffusion Decomposition */
        typedef WeightedEstimationBase                                          WeightedEstimationType;
        typedef WeightedEstimationType::Pointer                                 WeightedEstimationPointer;

        /** Containers definition */
        typedef btk::DiffusionSequence::GradientTable                           GradientTableType;
        typedef std::vector<float>                                             VectorType;


        itkNewMacro(Self);
        itkTypeMacro(WeightedEstimationFilter,itk::ImageToImageFilter);

        /**
         * @brief Set input diffusion dataset.
         * @param dataset Input diffusion dataset.
         */
        void SetDiffusionDataset(DatasetPointer dataset);

        /**
         * @brief Set Radius.
         * @param Radius The radius of Neighborhood search.
         */
        btkSetMacro(Radius,double);


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

        /** The radius of Neighborhood search */
        double  m_Radius;

         /** Input Diffusion dataset */
        DatasetPointer                    m_Dataset;

         /** Output sequence */
        Self::OutputImageType::Pointer    m_OutputSequence;

        /** Percentage of progression */
        float               m_Percent;

        /** Nb of loop */
        unsigned int        m_NbLoop;

        /** Size of the input sequence */
        DiffusionSequence::SizeType m_Size4D;


};

} // end namespace btk

#endif // BTKWEIGHTEDESTIMATIONFILTER_H
