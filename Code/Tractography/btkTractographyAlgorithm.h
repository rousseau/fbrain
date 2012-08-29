/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 22/08/2012
  Author(s): Julien Pontabry (pontabry@unistra.fr)
  
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

#ifndef BTK_TRACTOGRAPHY_ALGORITHM_H
#define BTK_TRACTOGRAPHY_ALGORITHM_H

// STL includes
#include "vector"

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkProcessObject.h"
#include "itkImage.h"
#include "itkPoint.h"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkAppendPolyData.h"

// Local includes
#include "btkMacro.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionModel.h"

namespace btk
{

/**
 * @brief Representation of an abstract tractography algorithm.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class TractographyAlgorithm : public itk::ProcessObject
{
    public:
        typedef TractographyAlgorithm           Self;
        typedef itk::ProcessObject              Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef itk::Image< short,3 > MaskImage;
        typedef itk::Image< short,3 > LabelImage;
        typedef itk::Point< float,3 > PhysicalPoint;

        itkTypeMacro(TractographyAlgorithm,itk::ProcessObject);

        btkSetMacro(DiffusionSequence, btk::DiffusionSequence::Pointer);
        btkGetMacro(DiffusionSequence, btk::DiffusionSequence::Pointer);

        btkSetMacro(DiffusionModel, btk::DiffusionModel::Pointer);
        btkGetMacro(DiffusionModel, btk::DiffusionModel::Pointer);

        btkSetMacro(Mask, Self::MaskImage::Pointer);
        btkGetMacro(Mask, Self::MaskImage::Pointer);

        btkSetMacro(RegionsOfInterest, Self::LabelImage::Pointer);
        btkGetMacro(RegionsOfInterest, Self::LabelImage::Pointer);

        btkSetMacro(SeedLabels, std::vector< unsigned short >);
        btkGetMacro(SeedLabels, std::vector< unsigned short >);

        btkSetMacro(SeedSpacing, float);
        btkGetMacro(SeedSpacing, float);

        /**
         * @brief Run the algorithm.
         */
        virtual void Update();

        /**
         * @brief Accessor to output estimated fibers
         * @param label Label corresponding to fibers
         * @return Estimated fibers which seeds have label label.
         */
        virtual std::vector< vtkSmartPointer< vtkAppendPolyData > > GetOutputFiber() const;

    protected:
        /**
         * @brief Constructor.
         */
        TractographyAlgorithm();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /**
         * @brief Resample the regions of interest to make a seed per voxels.
         */
        virtual void ResampleLabelImage();

        /**
         * @brief Propagate using the tractography algorithm at a seed point.
         * @param point Seed point.
         */
        virtual void PropagateSeed(Self::PhysicalPoint point) = 0;

    protected:
        /**
         * @brief Diffusion sequence used as measured observations in algorithm.
         */
        btk::DiffusionSequence::Pointer m_DiffusionSequence;

        /**
         * @brief Diffusion modeling used in algorithm.
         */
        btk::DiffusionModel::Pointer m_DiffusionModel;

        /**
         * @brief Mask image determining where the tractography algorithm can process.
         */
        Self::MaskImage::Pointer m_Mask;

        /**
         * @brief Space memory for current fiber under construction.
         */
        vtkSmartPointer< vtkPolyData > m_CurrentFiber;

    private:
        /**
         * @brief Image of regions of interest (seed regions for tractography algorithm).
         */
        Self::LabelImage::Pointer m_RegionsOfInterest;

        /**
         * @brief Labels corresponding to seed region in image of region of interest where the tractography algorithm can process.
         */
        std::vector< unsigned short > m_SeedLabels;

        /**
         * @brief Space between seeds (in mm).
         */
        float m_SeedSpacing;

        /**
         * @brief Estimated fibers for each label represented by VTK polydata lines.
         */
        std::vector< vtkSmartPointer< vtkAppendPolyData > > m_OutputFibers;

        /**
         * @brief Indices of labels corresponding to output solutions.
         */
        std::vector< int > m_OutputIndicesOfLabels;
};

} // namespace btk

#endif // BTK_TRACTOGRAPHY_ALGORITHM_H
