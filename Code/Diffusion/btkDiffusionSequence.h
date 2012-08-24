/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 05/07/2012
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

#ifndef BTK_DIFFUSION_SEQUENCE_H
#define BTK_DIFFUSION_SEQUENCE_H

// STL includes
#include "vector"

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkImage.h"

// Local includes
#include "btkMacro.h"
#include "btkGradientDirection.h"

namespace btk
{

/**
 * @brief Represent a diffusion weighted MRI dataset.
 * @author Julien Pontabry
 * @ingroup Diffusion
 */
class DiffusionSequence : public itk::Image< short,4 >
{
    public:
        typedef DiffusionSequence               Self;
        typedef itk::Image< short,4 >           Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef std::vector< GradientDirection > GradientTable;

        itkNewMacro(Self);

        itkTypeMacro(DiffusionSequence, itk::Image);

        btkGetMacro(GradientTable, std::vector< GradientDirection >);
        btkSetMacro(GradientTable, std::vector< GradientDirection >);

        btkGetMacro(BValues, std::vector< unsigned short >);
        btkSetMacro(BValues, std::vector< unsigned short >);

        // TODO : check how to make a good usage of gradient table
//        void UseWorldCoordinatesForGradientTable(); // Change gradient table to world coordinates if necessary
//        void UseImageCoordinatesForGradientTable(); // Change gradient table to image coordinates if necessary
//        void RotateGradientTable(Matrix m); // Rotate the gradient table with a rotation matrix (or transform)

    protected:
        /**
         * @brief Constructor.
         */
        DiffusionSequence();

        /**
         * @brief Destructor.
         */
        virtual ~DiffusionSequence();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:
        /** Gradient table of the diffusion sequence. */
        GradientTable m_GradientTable;

        /** B-values of the diffusion sequence. */
        std::vector< unsigned short > m_BValues;
};

} // namespace btk

#endif // BTK_DIFFUSION_SEQUENCE_H
