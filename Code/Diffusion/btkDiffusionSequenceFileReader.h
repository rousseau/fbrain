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

#ifndef BTK_DIFFUSION_SEQUENCE_FILE_READER_H
#define BTK_DIFFUSION_SEQUENCE_FILE_READER_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkImageFileReader.h"

// Local includes
#include "btkDiffusionSequence.h"
#include "btkGradientDirection.h"

namespace btk
{

/**
 * @brief Read a diffusion weighted MRI dataset
 * @author Julien Pontabry
 * @ingroup Diffusion
 */
class DiffusionSequenceFileReader : public itk::ImageFileReader< DiffusionSequence >
{
    public:
        typedef DiffusionSequenceFileReader               Self;
        typedef itk::ImageFileReader< DiffusionSequence > Superclass;
        typedef itk::SmartPointer< Self >                 Pointer;
        typedef itk::SmartPointer< const Self >           ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(DiffusionSequenceFileReader, itk::ImageFileReader);

        /**
         * @brief Update the process (read the diffusion weighted intensities, the gradient table and the b-values).
         * The radix of the name of the trhee files are supposed to be the same.
         */
        virtual void Update();

        /**
         * @brief Get the output diffusion sequence.
         * @return A pointer to the output diffusion sequence.
         */
        virtual DiffusionSequence *GetOutput();

    protected:
        /**
         * @brief Constructor.
         */
        DiffusionSequenceFileReader();

        /**
         * @brief Destructor.
         */
        virtual ~DiffusionSequenceFileReader();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:
        /** Gradient table of the diffusion sequence to be red. */
        std::vector< GradientDirection > m_GradientTable;

        /** B-values of the diffusion sequence to be red. */
        std::vector< unsigned short > m_BValues;
};

} // namespace btk

#endif // BTK_DIFFUSION_SEQUENCE_FILE_READER_H
