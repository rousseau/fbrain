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


#ifndef BTK_DIFFUSION_SEQUENCE_HELPER_H
#define BTK_DIFFUSION_SEQUENCE_HELPER_H

// STL includes
#include "string"

// Local includes
#include "btkDiffusionSequence.h"

namespace btk
{
    /**
     * Helper class for diffusion sequence management (read, write and create operations)
     * @author Julien Pontabry
     */
    class DiffusionSequenceHelper
    {
        public:

            /**
             * @brief Write a diffusion sequence.
             * @param sequence Diffusion sequence to write.
             * @param fileName File name of the diffusion sequence to write.
             */
            static void WriteSequence(btk::DiffusionSequence::Pointer sequence, const std::string &fileName);

            /**
             * @brief Write a vector of diffusion sequences.
             * @param sequences vector of diffusion sequences to write.
             * @param fileNames File names of the diffusion sequences to write.
             */
            static void WriteSequenceArray(std::vector< btk::DiffusionSequence::Pointer > &sequences, std::vector< std::string > &fileNames);

            /**
             * @brief Read a diffusion sequence.
             * @param fileName File name of the diffusion sequence to read.
             * @return A pointer to the diffusion sequence that have been red.
             */
            static btk::DiffusionSequence::Pointer ReadSequence(const std::string &fileName);

            /**
             * @brief Read a vector of diffusion sequences.
             * @param fileNames File names of the diffusion sequences to read.
             * @return A reference to a vector containing the diffusion sequences that have been red.
             */
            static std::vector< btk::DiffusionSequence::Pointer > &ReadSequenceArray(std::vector< std::string> &fileNames);

            /**
             * @brief Create a new image in the same physical space of a current diffusion sequence.
             * @param sequence Diffusion sequence of which physical space will be used for creation.
             * @return New image in the same physical space.
             */
            static btk::DiffusionSequence::Pointer CreateNewSequenceFromPhysicalSpaceOf(btk::DiffusionSequence::Pointer sequence);

            /**
             * @brief Create new images in the same physical space of current diffusion sequences.
             * @param sequences Vector of diffusion sequences of which physical space will be used for creation.
             * @return Vector of new images in the same physical space.
             */
            static std::vector< btk::DiffusionSequence::Pointer > &CreateNewSequenceFromPhysicalSpaceOf(std::vector< btk::DiffusionSequence::Pointer > &sequences);
    };

} // namespace btk

#include "btkDiffusionSequenceHelper.txx"

#endif // BTK_DIFFUSION_SEQUENCE_HELPER_H
