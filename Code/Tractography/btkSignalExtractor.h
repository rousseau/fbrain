/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

12 april 2010
< pontabry at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#ifndef BTK_SIGNALEXTRACTOR_H
#define BTK_SIGNALEXTRACTOR_H

    // STL includes
    #include "string"
    #include "vector"

    // Local includes
    #include "btkTypes.h"
    #include "btkDirection.h"


    namespace btk
    {

    /**
     * @class SignalExtractor
     * @brief Extract signal from a dwi dataset
     * @author Julien Pontabry
     */
    class SignalExtractor
    {
        public:
            /**
              * @brief Constructor
              * Extract signal from the data file using directions files. The signal is
              * normalized (division by the average reference image).
              * @param vectorsFileName Name of the file containing vectors
              * @param dataFileName Signal filename
              * @param maskFileName Image mask filename
              */
            SignalExtractor(const std::string &vectorsFileName, const std::string &dataFileName, const std::string &maskFileName, char displayMode);

            /**
              * @brief Destructor
              */
            ~SignalExtractor();

            /**
              * @brief Extract signal
              * Extract normalized signal from dwi sequence
              */
            void extract();

            /**
              * @brief Save dataset in some files
              * The signal is saved in an nifti image and the directions in spherical coordinates in
              * a text file.
              */
            void save();

            /**
             * @brief Get the gradient directions
             * @return A vector of gradient directions
             */
            std::vector<Direction> *GetDirections();

            /**
             * @brief Get the signal image
             * @return A signal sequence
             */
            Sequence::Pointer GetSignal();

            /**
             * @brief Get the signal's standard deviations
             * @return Signal's standard deviations
             */
            std::vector<Real> *GetSigmas();

            /**
             * @brief Get the mask image
             * @return Mask image
             */
            Mask::Pointer GetMask();

        private:
            /**
              * @brief Read some files
              * Files will be red
              * @param vectorsFileName Dwi vectors filename
              * @param dataFileName Dwi sequence filename
              * @param maskFileName Image mask filename
              */
            void readFiles(const std::string &vectorsFileName, const std::string &dataFileName, const std::string &maskFileName);

            /**
             * @brief Compute standard deviation
             */
            void computeSigmas();

        private:
            std::vector<Direction>    *m_directions;    /**< Set of directions */
            std::vector<unsigned int> *m_refIm;         /**< Index of references images */
            Sequence::Pointer          m_signal;        /**< Signal image */
            Sequence::Pointer          m_data;          /**< DWI sequence */
            Mask::Pointer              m_mask;          /**< Image mask */
            std::vector<Real>         *m_sigmas;        /**< Images' sigmas */

            char m_displayMode;
    };

    } // namespace btk

#endif // BTK_SIGNALEXTRACTOR_H

