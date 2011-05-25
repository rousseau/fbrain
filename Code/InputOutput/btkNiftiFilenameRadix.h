/*
 * Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 *
 * 23 may 2011
 * < pontabry at unistra dot fr >
 *
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-B license and that you accept its terms.
 */


#ifndef __BTK_NIFTI_FILENAME_RADIX_H__
#define __BTK_NIFTI_FILENAME_RADIX_H__

// STL includes
#include "string"


namespace btk
{

    /**
     * @brief Extract the radix in nifti filename.
     * Exemple: With filename equal to "data.nii.gz" or "data.nii", this function returns "data".
     * @param filename Nitfi filename
     * @return Radix of the Nitfi filename
     */
    std::string GetRadixOf(std::string filename)
    {
        std::string radix = "";

        unsigned int      i = 0;
        unsigned int    pos = 0;
        unsigned int length = filename.length();
        char          state = 0;

        while(i < length && state != 5)
        {
            switch(state)
            {
                case 0: // start
                    if(filename[i] == '.')
                        state = 1;
                    else
                        pos++;
                    break;

                case 1:
                    if(filename[i] == 'n')
                        state = 2;
                    else
                    {
                        state = 0;
                        pos += 2;
                    }
                    break;

                case 2:
                    if(filename[i] == 'i')
                        state = 3;
                    else
                    {
                        state = 0;
                        pos += 3;
                    }
                    break;

                case 3:
                    if(filename[i] == 'i')
                        state = 4;
                    else
                    {
                        state = 0;
                        pos += 4;
                    }
                    break;

                case 4: // final
                    if(filename[i] == '.')
                        state = 6;
                    else
                        state = 5;
                    break;

//                case 5: // error
//                    break;

                case 6:
                    if(filename[i] == 'g')
                        state = 7;
                    else
                        state = 5;
                    break;

                case 7:
                    if(filename[i] == 'z')
                        state = 8;
                    else
                        state = 5;
                    break;

//                case 8: // final
//                    break;

                default:
                    state = 5;
                    break;
            } // switch(state)

            i++;
        } // while

        if(state == 4 || state == 8)
            radix = filename.substr(0,pos);

        return radix;
    }

} // namespace btk

#endif /*__BTK_NIFTI_FILENAME_RADIX_H__*/
