/*
 * Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 *
 * 15/03/2012
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


#ifndef __BTK_FILENAME_TOOLS_H__
#define __BTK_FILENAME_TOOLS_H__

// STL includes
#include "string"

namespace btk
{

    /**
     * @brief Get the position of the first point of the extension
     * @param filename Filename
     * @return Positopn of the first point of the extension
     * @ingroup InputOutput
     */
    unsigned int GetExtensionPosition(std::string filename)
    {
        unsigned int      i = 0;
        unsigned int    pos = 0;
        unsigned int length = filename.length();
        char          state = 0;
        bool           pass;

        while(i < length)
        {
            switch(state)
            {
                case 0: // start
                    if(filename[i] == '.')
                        state = 1;
                    else
                        pos++;
                    break;

                case 1: // test extensions
                    pass = false;

                    if(i+6 == length)
                    {
                        std::string tmpString = filename.substr(i,i+5);
                        if(tmpString == "nii.gz")
                            pass = true;
                    }
                    else if(i+4 == length)
                    {
                        std::string tmpString = filename.substr(i,i+3);
                        if(tmpString == "nhdr" || tmpString == "nrrd")
                            pass = true;
                    }
                    else if(i+3 == length)
                    {
                        std::string tmpString = filename.substr(i,i+2);
                        if(tmpString == "nii")
                            pass = true;
                    }

                    if(pass)
                    {
                        i = length-1;
                    }
                    else
                    {
                        state = 0;
                        pos  += 2;
                    }
                    break;

                default: // error
                    state = 0;
                    break;
            } // switch

            i++;
        } // while

        return pos;
    }

    /**
     * @brief Extract the radix filename.
     * Exemple: With filename equal to "xxx.yy.zz" or "xxx.yy", this function returns "xxx".
     * @param filename Filename
     * @return Radix of the filename
     * @ingroup InputOutput
     */
    std::string GetRadixOf(std::string filename)
    {
        std::string radix = "";

        unsigned int pos = GetExtensionPosition(filename);
        radix            = filename.substr(0,pos);

        return radix;
    }

    /**
     * @brief Extract the extension of the filename.
     * Exemple: With filename equal to "xxx.yy.zz" or "xxx.yy", this function returns ".yy.zz" or ".yy".
     * @param filename Filename
     * @return Extension of the filename
     * @ingroup InputOutput
     */
    std::string GetExtensionOf(std::string filename)
    {
        std::string extension = "";

        unsigned int pos = GetExtensionPosition(filename);
        extension        = filename.substr(pos,filename.length()-1);

        return extension;
    }

} // namespace btk

#endif /*__BTK_FILENAME_TOOLS_H__*/
