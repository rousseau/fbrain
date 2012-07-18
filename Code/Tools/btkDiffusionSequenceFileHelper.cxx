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

#include "btkDiffusionSequenceFileHelper.h"


// STL includes
#include "iostream"
#include "fstream"

// Local includes
#include "btkMacro.h"


namespace btk
{

std::vector< GradientDirection > &DiffusionSequenceFileHelper::ReadGradientTable(const std::string &gradientTableFileName)
{
    std::cout << "Reading file \"" << gradientTableFileName << "\"..." << std::flush;

    std::vector< GradientDirection > *gradientTable_ptr = new std::vector< GradientDirection >;
    std::vector< GradientDirection > &gradientTable = *gradientTable_ptr;
    std::fstream file(gradientTableFileName.c_str(), std::fstream::in);

    if(!file.is_open())
    {
        // FIXME Use exceptions
        btkCerrMacro("Error: unable to read file !");
        std::exit(EXIT_FAILURE);
    }
    else // file is open
    {
        float value = 0;
        std::vector< float > values;

        while((file >> value))
            values.push_back(value);

        unsigned int  nbVectors = values.size()/3;
        unsigned int nbVectors2 = nbVectors + nbVectors;

        // TODO Define GradientDirection
        for(unsigned int i = 0; i < nbVectors; i++)
            gradientTable.push_back(GradientDirection(values[i], values[i+nbVectors], values[i+nbVectors2]));

        file.close();
    }

    std::cout << "done." << std::endl;

    return gradientTable;
}

void DiffusionSequenceFileHelper::WriteGradientTable(std::vector< GradientDirection > &gradientTable, const std::string &gradientTableFileName)
{
    std::cout << "Writing file \"" << gradientTableFileName << "\"..." << std::flush;

    std::fstream file(gradientTableFileName.c_str(), std::fstream::out);

    if(!file.is_open())
    {
        // FIXME Use exceptions
        btkCerrMacro("Error: unable to write file !");
        std::exit(EXIT_FAILURE);
    }
    else // file is open
    {
        // FIXME : set the precision for writing floating numbers
        unsigned int nbGradients = gradientTable.size();

        // X coordinate
        for(unsigned int i = 0; i < nbGradients; i++)
            file << gradientTable[i][0] << " ";

        file << "\n";

        // Y coordinate
        for(unsigned int i = 0; i < nbGradients; i++)
            file << gradientTable[i][1] << " ";

        file << "\n";

        // Z coordinate
        for(unsigned int i = 0; i < nbGradients; i++)
            file << gradientTable[i][2] << " ";

        file << "\n";

        file.close();
    }

    std::cout << "done." << std::endl;
}

std::vector< unsigned short > &DiffusionSequenceFileHelper::ReadBValues(const std::string &bValuesFileName)
{
    std::cout << "Reading file \"" << bValuesFileName << "\"..." << std::flush;

    std::vector< unsigned short > *bValues_ptr = new std::vector< unsigned short >;
    std::vector< unsigned short >     &bValues = *bValues_ptr;
    std::fstream file(bValuesFileName.c_str(), std::fstream::in);

    if(!file.is_open())
    {
        // FIXME Use exceptions
        btkCerrMacro("Error: unable to read file !");
        std::exit(EXIT_FAILURE);
    }
    else // file is open
    {
        unsigned int value = 0;

        while((file >> value))
            bValues.push_back(value);

        file.close();
    }

    std::cout << "done." << std::endl;

    return bValues;
}

void DiffusionSequenceFileHelper::WriteBValues(std::vector< unsigned short > &bValues, const std::string &bValuesFileName)
{
    std::cout << "Writing file \"" << bValuesFileName << "\"..." << std::flush;

    std::fstream file(bValuesFileName.c_str(), std::fstream::out);

    if(!file.is_open())
    {
        // FIXME Use exceptions
        btkCerrMacro("Error: unable to write file !");
        std::exit(EXIT_FAILURE);
    }
    else // file is open
    {
        unsigned int nbValues = bValues.size();

        for(unsigned int i = 0; i < nbValues; i++)
            file << bValues[i] << " ";

        file << "\n";

        file.close();
    }

    std::cout << "done." << std::endl;
}

} // namespace btk
