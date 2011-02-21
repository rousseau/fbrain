/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 18/02/2011
  Author(s): Julien Pontabry < pontabry at unistra dot fr >

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


#include "btkLandmarksFileReader.h"

// STL includes
#include "fstream"
#include "cctype"
#include "cstdlib"


namespace btk
{

btkLandmarksFileReader::btkLandmarksFileReader()
{
    m_filename = "";

    for(unsigned int i=0; i<4; i++)
    {
        m_lpt[i] = 0;
        m_rpt[i] = 0;
        m_ppt[i] = 0;
        m_apt[i] = 0;
    }
}

void btkLandmarksFileReader::SetInputFileName(const std::string &filename)
{
    m_filename = filename;
}

void btkLandmarksFileReader::Update()
{
    std::fstream file(m_filename.c_str(), std::fstream::in);

    if(!file.is_open())
    {
        std::cout << "Error: Unable to read the landmarks file !" << std::endl;
    }
    else // file is open
    {
        std::string s;
        unsigned int found = 0;

        while(found < 4 && std::getline(file,s))
        {
            if(s.size() >= 1 && s[0] != '#')
            {
                std::vector<std::string> splitted;

                if(this->SplitString(s, ',', splitted) >= 3)
                {
                    unsigned int triplet = 0;
                    unsigned int i       = 0;

                    while(triplet < 3 && i<splitted.size())
                    {
                        if(this->IsNumber(splitted[i]))
                        {
                            if(found == 0)
                                m_lpt[triplet] = atof(splitted[i].c_str());
                            else if(found == 1)
                                m_rpt[triplet] = atof(splitted[i].c_str());
                            else if(found == 2)
                                m_ppt[triplet] = atof(splitted[i].c_str());
                            else
                                m_apt[triplet] = atof(splitted[i].c_str());

                            triplet++;
                        }

                        i++;
                    }

                    if(triplet != 3)
                    {
                        std::cout << "Error: There is a mistake in the landmarks file !" << std::endl;
                    }
                }
            }
        }

        if(found != 4)
        {
            std::cout << "Error: There is a mistake in the landmarks file !" << std::endl;
        }

        file.close();
    }
}

double *btkLandmarksFileReader::GetOutputLPT()
{
    return m_lpt;
}

double *btkLandmarksFileReader::GetOutputRPT()
{
    return m_rpt;
}

double *btkLandmarksFileReader::GetOutputPPT()
{
    return m_ppt;
}

double *btkLandmarksFileReader::GetOutputAPT()
{
    return m_apt;
}

int btkLandmarksFileReader::SplitString(std::string &string, char separator, std::vector<std::string> &strings)
{
    strings.clear();

    std::string::size_type spPos = string.find(separator);

    while(spPos != std::string::npos)
    {
        strings.push_back(string.substr(0,spPos));
        string = string.substr(spPos+1);
        spPos  = string.find(separator);
    }

    strings.push_back(string);

    return strings.size();
}

bool btkLandmarksFileReader::IsNumber(const std::string &string)
{
    bool isNumber  = true;
    unsigned int i = 0;

    while(isNumber && i<string.size())
    {
        if(!std::isdigit(string[i]) && string[i] != '.')
            isNumber = false;

        i++;
    }

    return isNumber;
}

} // namespace btk
