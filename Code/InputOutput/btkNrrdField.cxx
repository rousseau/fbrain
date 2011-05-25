/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

19 may 2011
< pontabry at unistra dot fr >

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
*/


#include "btkNrrdField.h"


// STL includes
#include "sstream"


namespace btk
{

class btkNrrdField::btkNrrdFieldPriv
{
    public:
        std::string field;
        std::string   key;
        std::string value;

    public:
        void analyze();
};

void btkNrrdField::btkNrrdFieldPriv::analyze()
{
    if(!field.empty())
    {
        unsigned int i = 0;

        // state A
        while(i < field.length() && field[i] != ':' && field[i] != '=')
            i++;


        if(i < field.length() && field[i] == ':')
        {
            // get the key string before going on state B
            key = field.substr(0,i++);

            // state B
            if(i < field.length() && field[i] == '=')
            {
                // get the position of the first character of the value string
                unsigned int pos = ++i;

                // state D
                while(i < field.length() && field[i] != ':' && field[i] != '=')
                    i++;

                // get the value string
                value = field.substr(pos,i-1);
            }
        }
    }
}

btkNrrdField::btkNrrdField(std::string field) : m(new btkNrrdFieldPriv)
{
    m->field = field;
    m->analyze();

    std::stringstream st;
    st << m->key;
    m->key.clear();
    st >> m->key;
}

btkNrrdField::~btkNrrdField()
{
    delete m;
}

std::string btkNrrdField::GetKey() const
{
    return m->key;
}

std::string btkNrrdField::GetValue() const
{
    return m->value;
}

} // namespace btk
