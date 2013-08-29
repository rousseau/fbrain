/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 12/02/2010
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

#include "btkMathFunctions.h"

namespace btk
{

unsigned int MathFunctions::factorial(unsigned int n)
{
    if(n < 2)
        return 1;
    else
        return factorial(n-1) * n;
}
//-------------------------------------------------------------------------------------------------
double MathFunctions::RadiansToDegrees(double rad)
{
    return (rad * 180/M_PI);
}
//-------------------------------------------------------------------------------------------------
double MathFunctions::DegreesToRadians(double deg)
{
    return(deg * M_PI/180);
}
//-------------------------------------------------------------------------------------------------
double MathFunctions::Random()
{
    return static_cast< double >(rand()) / static_cast< double >(RAND_MAX);
}
//-------------------------------------------------------------------------------------------------
double MathFunctions::Random(double min , double max)
{
    return static_cast< double >(rand())/(static_cast< double >(RAND_MAX)/std::abs(max - min)) - std::abs(min);
}
//-------------------------------------------------------------------------------------------------
double MathFunctions::Round(double value)
{
    return floor(value + 0.5);
}
//-------------------------------------------------------------------------------------------------
double MathFunctions::Sinc(double _x)
{
    if(btk::btkFloatingEqual(_x, 0.0))
    {
        return 1.0;
    }
    else
    {
        return (sin(_x)/_x);
    }
}
//-------------------------------------------------------------------------------------------------
float MathFunctions::Sinc(float _x)
{
    if(btk::btkFloatingEqual((float)_x, (float)0.0))
    {
        return 1.0;
    }
    else
    {
        return (sin(_x)/_x);
    }
}
//-------------------------------------------------------------------------------------------------
int MathFunctions::Sinc(int _x)
{
    if(_x == 0)
    {
        return 1.0;
    }
    else
    {
        return (sin(_x)/_x);
    }
}
} // namespace btk
