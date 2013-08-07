/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/07/2013
  Author(s): Aïcha Bentaieb (abentaieb@unistra.fr)

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

#ifndef BTKREGIONGROW_H
#define BTKREGIONGROW_H

#include <vector>
#include "vtkIdList.h"
#include "vtkPolyData.h"
#include "vtkPolyDataAlgorithm.h"
#define PI 3.14159265

class btkRegionGrow : public vtkPolyDataAlgorithm
{

    vtkTypeRevisionMacro(btkRegionGrow,vtkPolyDataAlgorithm);

    typedef std::vector<vtkIdType >   IdVector;

    public:
        static btkRegionGrow *New();
         void setPolydata(vtkPolyData * polydata)
        {
            m_Polydata = polydata;
        }

        void setSeedPoint(vtkIdType seedPoint)
        {
            m_SeedPoint = seedPoint;
        }

        void setNumberOfPoints (int n)
        {
            m_numberOfPoints = n;
        }

        IdVector GetOutput()
        {

            return m_Vec;
        }

        void setNumberIterations(unsigned int Nb)
        {
            m_Nb_Loop = Nb;
        }

        void Update();

    protected:
        btkRegionGrow();



    private:
//*************************************
//AICHA : A QUOI CORRESPONDENT CES VARIABLES ?
//*************************************    
    
        IdVector    m_Vec;
        unsigned int m_Size;
        vtkPolyData *m_Polydata;
        vtkIdType   m_SeedPoint;
        int m_numberOfPoints;
        unsigned int m_Nb_Loop;


        void growing(vtkIdType m_SeedPoint, int m_numberOfPoints);
        IdVector getNeighbours(vtkIdType point);
        std::vector<double> getNormalVector(vtkIdType point);
        double getCurvValue(vtkIdType pointID);
        bool isInside(vtkIdType value, IdVector vec);
};

#endif // BTKREGIONGROW_H
