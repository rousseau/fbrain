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
        IdVector    m_Vec;
        unsigned int m_Size;
        vtkPolyData *m_Polydata;
        vtkIdType   m_SeedPoint;
        unsigned int m_Nb_Loop;


        void growing(vtkIdType m_SeedPoint);
        IdVector getNeighbours(vtkIdType point);
        std::vector<double> getNormalVector(vtkIdType point);
        double getCurvValue(vtkIdType pointID);
        bool isInside(vtkIdType value, IdVector vec);
};

#endif // BTKREGIONGROW_H
