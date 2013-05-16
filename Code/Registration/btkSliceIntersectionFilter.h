/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 26/04/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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

#ifndef BTKSLICEINTERSECTIONFILTER_H
#define BTKSLICEINTERSECTIONFILTER_H

/* ITK */
#include "itkObject.h"
#include "itkTransform.h"
#include "itkSmartPointer.h"

/* VTK */
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"



/* BTK */
#include "btkMacro.h"

namespace btk
{
/**
 * This class compute the intersection line between two slices
 * The output is the starting Point and the ending point of the line
 * In order to display movement user must provide an itk::Transform, even if transform is identity.
 *
 * @author Marc Schweitzer
 * \ingroup ImageFilters
 */
template <typename TImage>
class SliceIntersectionFilter : public itk::Object
{
    public:
        /** Typedefs */
        typedef btk::SliceIntersectionFilter<TImage> Self;
        typedef itk::Object Superclass;
        typedef itk::SmartPointer< Self > Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;
        typedef TImage ImageType;

        typedef itk::Transform<double, 3,3> Transform;


        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Update Method */
        inline virtual void Update();

        /** Get/Set input Images */
        btkSetMacro(Image1,typename ImageType::Pointer);
        btkGetMacro(Image1,typename ImageType::Pointer);
        btkSetMacro(Image2,typename ImageType::Pointer);
        btkGetMacro(Image2,typename ImageType::Pointer);

        /** Get/Set inputs slices  */
        btkSetMacro(Slice1,unsigned int);
        btkGetMacro(Slice1,unsigned int);
        btkSetMacro(Slice2,unsigned int);
        btkGetMacro(Slice2,unsigned int);

        /** Set/Get transform */
        btkSetMacro(Transform1, Transform*);
        btkGetMacro(Transform1, Transform*);
        btkSetMacro(Transform2, Transform*);
        btkGetMacro(Transform2, Transform*);

        /** Get if there is an intersection */
        btkGetMacro(IsIntersection,bool);

        /** Get Points */
        btkGetMacro(Point1,typename ImageType::PointType);
        btkGetMacro(Point2,typename ImageType::PointType);



    protected:
        SliceIntersectionFilter();
        virtual ~SliceIntersectionFilter(){}
        inline virtual void Initialize();
        inline virtual void CreatePlane1();
        inline virtual void CreatePlane2();

    private:
         unsigned int m_Slice1;
         unsigned int m_Slice2;

         Transform* m_Transform1;
         Transform* m_Transform2;

         typename ImageType::Pointer m_Image1;
         typename ImageType::Pointer m_Image2;

         typename ImageType::PointType m_Point1;
         typename ImageType::PointType m_Point2;

         bool m_IsIntersection;



         vtkSmartPointer< vtkPolyData > m_Plane1;
         vtkSmartPointer< vtkPolyData > m_Plane2;

         vtkSmartPointer< vtkPolyData > m_Line;

         std::vector< std::vector< double > > m_LinePoints;




};
} //btk

#include "btkSliceIntersectionFilter.txx"
#endif // BTKSLICEINTERSECTIONFILTER_H
