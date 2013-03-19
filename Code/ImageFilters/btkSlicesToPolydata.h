/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 08/03/2013
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

#ifndef BTKSLICESTOPOLYDATA_H
#define BTKSLICESTOPOLYDATA_H

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

template <typename TImage>
class SlicesToPolyData : public itk::Object
{

public:
        typedef btk::SlicesToPolyData<TImage> Self;
        typedef itk::Object Superclass;
        typedef itk::SmartPointer< Self > Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;
        typedef TImage ImageType;

        typedef itk::Transform<double, 3,3> Transform;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);
        /** Update Method */
        virtual void Update();

        /** Set/Get input image*/
        btkSetMacro(Input,typename ImageType::Pointer);
        btkGetMacro(Input, typename ImageType::Pointer);

        /** Set/Get transform */
        btkSetMacro(Transform, Transform*);
        btkGetMacro(Transform, Transform*);

        /** Get Output PolyData */
        btkGetMacro(Output,vtkSmartPointer< vtkPolyData > );


    protected :
        /** Initialize Method to call before Update() */
        virtual void Initialize();
        SlicesToPolyData();
        virtual ~SlicesToPolyData();

    private :
        typename ImageType::Pointer m_Input;
        Transform* m_Transform;

        vtkSmartPointer< vtkPolyData > m_Output;


};

} //end namespace

#include "btkSlicesToPolydata.txx"

#endif // BTKSLICESTOPOLYDATA_H
