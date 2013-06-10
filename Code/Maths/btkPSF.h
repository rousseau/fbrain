/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 27/05/2013
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

#ifndef BTKPSF_H
#define BTKPSF_H

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkFixedArray.h"
#include "itkPoint.h"
#include "itkMatrix.h"
#include "itkVector.h"
#include "itkSize.h"
#include "itkImage.h"

#include "btkMacro.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_math.h"

namespace btk
{

class PSF: public itk::Object
{
    public:
        typedef btk::PSF                        Self;
        typedef itk::Object                     Superclass;

        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;
        typedef double                          OutputType;
        typedef itk::Point< double, 3 >         InputType;

        /** Direction type. */
        typedef itk::Matrix< double,3,3 >  DirectionType;

        /** Point type */
        typedef itk::Point< double,3 >  PointType;

        /** Spacing type */
        typedef itk::Vector< double,3 >  SpacingType;

        typedef itk::Image< float, 3 >  ImageType;

        typedef itk::Size < 3 >         SizeType;

        /** Method for creation through the object factory. */
        //itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(btk::PSF, itk::Object);

        virtual OutputType Evaluate(const InputType & position) const = 0;

        virtual void ConstructImage(SizeType _size) = 0;

        /** Sets the position of the PSF. */
        virtual void SetCenter(PointType center)
        {
          m_Center = center.GetVnlVector();
        }

        /** Sets direction of the PSF. */
        virtual void SetDirection(DirectionType direction)
        {
            m_Direction = direction.GetVnlMatrix();
            m_idir = m_Direction.get_column(0);
            m_jdir = m_Direction.get_column(1);
            m_kdir = m_Direction.get_column(2);

        }
        /** Sets spacing. This method changes the standard deviations of the Gaussian
         * function accordingly. */
        virtual void SetSpacing(SpacingType spacing)
        {
          m_Spacing = spacing.GetVnlVector();
        }

        btkGetMacro(PsfImage,ImageType::Pointer);

    protected:
        PSF();
        virtual ~PSF(){}
        void PrintSelf(std::ostream& os, itk::Indent indent) const;

        vnl_matrix< double > m_Direction;
        vnl_vector< double > m_Center;
        vnl_vector< double > m_Spacing;


        vnl_vector<double> m_idir;
        vnl_vector<double> m_jdir;
        vnl_vector<double> m_kdir;

        ImageType::Pointer m_PsfImage;





};

}
#endif // BTKPSF_H
