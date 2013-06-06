/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 
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

#ifndef BTKGAUSSIANPSF_H
#define BTKGAUSSIANPSF_H

#include "btkPSF.h"

#include "itkGaussianSpatialFunction.h"
#include "itkFixedArray.h"

namespace btk
{


class GaussianPSF : public btk::PSF
{
    public:
        typedef btk::GaussianPSF                Self;
        typedef btk::PSF                        Superclass;

        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef Superclass::OutputType          OutputType;
        typedef Superclass::InputType           InputType;
        typedef Superclass::DirectionType       DirectionType;
        typedef Superclass::SpacingType         SpacingType;
        typedef Superclass::PointType           PointType;

        /** Gaussian function type */
        typedef itk::GaussianSpatialFunction< double,
                                         3,
                                         PointType > GaussianFunctionType;

        /** Array type */
        typedef itk::FixedArray< double,3 > ArrayType;


        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(btk::GaussianPSF, btk::PSF);

        virtual OutputType Evaluate(const InputType & position) const;


        virtual void SetSpacing(SpacingType spacing)
        {
          m_Spacing = spacing.GetVnlVector();

          m_Sigma[0] = m_Spacing[0] / 2.3548;
          m_Sigma[1] = m_Spacing[1] / 2.3548;
          m_Sigma[2] = m_Spacing[2] / 2.3548;

//          sigma[0] = sqrt(m_Spacing[0]*m_Spacing[0]/(8*log(2)));
//          sigma[1] = sqrt(m_Spacing[1]*m_Spacing[1]/(8*log(2)));
//          sigma[2] = sqrt(m_Spacing[2]*m_Spacing[2]/(8*log(2)));

      //    sigma[0] = 0.5*m_Spacing[0];
      //    sigma[1] = 0.5*m_Spacing[1];
      //    sigma[2] = 0.5*m_Spacing[2];



          //m_Gaussian -> SetSigma( m_Sigma );

        }


    protected:
        GaussianPSF();
        virtual ~GaussianPSF(){}
        void PrintSelf(std::ostream& os, itk::Indent indent) const;


    private:

        GaussianFunctionType::Pointer m_Gaussian;

        ArrayType m_Sigma;

};

}//end namespace
#endif // BTKGAUSSIANPSF_H
