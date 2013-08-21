/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 28/05/2013
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

#ifndef BTKHYBRIDPSF_H
#define BTKHYBRIDPSF_H

#include "btkPSF.h"
#include "btkMacro.h"


#include "itkFixedArray.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"

namespace btk
{
/**
 * @brief Hybrid PSF represent several PSF fonction, one fonction on each axes.
 * For example a boxcar on x, and a gaussian on y and z.
 * * A common use of this class is :
 * - first construct the image with the correct parameters (center, size, spacing, fonction for each axis)
 * - Get the Psf Image previously constructed
 * - iterate over the image to get the values
 *
 * NOTE : The Evaluate Method is currently not correct.
 * TODO : An Evaluate Method for compute the value of a point
 * (for example iterate over the image and break when the iterated point == point we wanted, then return the value of this point )
 */
class HybridPSF : public btk::PSF
{
    public:
        typedef btk::HybridPSF                    Self;
        typedef btk::PSF                        Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef Superclass::OutputType          OutputType;
        typedef Superclass::InputType           InputType;
        typedef Superclass::DirectionType       DirectionType;
        typedef Superclass::SpacingType         SpacingType;
        typedef Superclass::PointType           PointType;

        typedef itk::Image < float, 3 >         ImageType;
        typedef itk::ImageRegionIteratorWithIndex< ImageType > itkIteratorWithIndex;
        typedef itk::ContinuousIndex<double,3>     itkContinuousIndex;

        /** Array type */
        typedef itk::FixedArray< double,3 > ArrayType;

        /** typedef of a pointer of function (param float and return float) */
        typedef float(HybridPSF::*function_type)(float);

        /** Type of fonction */
        enum FUNCTION_TYPE
        {
            BOXCAR = 0,
            GAUSSIAN,
            SINC
        };

        /** AXIS (X,Y,Z) */
        enum AXIS
        {
            X = 0,
            Y = 1,
            Z = 2
        };

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(btk::HybridPSF, btk::PSF);

        /** Construct Image */
        virtual void ConstructImage();

        /** Update Sigma with spacing and size (for the gaussian function) */
        virtual void UpdateSigma()
        {
            m_Sigma[0] = m_Spacing[0] * m_Size[0] / 2.3548;
            m_Sigma[1] = m_Spacing[1] * m_Size[1] / 2.3548;
            m_Sigma[2] = m_Spacing[2] * m_Size[2] / 2.3548;
        }
        /** Set Size */
        void SetSize(SizeType _size)
        {
            m_Size =  _size;
            this->UpdateSigma();
        }
        /** Set Spacing */
        virtual void SetSpacing(SpacingType spacing)
        {
            m_Spacing = spacing.GetVnlVector();
            this->UpdateSigma();
        }


        /** Set X axis fonction, _f is a FONCTION_TYPE */
        void SetXFunction(FUNCTION_TYPE _f)
        {
            switch(_f)
            {
                case BOXCAR:
                    m_Functions[X] = &HybridPSF::functionBoxCar;
                    break;

                case GAUSSIAN:
                    m_Functions[X] = &HybridPSF::functionGauss;
                    break;

                case SINC:
                    m_Functions[X] = &HybridPSF::functionSinc;
                    break;

                default:
                    m_Functions[X] = &HybridPSF::functionBoxCar;
                    break;
            };


        }
        /** Set Y axis fonction, _f is a FONCTION_TYPE */
        void SetYFunction(FUNCTION_TYPE _f)
        {
            switch(_f)
            {
                case BOXCAR:
                    m_Functions[Y] = &HybridPSF::functionBoxCar;
                    break;

                case GAUSSIAN:
                    m_Functions[Y] = &HybridPSF::functionGauss;
                    break;

                case SINC:
                    m_Functions[Y] = &HybridPSF::functionSinc;
                    break;

                default:
                    m_Functions[1] = &HybridPSF::functionBoxCar;
                    break;
            };
        }
        /** Set Z axis fonction, _f is a FONCTION_TYPE */
        void SetZFunction(FUNCTION_TYPE _f)
        {
            switch(_f)
            {
                case BOXCAR:
                    m_Functions[Z] = &HybridPSF::functionBoxCar;
                    break;

                case GAUSSIAN:
                    m_Functions[Z] = &HybridPSF::functionGauss;
                    break;

                case SINC:
                    m_Functions[Z] = &HybridPSF::functionSinc;
                    break;

                default:
                    m_Functions[Z] = &HybridPSF::functionBoxCar;
                    break;
            };
        }

    protected:
        /** Constructor */
        HybridPSF();
        /** Destructor */
        virtual ~HybridPSF(){}
        /** Print */
        void PrintSelf(std::ostream& os, itk::Indent indent) const;
        /** Intitialize (empty) */
        virtual void Initialize(){}
        /** Sinc fonction */
        float functionSinc(float _x) ;
        /** Gaussian fonction */
        float functionGauss(float _x) ;
        /** BoxCar fonction */
        float functionBoxCar(float _x);

    private:

        ArrayType m_Sigma; /** Sigma */
        std::vector< function_type > m_Functions; /** Vector of fonctions, one for each axis */
        mutable unsigned int m_Axis; /** Axis */
};
}
#endif // BTKHYBRIDPSF_H
