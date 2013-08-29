/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 19/08/2013
  Author(s):Frederic Champ (champ(at)unistra.fr)
  
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

#ifndef BTKDIFFUSIONSLICE_H
#define BTKDIFFUSIONSLICE_H

// ITK includes
#include "itkAffineTransform.h"
#include "itkImage.h"

// Local includes
#include "btkGradientDirection.h"


namespace btk
{

/**
 * @brief Represent a diffusion weighted slice of MRI dataset.
 *
 * it's an tk::Image< short, 3> with:
 *      - A gradient direction
 *      - A B-value
 *      - An outliers status
 *
 * @author Frederic Champ
 * @ingroup Diffusion
 */

class DiffusionSlice : public itk::Image< short, 3 >
{
    public:
        typedef DiffusionSlice                  Self;
        typedef itk::Image< short, 3>           Superclass;
        typedef itk::SmartPointer< Self>        Pointer;
        typedef itk::SmartPointer< const Self>  ConstPointer;

        /** Transform typedefs; */
        typedef itk::AffineTransform<double,3>  TransformType;
        typedef TransformType::Pointer          TransformPointer;

         /** Vector typedefs; */
        typedef itk::Vector<double, 3>          VectorType;

        /** Diffusion parameters typedef.; */
        typedef unsigned int                    BValueType;
        typedef GradientDirection               GradType;

        itkNewMacro(Self);
        itkTypeMacro(DiffusionSlice, itk::Image);


        /**
         * @brief Set Slice image.
         * @param image Input slice image.
         */
        void SetImage(Superclass::ConstPointer image);

        /**
         * @brief Set the gradient direction of the slice.
         * @param gradientDirection The Gradient direction of the slice.
         */
        void SetGradientDirection(GradientDirection gradientDirection)
        {
            m_GradientDirection = gradientDirection;
        }

        /**
         * @brief Get the gradient direction of the slice.
         * @return gradientDirection The Gradient direction of the slice.
         */
        GradientDirection GetGradientDirection()
        {
            return m_GradientDirection;
        }

        /**
         * @brief Set/Get the B-value of the slice.
         */
        btkSetMacro(BValue, BValueType);
        btkGetMacro(BValue, BValueType);

        /**
         * @brief Set Outliers status of the slice.
         * @param OutlierStatus Boolean of outliers status.
         */
        btkSetMacro(OutlierStatus, bool);

        /**
         * @brief Get Outliers status of the slice.
         * @return OutlierStatus Boolean of outliers status.
         */
        bool IsOutlier()
        {
            return m_OutlierStatus;
        }

        /**
         * @brief Transform the slice.
         * @return transform Affine transformation.
         */
        void Transform(TransformPointer transform);

    protected:

        /**
         * @brief Constructor.
         */
        DiffusionSlice();

        /**
         * @brief Destructor.
         */
        virtual ~DiffusionSlice();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

    private:
        /** Slice gradient direction */
        GradientDirection           m_GradientDirection;

        /** Slice B-value */
        BValueType                  m_BValue;

        /** Outliers status */
        bool                        m_OutlierStatus;

};

} //end namespace btk
#endif // BTKDIFFUSIONSLICE_H
