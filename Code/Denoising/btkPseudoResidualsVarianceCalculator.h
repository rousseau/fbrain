/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 25/01/2011
  Author(s): François Rousseau (rousseau@unistra.fr)
             Julien Pontabry   (pontabry@unistra.fr)
             Marc Schweitzer   (marc.schweitzer@unistra.fr)

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


#ifndef BTK_PSEUDO_RESIDUALS_VARIANCE_CALCULATOR_H
#define BTK_PSEUDO_RESIDUALS_VARIANCE_CALCULATOR_H

// ITK includes
#include "itkMacro.h"
#include "itkObject.h"
#include "itkImage.h"

// Local includes
#include "btkMacro.h"
#include "btkDiffusionSignal.h"


namespace btk
{

// TODO : template this class to accept multiple 3D image type (vector image and image).
/**
 * @class PseudoResidualsVarianceCalculator
 * @author François Rousseau
 * @author Julien Pontabry
 * @author Marc Schweitzer
 * This program implements a method to compute the pseudo residuals variances based on the work of Coupé et al. described in : Coupé, P., Yger, P., Prima, S., Hellier, P., Kervrann, C., Barillot, C., 2008. An optimized blockwise nonlocal means denoising filter for 3-D magnetic resonance images. IEEE Transactions on Medical Imaging 27 (4), 425–441.
 */
class PseudoResidualsVarianceCalculator : public itk::Object
{
    public:
        typedef PseudoResidualsVarianceCalculator Self;
        typedef itk::Object                       Superclass;
        typedef itk::SmartPointer< Self >         Pointer;
        typedef itk::SmartPointer< const Self >   ConstPointer;

        typedef itk::Image< unsigned char,3 > MaskImage;
        typedef itk::VectorImage< double,3 >  InternalImage;

        itkNewMacro(Self);
        itkTypeMacro(PseudoResidualsVarianceCalculator,itk::Object);

        btkSetMacro(MaskImage, MaskImage::Pointer);
        btkGetMacro(MaskImage, MaskImage::Pointer);

        btkGetMacro(PseudoResidualsVariance, std::vector< double >);

        btkGetMacro(PseudoResidualsStdDeviation, std::vector< double >);


        /**
         * @brief Set input image.
         * @param signal Input signal image.
         */
        void SetInputImage(DiffusionSignal::Pointer signal);

        /**
         * @brief Compute the pseudo residuals variance.
         */
        void Compute();


    protected:
        /**
         * @brief Constructor.
         */
        PseudoResidualsVarianceCalculator();

        /**
         * @brief Destructor.
         **/
        virtual ~PseudoResidualsVarianceCalculator();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;


    private:
        PseudoResidualsVarianceCalculator(const Self &); // Purposely not implemented
        void operator=(const Self &);                    // Purposely not implemented

        /**
         * @brief Generate data (compute the variance of the pseudo residuals).
         */
        void GenerateData();


    private:
        /**
         * @brief m_InputImage
         */
        InternalImage::Pointer m_InputImage;

        /**
         * @brief m_MaskImage
         */
        MaskImage::Pointer m_MaskImage;

        /**
         * @brief Pseudo residuals variance of the input image.
         */
        std::vector< double > m_PseudoResidualsVariance;

        /**
         * @brief Pseudo residuals standard deviation of the input image.
         */
        std::vector< double > m_PseudoResidualsStdDeviation;

        /**
         * @brief Precomputed constant.
         */
        static const double m_1over6 = 0.166666666666667;

        /**
         * @brief Precomputed constant.
         */
        static const double m_sqrt6over7 = 0.925820099772551;
};

}

#endif // BTK_PSEUDO_RESIDUALS_VARIANCE_CALCULATOR_H
