/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 22/08/2012
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

#ifndef BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_H
#define BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_H

// ITK includes
#include "itkMacro.h"
#include "itkSmartPointer.h"

// Local includes
#include "btkMacro.h"
#include "btkTractographyAlgorithm.h"

namespace btk
{

/**
 * @brief Define a deterministic streamline propagation algorithm.
 * @author Julien Pontabry
 * @ingroup Tractography
 */
class StreamlineTractographyAlgorithm : public btk::TractographyAlgorithm
{
    public:
        typedef StreamlineTractographyAlgorithm Self;
        typedef btk::TractographyAlgorithm      Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef Superclass::PhysicalPoint PhysicalPoint;

        itkNewMacro(Self);

        itkTypeMacro(StreamlineTractographyAlgorithm,btk::TractographyAlgorithm);

        btkSetMacro(StepSize,float);
        btkGetMacro(StepSize,float);

        void UseRungeKuttaOrder4(bool arg)
        {
            m_UseRungeKuttaOrder4 = arg;
        }

    protected:
        /**
         * @brief Constructor.
         */
        StreamlineTractographyAlgorithm();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;

        /**
         * @brief Propagate using the tractography algorithm at a seed point.
         * @param point Seed point.
         */
        virtual void PropagateSeed(Self::PhysicalPoint point);

    private:
        /**
         * @brief Step size between two points of output.
         */
        float m_StepSize;

        /**
         * @brief True if the RK4 method is used or false if the RK1 (Euler) method is used.
         */
        bool m_UseRungeKuttaOrder4;
};

} // namespace btk

#endif // BTK_STREAMLINE_TRACTOGRAPHY_ALGORITHM_H
