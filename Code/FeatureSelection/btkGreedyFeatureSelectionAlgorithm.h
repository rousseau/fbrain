/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 15/01/2013
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

#ifndef BTK_GREEDY_FEATURE_SELECTION_H
#define BTK_GREEDY_FEATURE_SELECTION_H

// ITK includes
#include "itkSmartPointer.h"
#include "itkMacro.h"
#include "itkProcessObject.h"

// Local includes
#include "btkFeatureSelectionAlgorithm.h"


namespace btk
{

/**
 * @brief Greedy algorithm used for feature selection.
 * @author Julien Pontabry
 * @date 15/01/2013
 * @ingroup FeatureSelection
 *
 * This is an implementation of the sequential floating search forward-backward (SFS)
 * from Pudil et al. (1994). The initial algorithm has been adapted for our purpose
 * (search minimal set of features based on regression reconstruction error).
 *
 * His complexity is in O(np), where n and p denotes repsectively the max number of features
 * to select and the initial number of features.
 *
 * Pudil et al. (1994) Floating search methods in feature selection. In: Pattern Recognition
 * Letters 15.11, p. 1119-1125.
 */
class GreedyFeatureSelectionAlgorithm : public FeatureSelectionAlgorithm
{
    public:
        typedef GreedyFeatureSelectionAlgorithm  Self;
        typedef FeatureSelectionAlgorithm        Superclass;
        typedef itk::SmartPointer< Self >        Pointer;
        typedef itk::SmartPointer< const Self >  ConstPointer;

        itkTypeMacro(GreedyFeatureSelectionAlgorithm,FeatureSelectionAlgorithm);
        itkNewMacro(Self);


        /**
         * @brief Run the algorithm.
         * @warning Some parts of the main loop are parallelized using OpenMP.
         */
        void Update();

    protected:
        /**
         * @brief Constructor.
         */
        GreedyFeatureSelectionAlgorithm();

        /**
         * @brief Print a message on output stream.
         * @param os Output stream where the message is printed.
         * @param indent Indentation.
         */
        virtual void PrintSelf(std::ostream &os, itk::Indent indent) const;
};

} // namespace btk

#endif // BTK_GREEDY_FEATURE_SELECTION_H
