/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 19/02/2013
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

#ifndef BTK_SPACE_REDUCTION_COMMAND_ITERATION_UPDATE_H
#define BTK_SPACE_REDUCTION_COMMAND_ITERATION_UPDATE_H

// ITK includes
#include "itkCommand.h"

// Local includes
#include "btkFeatureSelectionAlgorithm.h"


namespace btk
{

/**
 * @brief Commandline iteration observer for feature selection algorithm.
 * @author Julien Pontabry
 * @date 19/02/2013
 * @ingroup FeatureSelection
 */
class FeatureSelectionCommandIterationUpdate : public itk::Command
{
    public:
        typedef FeatureSelectionCommandIterationUpdate Self;
        typedef itk::Command                           Superclass;
        typedef itk::SmartPointer< Self >              Pointer;

        itkNewMacro(Self);

    public:
        /**
         * @brief Execute an event from a caller.
         * @param caller Object which raised an event
         * @param event Event which was raised
         */
        void Execute(itk::Object *caller, const itk::EventObject & event)
        {
            Execute((const itk::Object *)caller, event);
        }

        /**
         * @brief Execute an event from a caller.
         * @param caller Object which raised an event
         * @param event Event which was raised
         */
        void Execute(const itk::Object * object, const itk::EventObject & event)
        {
            FeatureSelectionAlgorithm::ConstPointer optimizer = dynamic_cast< const FeatureSelectionAlgorithm * >(object);

            if(!itk::IterationEvent().CheckEvent(&event))
            {
                return;
            }

            std::cout << "-> Iteration "<< optimizer->GetCurrentIteration() << std::endl;
            std::cout << '\t' << optimizer->GetCurrentMessage() << std::endl;
        }

    protected:
        FeatureSelectionCommandIterationUpdate(){}
};

}// Namespace btk

#endif // BTK_SPACE_REDUCTION_COMMAND_ITERATION_UPDATE_H
