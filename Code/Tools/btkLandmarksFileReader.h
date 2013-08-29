/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 18/02/2011
  Author(s): Julien Pontabry < pontabry at unistra dot fr >

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


#ifndef __BTK_LANDMARKS_FILE_READER_H__
#define __BTK_LANDMARKS_FILE_READER_H__


// ITK includes
#include "itkMacro.h"
#include "itkObjectFactory.h"
#include "itkLightObject.h"
#include "itkSmartPointer.h"

// STL includes
#include "string"
#include "vector"


namespace btk
{
    /**
     * @class btkLandmarksFileReader
     * @brief The btkLandmarksFileReader class
     * @author Julien Pontabry
     * @ingroup Tools
     */
    class btkLandmarksFileReader : public itk::LightObject
    {
        public:
        typedef btkLandmarksFileReader        Self;
        typedef itk::LightObject              Superclass;
        typedef itk::SmartPointer<Self>       Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

        void SetInputFileName(const std::string &filename);
        void Update();
        double *GetOutputLPT();
        double *GetOutputRPT();
        double *GetOutputPPT();
        double *GetOutputAPT();

        protected:
            btkLandmarksFileReader();

        private:
            inline int SplitString(std::string &string, char separator, std::vector<std::string> &strings);
            inline bool IsNumber(const std::string &string);

        private:
            std::string m_filename;

            double m_lpt[3];
            double m_rpt[3];
            double m_ppt[3];
            double m_apt[3];
    };

} // namespace btk

#endif // __BTK_LANDMARKS_FILE_READER_H__
