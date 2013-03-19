/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/04/2012
  Author(s): Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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

#ifndef __BTK_IOTRANSFORMHELPER_H__
#define __BTK_IOTRANSFORMHELPER_H__

#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

namespace btk
{
/**
 * Helper class for Reading and Writing Transforms
 * This Class is templated over Transform Type, it is used only for reading (the read class while be cast into Templated argument)
 * @author Marc Schweitzer
 *
 */
template <class TTransform /*= itk::Transform<double, 3 >*/ >
class IOTransformHelper
{
public :
    /**
     * @brief Templated transform
     */
    typedef TTransform TransformType;
    /**
     * @brief Pointer of templated transform
     */
    typedef typename TransformType::Pointer TransformPointerType;
    /**
     * @brief Transform reader type
     */
    typedef itk::TransformFileReader itkTReader;
    /**
     * @brief Transform writer type
     */
    typedef itk::TransformFileWriter itkTWriter;
    /**
     * @brief Read a transform
     * @param fileName of the transform to read
     * @return A pointer to the transform (templated)
     */
    static TransformPointerType ReadTransform(std::string fileName);
    /**
     * @brief Read a  vector of transforms
     * @param fileNames of the transforms to read
     * @return A array of pointer of transforms (templated)
     */
    static std::vector< TransformPointerType > & ReadTransform(std::vector< std::string > &fileNames );
    /**
     * @brief Write a transform
     * @param A pointer to the transform
     * @param The name of the transform file
     */
    static void WriteTransform(TransformPointerType transform, const std::string fileName );
    /**
     * @brief Write a array of transforms
     * @param A array of pointer of transform
     * @param A array of name of transform file
     */
    static void WriteTransform(std::vector< TransformPointerType > &transforms, std::vector< std::string > &fileNames);

};

}


#include "btkIOTransformHelper.txx"
#endif
