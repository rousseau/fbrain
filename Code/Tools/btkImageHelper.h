/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 11/04/2012
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


#ifndef __BTK_IMAGEHELPER_H__
#define __BTK_IMAGEHELPER_H__

// STL includes
#include "string"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


namespace btk
{
    /**
     * Helper class for image management (read, write and create operations)
     * @author Marc Schweitzer, Julien Pontabry
     */
    template < class TImage >
    class ImageHelper
    {
        public:

            /**
             * @brief Image reader type.
             */
            typedef itk::ImageFileReader< TImage > ImageReader;

            /**
             * @brief Image writer type.
             */
            typedef itk::ImageFileWriter< TImage > ImageWriter;


            /**
             * @brief Write an image.
             * @param image Image to write.
             * @param fileName File name of the image to write.
             */
            static void WriteImage(typename TImage::Pointer image, std::string &fileName);

            /**
             * @brief Write a vector of images.
             * @param images vector of images to write.
             * @param fileNames File names of the images to write.
             */
            static void WriteImageArray(std::vector< typename TImage::Pointer > &images, std::vector< std::string > &fileNames);

            /**
             * @brief Read an image.
             * @param fileName File name of the image to read.
             * @return A pointer to the image that have been red.
             */
            static typename TImage::Pointer ReadImage(std::string &fileName);

            /**
             * @brief Read a vector of images.
             * @param fileNames File names of the images to read.
             * @return A reference to a vector containing the images that have been red.
             */
            static std::vector< typename TImage::Pointer > &ReadImageArray(std::vector< std::string> &fileNames);

            /**
             * @brief Create a new image in the same physical space of a current image.
             * @param image Image of which physical space will be used for creation.
             * @return New image in the same physical space.
             */
            static typename TImage::Pointer CreateNewFromSpaceOf(typename TImage::Pointer image);
    };

} // namespace btk

#include "btkImageHelper.txx"

#endif // __BTK_IMAGEHELPER_H__
