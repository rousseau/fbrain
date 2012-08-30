/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 11/04/2012
  Author(s): Marc Schweitzer (marc.schweitzer(at)unistra.fr)
             Julien Pontabry (pontabry@unistra.fr)

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

#ifndef BTK_IMAGE_HELPER_TXX
#define BTK_IMAGE_HELPER_TXX

#include "btkImageHelper.h"


// STL includes
#include "iostream"


namespace btk
{

template < class TImageInput, class TImageOutput >
void ImageHelper< TImageInput, TImageOutput >::WriteImage(typename TImageInput::Pointer image, const std::string &fileName)
{
    std::cout << "Writing \"" << fileName << "\"... " << std::flush;

    typename ImageWriter::Pointer writer = ImageWriter::New();
    writer->SetFileName(fileName);
    writer->SetInput(image);
    writer->Update();

    std::cout << "done." << std::endl;
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
void ImageHelper< TImageInput, TImageOutput >::WriteImage(std::vector< typename TImageInput::Pointer > &images, std::vector<std::string> &fileNames)
{

    if(images.size() == fileNames.size())
    {
        for(int i = 0; i < images.size(); i++)
        {
            WriteImage(images[i], fileNames[i]);
        }
    }
    else
    {
        std::string err("vector of images and vector of names have not the same size !");
        throw err;

    }
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
typename TImageInput::Pointer ImageHelper< TImageInput, TImageOutput >::ReadImage(const std::string &fileName)
{
    std::cout << "Reading image \"" << fileName << "\"... " << std::flush;

    typename ImageReader::Pointer reader = ImageReader::New();
    reader->SetFileName(fileName);
    reader->Update();

    std::cout << "done." << std::endl;

    return reader->GetOutput();
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
std::vector< typename TImageInput::Pointer > &ImageHelper< TImageInput, TImageOutput >::ReadImage(std::vector<std::string> &fileNames)
{
    std::vector< typename TImageInput::Pointer > *ptrImages = new std::vector< typename TImageInput::Pointer >;
    std::vector< typename TImageInput::Pointer > &images = *ptrImages;
    images.resize(fileNames.size());

    for(int i = 0; i < fileNames.size(); i++)
    {
        images[i] = ReadImage(fileNames[i]);
    }

    return images;
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
typename TImageOutput::Pointer ImageHelper< TImageInput, TImageOutput >::CreateNewImageFromPhysicalSpaceOf(typename TImageInput::Pointer image)
{
    typename TImageOutput::Pointer newImage = TImageOutput::New();
    newImage->SetRegions(image->GetLargestPossibleRegion());
    newImage->SetOrigin(image->GetOrigin());
    newImage->SetSpacing(image->GetSpacing());
    newImage->SetDirection(image->GetDirection());
    newImage->Allocate();
//    newImage->FillBuffer(0);

    return newImage;
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
std::vector< typename TImageOutput::Pointer > &ImageHelper< TImageInput, TImageOutput >::CreateNewImageFromPhysicalSpaceOf(std::vector< typename TImageInput::Pointer > &images)
{
    std::vector< typename TImageOutput::Pointer > *ptrNewImages = new std::vector< typename TImageOutput::Pointer >;
    std::vector< typename TImageOutput::Pointer > &newImages = *ptrNewImages;

    for(typename std::vector< typename TImageInput::Pointer >::iterator it = images.begin(); it != images.end(); it++)
    {
        newImages.push_back(CreateNewImageFromPhysicalSpaceOf(*it));
    }

    return newImages;
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
typename TImageOutput::Pointer ImageHelper< TImageInput, TImageOutput >::ReadOrCreateImage(const std::string &fileName, typename TImageInput::Pointer image, typename TImageOutput::PixelType value)
{
  typename TImageOutput::Pointer newImage = TImageOutput::New();

  if(fileName != "")
  {
    newImage = ReadImage(fileName);
  }
  else
  {
    std::cout << "Creating new image with pixel value set to " << value << std::endl;
    newImage = CreateNewImageFromPhysicalSpaceOf(image); 
    newImage->FillBuffer(value);   
  }   
  return newImage; 
}



} // namespace btk

#endif // BTK_IMAGE_HELPER_TXX
