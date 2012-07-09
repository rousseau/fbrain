#ifndef __BTK_IMAGEHELPER_TXX__
#define __BTK_IMAGEHELPER_TXX__

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
void ImageHelper< TImageInput, TImageOutput >::WriteImageArray(std::vector< typename TImageInput::Pointer > &images, std::vector<std::string> &fileNames)
{
    int i = images.size();

    if(images.size() == fileNames.size())
    {
        for(i = 0; i < images.size(); i++)
        {
            std::cout << "Writing \"" << fileNames[i] << "\"... " << std::flush;

            typename ImageWriter::Pointer writer = ImageWriter::New();
            writer->SetFileName(fileNames[i]);
            writer->SetInput(images[i]);
            writer->Update();

            std::cout << "done." << std::endl;
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
std::vector< typename TImageInput::Pointer > &ImageHelper< TImageInput, TImageOutput >::ReadImageArray(std::vector<std::string> &fileNames)
{
    int i = fileNames.size();
    std::vector< typename TImageInput::Pointer > *ptrImages = new std::vector< typename TImageInput::Pointer >;
    std::vector< typename TImageInput::Pointer > &images = *ptrImages;
    images.resize(i);

    for(i = 0; i < fileNames.size(); i++)
    {
        std::cout << "Reading image \"" << fileNames[i] << "\"... " << std::flush;

        typename ImageReader::Pointer reader = ImageReader::New();
        reader->SetFileName(fileNames[i]);
        reader->Update();

        std::cout << "done." << std::endl;

        images[i] = reader->GetOutput();
    }

    return images;
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
typename TImageOutput::Pointer ImageHelper< TImageInput, TImageOutput >::CreateNewFromPhysicalSpaceOf(typename TImageInput::Pointer image)
{
    typename TImageOutput::Pointer newImage = TImageOutput::New();
    newImage->SetRegions(image->GetLargestPossibleRegion());
    newImage->SetOrigin(image->GetOrigin());
    newImage->SetSpacing(image->GetSpacing());
    newImage->SetDirection(image->GetDirection());
    newImage->Allocate();
    newImage->FillBuffer(0);

    return newImage;
}

//----------------------------------------------------------------------------------------

template < class TImageInput, class TImageOutput >
std::vector< typename TImageOutput::Pointer > &ImageHelper< TImageInput, TImageOutput >::CreateNewFromPhysicalSpaceOf(std::vector< typename TImageInput::Pointer > &images)
{
    std::vector< typename TImageOutput::Pointer > *ptrNewImages = new std::vector< typename TImageOutput::Pointer >;
    std::vector< typename TImageOutput::Pointer > &newImages = *ptrNewImages;

    for(typename std::vector< typename TImageInput::Pointer >::iterator it = images.begin(); it != images.end(); it++)
    {
        newImages.push_back(CreateNewFromPhysicalSpaceOf(*it));
    }

    return newImages;
}

} // namespace btk

#endif // __BTK_IMAGEHELPER_TXX__
