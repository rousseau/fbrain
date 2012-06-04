#ifndef __BTK_IMAGEHELPER_TXX__
#define __BTK_IMAGEHELPER_TXX__

#include "btkImageHelper.h"

// STL includes
#include "iostream"

namespace btk
{

    template < class TImage >
    void ImageHelper< TImage >::WriteImage(typename TImage::Pointer image, std::string &fileName)
    {
        std::cout << "Writing \"" << fileName << "\"... " << std::flush;

        typename ImageWriter::Pointer writer = ImageWriter::New();
        writer->SetFileName(fileName);
        writer->SetInput(image);
        writer->Update();

        std::cout << "done." << std::endl;
    }

    //----------------------------------------------------------------------------------------

    template < class TImage >
    void ImageHelper< TImage >::WriteImageArray(std::vector< typename TImage::Pointer > &images, std::vector<std::string> &fileNames)
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

    template < class TImage >
    typename TImage::Pointer ImageHelper<TImage>::ReadImage(std::string &fileName)
    {
        std::cout << "Reading image \"" << fileName << "\"... " << std::flush;

        typename ImageReader::Pointer reader = ImageReader::New();
        reader->SetFileName(fileName);
        reader->Update();

        std::cout << "done." << std::endl;

        return reader->GetOutput();
    }

    //----------------------------------------------------------------------------------------

    template < class TImage >
    std::vector< typename TImage::Pointer > &ImageHelper<TImage>::ReadImageArray(std::vector<std::string> &fileNames)
    {
        int i = fileNames.size();
        std::vector< typename TImage::Pointer > *ptrImages = new std::vector< typename TImage::Pointer >;
        std::vector< typename TImage::Pointer > &images = *ptrImages;
        images.resize(i);

        for(i = 0; i < fileNames.size(); i++)
        {
            std::cout << "Reading image \"" << fileNames[i] << "\"... " << std::flush;

            typename ImageWriter::Pointer reader = ImageWriter::New();
            reader->SetFileName(fileNames[i]);
            reader->Update();

            std::cout << "done." << std::endl;

            images[i] = reader->GetOutput();
        }

        return images;
    }

    //----------------------------------------------------------------------------------------

    template < class TImage >
    typename TImage::Pointer ImageHelper<TImage>::CreateNewFromSpaceOf(typename TImage::Pointer image)
    {
        typename TImage::Pointer newImage = TImage::New();
        newImage->SetRegions(image->GetLargestPossibleRegion());
        newImage->SetOrigin(image->GetOrigin());
        newImage->SetSpacing(image->GetSpacing());
        newImage->SetDirection(image->GetDirection());
        newImage->Allocate();
        newImage->FillBuffer(0);

        return newImage;
    }

} // namespace btk

#endif // __BTK_IMAGEHELPER_TXX__
