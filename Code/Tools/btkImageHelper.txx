#ifndef __BTK_IMAGEHELPER_TXX__
#define __BTK_IMAGEHELPER_TXX__

#include "btkImageHelper.h"
#include "iostream"

namespace btk
{
//----------------------------------------------------------------------------------------
template <class TImage >
void ImageHelper< TImage >::WriteImage(itkImagePointer image, std::string &fileName)
{
    std::cout<<"Writing : "<<fileName;

    itkWriterPointer writer = itkWriter::New();
    writer->SetFileName(fileName);
    writer->SetInput(image);
    writer->Update();
    std::cout<<" ...done !"<<std::endl;

}
//----------------------------------------------------------------------------------------
template <class TImage >
void ImageHelper< TImage >::WriteImageArray(std::vector<itkImagePointer> &images, std::vector<std::string> &fileNames)
{
    int i = images.size();

    if(images.size() == fileNames.size())
    {
        for(i = 0; i < images.size(); i++)
        {
            std::cout<<"Writing : "<<fileNames[i];
            itkWriterPointer writer = itkWriter::New();
            writer->SetFileName(fileNames[i]);
            writer->SetInput(images[i]);
            writer->Update();
            std::cout<<" ...done !"<<std::endl;
        }
    }
    else
    {
        std::string err("vector of images and vector of names have not the same size !");
        throw err;

    }
}
//----------------------------------------------------------------------------------------
template <class TImage >
typename ImageHelper<TImage>::itkImagePointer ImageHelper<TImage>::ReadImage(std::string &fileName)
{
    std::cout<<"Reading image : "<<fileName<<std::endl;
    itkReaderPointer reader = itkReader::New();
    reader->SetFileName(fileName);
    reader->Update();
    itkImagePointer img = reader->GetOutput();
    return img;

}
//----------------------------------------------------------------------------------------
template <class TImage >
void ImageHelper< TImage >::ReadImageArray(std::vector< itkImagePointer > &images, std::vector<std::string> &fileNames)
{
    int i = fileNames.size();
    images.resize(i);

    for(i = 0; i < fileNames.size(); i++)
    {
        std::cout<<"Reading image : "<<fileNames[i]<<std::endl;
        itkReaderPointer reader = itkReader::New();
        reader->SetFileName(fileNames[i]);
        reader->Update();

        images[i] = reader->GetOutput();
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------



}


#endif
