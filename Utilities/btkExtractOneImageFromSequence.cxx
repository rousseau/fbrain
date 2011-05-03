/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

21 february 2011
< pontabry at unistra dot fr >
23 March 2011
rousseau at unistra dot fr

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
*/


// TCLAP include
#include <tclap/CmdLine.h>

// STL includes
#include "iostream"
#include "string"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"


template<typename T>
int ImageExtractor(std::string inputFileName, std::string outputFileName, unsigned int imageIndex)
{
  //typedef short ScalarType;
  const unsigned int SequenceDimension = 4;
  typedef itk::Image<T,SequenceDimension> Sequence;
  typedef itk::ImageFileReader<Sequence> SequenceFileReader;

  const unsigned int ImageDimension = 3;
  typedef itk::Image<T,ImageDimension> Image;
  typedef itk::ImageFileWriter<Image> ImageFileWriter;

  typedef itk::ExtractImageFilter<Sequence,Image> ExtractImageFilter;

  // Read DWI sequence
  typename SequenceFileReader::Pointer reader = SequenceFileReader::New();
  reader->SetFileName(inputFileName);
  reader->Update();

  // Extract image filter
  typename ExtractImageFilter::Pointer filter = ExtractImageFilter::New();

  filter->SetInput(reader->GetOutput());
  typename Sequence::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
  region.SetSize(3,0);
  region.SetIndex(3,imageIndex);
  filter->SetExtractionRegion(region);
  filter->Update();

  // Write image

  typename ImageFileWriter::Pointer writer = ImageFileWriter::New();
  writer->SetInput(filter->GetOutput());
  writer->SetFileName(outputFileName);

  try
  {
    writer->Update();
    std::cout << "done." << std::endl;
  }
  catch(itk::ImageFileWriterException &err)
  {
    std::cout << "Error: " << std::endl;
    std::cout << err << std::endl;
  }
  return 0;
}




int main(int argc, char *argv[])
{
    // Command line variables
    std::string inputFileName;
    std::string outputFileName;
    unsigned int imageIndex;

    //
    // Parse program arguments
    //

    try
    {
        // Defines command line parser
        TCLAP::CmdLine cmd("BTK DWI baseline image extractor", ' ', "0.1");

        // Defines arguments
        TCLAP::ValueArg<std::string> inputArg("i", "input", "Input DWI sequence", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputArg("o", "output", "Output image", false, "B0.nii.gz", "string", cmd);
        TCLAP::ValueArg<unsigned int> imageIndexArg("", "image_index", "Index of the image to extract", false, 0, "unsigned int", cmd);

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back argument's values
        inputFileName  = inputArg.getValue();
        outputFileName = outputArg.getValue();
        imageIndex     = imageIndexArg.getValue();
    }
    catch(TCLAP::ArgException &e)
    {
        std::cout << "TCLAP error: " << e.error() << " for argument " << e.argId() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    //Guessing the component type
    typedef itk::ImageIOBase::IOComponentType ScalarPixelType;
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFileName.c_str(), itk::ImageIOFactory::ReadMode);
 
    imageIO->SetFileName(inputFileName);
    imageIO->ReadImageInformation();
    const ScalarPixelType pixelType = imageIO->GetComponentType();
    std::cout << "Pixel Type is " << imageIO->GetComponentTypeAsString(pixelType) // 'double'
              << std::endl;

    switch (pixelType){
      case itk::ImageIOBase::SHORT:  ImageExtractor<short>(inputFileName, outputFileName, imageIndex); break;
      case itk::ImageIOBase::USHORT: ImageExtractor<unsigned short>(inputFileName, outputFileName, imageIndex); break;
      case itk::ImageIOBase::FLOAT:  ImageExtractor<float>(inputFileName, outputFileName, imageIndex); break;
      case itk::ImageIOBase::DOUBLE: ImageExtractor<double>(inputFileName, outputFileName, imageIndex); break;

      default:
        std::cerr << "Pixel Type ("
                  << imageIO->GetComponentTypeAsString(pixelType)
                  << ") not supported. Exiting." << std::endl;
        return -1;
    }

    return EXIT_SUCCESS;
}
