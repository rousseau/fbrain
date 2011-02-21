/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

21 februar 2011
< pontabry at unistra dot fr >

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
#include "CmdLine.h"

// STL includes
#include "iostream"
#include "string"

// ITK includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"


typedef double ScalarType;

const unsigned int ImageDimension = 3;
typedef itk::Image<ScalarType,ImageDimension> Image;
typedef itk::ImageFileReader<Image> ImageFileReader;
typedef itk::ImageFileWriter<Image> ImageFileWriter;


int main(int argc, char *argv[])
{
    // Command line variables
    std::string inputFileName;
    std::string outputFileName;


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

        // Parsing arguments
        cmd.parse(argc, argv);

        // Get back argument's values
        inputFileName  = inputArg.getValue();
        outputFileName = outputArg.getValue();
    }
    catch(TCLAP::ArgException &e)
    {
        std::cout << "TCLAP error: " << e.error() << " for argument " << e.argId() << std::endl;
        std::exit(EXIT_FAILURE);
    }


    // Read DWI sequence

    std::cout << "Reading DWI sequence..." << std::flush;

    ImageFileReader::Pointer reader = ImageFileReader::New();
    reader->SetFileName(inputFileName);
    reader->Update();

    std::cout << "done." << std::endl;


    // Write image

    std::cout << "Writing baseline image..." << std::flush;

    ImageFileWriter::Pointer writer = ImageFileWriter::New();
    writer->SetInput(reader->GetOutput());
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


    return EXIT_SUCCESS;
}
