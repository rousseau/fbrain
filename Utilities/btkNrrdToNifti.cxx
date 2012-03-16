/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

15 februar 2011
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
#include "tclap/CmdLine.h"

// STL includes
#include "string"
#include "cstdlib"
#include "fstream"

// ITK includes
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkMatrix.h"
#include "vnl/vnl_inverse.h"

// BTK includes
#include "btkNrrdField.h"
#include "btkFileNameTools.h"


// Precision of images' intensities
typedef short PixelType;

// Input sequence (Nrrd -- nhdr)
typedef itk::VectorImage<PixelType,3>           InputSequence;
typedef itk::ImageFileReader<InputSequence>     InputSequenceFileReader;
typedef itk::ImageFileWriter<InputSequence>     InputSequenceFileWriter;
typedef itk::ImageRegionIterator<InputSequence> InputSequenceIterator;

// Output sequence (Nifti -- .nii.gz)
typedef itk::Image<PixelType,4>                  OutputSequence;
typedef itk::ImageFileReader<OutputSequence>     OutputSequenceFileReader;
typedef itk::ImageFileWriter<OutputSequence>     OutputSequenceFileWriter;
typedef itk::ImageRegionIterator<OutputSequence> OutputSequenceIterator;

// Anatomical volume (3D)
typedef itk::Image<PixelType,3>         Image;
typedef itk::ImageFileReader<Image>     ImageFileReader;
typedef itk::ImageFileWriter<Image>     ImageFileWriter;
typedef itk::ImageRegionIterator<Image> ImageIterator;


int main(int argc, char *argv[])
{
    //
    // Parse program's arguments
    //

    // Define command line parameters
    std::string inputFileName;
    std::string outputFileName;
    bool dwiConversion;


    // Define command line parser
    TCLAP::CmdLine cmd("Nrrd to Nifti Converter", ' ', "0.1");

    // Define command line arguments
    TCLAP::ValueArg<std::string> inputArg("i", "input", "Input image", true, "", "string", cmd);
    TCLAP::ValueArg<std::string> outputArg("o", "output", "Output image (default: \"name.nii.gz\" -- input name with nifti extension)", false, "", "string", cmd);

    TCLAP::SwitchArg dwiConversionArg("", "dwi", "Proceed to a conversion diffusion weighted image conversion (default: false -- anatomical 3D volume)", cmd, false);

    // Parse arguments
    cmd.parse(argc, argv);

    // Get back arguments' values
    inputFileName   = inputArg.getValue();
    outputFileName  = outputArg.getValue();
    dwiConversion   = dwiConversionArg.getValue();


    // If no filename is given for output, set it up with input name
    if(outputFileName.empty())
    {
        outputFileName = btk::GetRadixOf(inputFileName);
    }
    else // !outputFileName.empty()
    {
        outputFileName = btk::GetRadixOf(outputFileName);
    }


    //
    // Main part
    //

    try
    {
        //
        // Verifications before processing
        //

        // Only .nhdr and .nrrd file formats are supported
        std::string inputFormat = btk::GetExtensionOf(inputFileName);

        if(inputFormat != ".nhdr" && inputFormat != ".nrrd")
            throw std::string("Only NRRD file formats, such as .nhdr/.raw(.gz) and .nrrd, are supported !");


        //
        // Processing
        //

        if(dwiConversion)
        {
            //
            // Read input image (Nrrd file format)
            //

            std::cout << "Reading file \"" << inputFileName << "\"... " << std::flush;

            // Read raw data
            InputSequenceFileReader::Pointer reader = InputSequenceFileReader::New();
            reader->SetFileName(inputFileName.c_str());
            reader->Update();

            // Read header data
            std::fstream headerFile(inputFileName.c_str(), std::fstream::in);

            std::string token;
            std::string key;
            std::string value;
            unsigned int bvalue = 0;
            bool vectorFound = false;

            char buffer[256];
            headerFile.getline(buffer, 256);
            token = buffer;

            while(!vectorFound && !headerFile.eof() && !headerFile.fail() && !headerFile.bad() && !token.empty())
            {
                btk::btkNrrdField currentField(token);
                key   = currentField.GetKey();
                value = currentField.GetValue();

                if(key == "DWMRI_b-value")
                {
                    std::stringstream stream;
                    stream << value;
                    stream >> bvalue;
                }
                else if(key == "DWMRI_gradient_0000")
                {
                    vectorFound = true;
                }

                headerFile.getline(buffer,256);
                token = buffer;
            } // while file is not empty

            std::vector<double> vx;
            std::vector<double> vy;
            std::vector<double> vz;
            std::vector<unsigned int> bvals;

            if(vectorFound == true)
            {
                double x,y,z;

                std::stringstream stream;
                stream << value;
                stream >> x; vx.push_back(x);
                stream >> y; vy.push_back(y);
                stream >> z; vz.push_back(z);

                if(x == 0 && y == 0 && z == 0)
                    bvals.push_back(0);
                else
                    bvals.push_back(bvalue);

                while(!headerFile.eof() && !headerFile.fail() && !headerFile.bad() && !token.empty())
                {
                    btk::btkNrrdField field(token);
                    key   = field.GetKey();
                    value = field.GetValue();

                    stream.clear();
                    stream << value;
                    stream >> x; vx.push_back(x);
                    stream >> y; vy.push_back(y);
                    stream >> z; vz.push_back(z);

                    if(x == 0 && y == 0 && z == 0)
                        bvals.push_back(0);
                    else
                        bvals.push_back(bvalue);

                    headerFile.getline(buffer,256);
                    token = buffer;
                } // while file is not empty
            } // if vector found

            headerFile.close();

            std::cout << "done." << std::endl;


            //
            // Convert file
            //

            std::cout << "Converting file... " << std::flush;

            // Get input image's properties
            InputSequence::RegionType inputRegion = reader->GetOutput()->GetLargestPossibleRegion();
            unsigned int componentsNumber         = reader->GetOutput()->GetVectorLength();

            InputSequence::SizeType inputSize = inputRegion.GetSize();
            OutputSequence::SizeType outputSize;
            outputSize[0] = inputSize[0];
            outputSize[1] = inputSize[1];
            outputSize[2] = inputSize[2];
            outputSize[3] = componentsNumber;

            InputSequence::PointType inputOrigin = reader->GetOutput()->GetOrigin();
            OutputSequence::PointType outputOrigin;
            outputOrigin[0] = inputOrigin[0];
            outputOrigin[1] = inputOrigin[1];
            outputOrigin[2] = inputOrigin[2];
            outputOrigin[3] = 0;

            InputSequence::SpacingType inputSpacing = reader->GetOutput()->GetSpacing();
            OutputSequence::SpacingType outputSpacing;
            outputSpacing[0] = inputSpacing[0];
            outputSpacing[1] = inputSpacing[1];
            outputSpacing[2] = inputSpacing[2];
            outputSpacing[3] = 1;

            itk::Matrix<double,3,3> inputDirection = reader->GetOutput()->GetDirection();
            itk::Matrix<double,4,4> outputDirection;
            outputDirection(0,0) = inputDirection(0,0); outputDirection(0,1) = inputDirection(0,1); outputDirection(0,2) = inputDirection(0,2); outputDirection(0,3) = 0;
            outputDirection(1,0) = inputDirection(1,0); outputDirection(1,1) = inputDirection(1,1); outputDirection(1,2) = inputDirection(1,2); outputDirection(1,3) = 0;
            outputDirection(2,0) = inputDirection(2,0); outputDirection(2,1) = inputDirection(2,1); outputDirection(2,2) = inputDirection(2,2); outputDirection(2,3) = 0;
            outputDirection(3,0) = 0;                   outputDirection(3,1) = 0;                   outputDirection(3,2) = 0;                   outputDirection(3,3) = 1;

            // Create and initialize a new sequence with the same properties as the original
            OutputSequence::Pointer newSequence = OutputSequence::New();
            newSequence->SetRegions(outputSize);
            newSequence->SetOrigin(outputOrigin);
            newSequence->SetSpacing(outputSpacing);
            newSequence->SetDirection(outputDirection);
            newSequence->Allocate();
            newSequence->FillBuffer(0);


            // Conversion loop
            InputSequenceIterator inputIt(reader->GetOutput(), inputRegion);

            for(unsigned int k=0; k<componentsNumber; k++)
            {
                // Set region in output sequence for current gradient image
                OutputSequence::IndexType outputRegionCorner;
                outputRegionCorner[0] = 0; outputRegionCorner[1] = 0; outputRegionCorner[2] = 0; outputRegionCorner[3] = k;

                OutputSequence::SizeType outputRegionSize = outputSize;
                outputRegionSize[3] = 1;

                OutputSequence::RegionType outputRegion;
                outputRegion.SetIndex(outputRegionCorner);
                outputRegion.SetSize(outputRegionSize);

                // Set region iterator for current gradient image
                OutputSequenceIterator outputIt(newSequence, outputRegion);

                // Loop
                for(inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd() && !outputIt.IsAtEnd(); ++inputIt, ++outputIt)
                {
                    outputIt.Set(inputIt.Get()[k]);
                } // for each voxel
            } // for each component


//////////////// FIXME : à corriger par une version ITK
//            // Conversion loop for gradient vectors (world coordinates to image coordinates)
//            vnl_vector<double> worldCoord(3);
//            vnl_vector<double> imgCoord(3);
//            vnl_matrix<double> wcToImg = vnl_inverse(inputDirection.GetVnlMatrix());
//            std::vector<double>::iterator ivx;
//            std::vector<double>::iterator ivy;
//            std::vector<double>::iterator ivz;
//            for(ivx=vx.begin(), ivy=vy.begin(), ivz=vz.begin(); ivx!=vx.end() && ivy!=vy.end() && ivz!=vz.end(); ivx++, ivy++, ivz++)
//            {
//                worldCoord(0) = -(*ivx); worldCoord(1) = -(*ivy); worldCoord(2) = *ivz;
//                imgCoord = wcToImg * worldCoord;
//                *ivx = imgCoord(0); *ivy = imgCoord(1); *ivz = imgCoord(2);
//            }
///////////////
            std::cout << "done." << std::endl;


            //
            // Write output image (Nifti file format)
            //

            std::cout << "Writing files \"" << outputFileName << "{.nii.gz|.bval|.bvec}\"... " << std::flush;

            // Write raw data
            OutputSequenceFileWriter::Pointer writer = OutputSequenceFileWriter::New();
            writer->SetFileName(outputFileName+".nii.gz");
            writer->SetInput(newSequence);
            writer->Update();

            // Write gradient table of vectors
            std::fstream bvecFile((outputFileName+".bvec").c_str(), std::fstream::out);

            for(std::vector<double>::iterator it=vx.begin(); it != vx.end(); it++)
                bvecFile << *it << " ";

            bvecFile << std::endl;

            for(std::vector<double>::iterator it=vy.begin(); it != vy.end(); it++)
                bvecFile << *it << " ";

            bvecFile << std::endl;

            for(std::vector<double>::iterator it=vz.begin(); it != vz.end(); it++)
                bvecFile << *it << " ";

            bvecFile << std::endl;

            bvecFile.close();

            // Write b-values
            std::fstream bvalFile((outputFileName+".bval").c_str(), std::fstream::out);

            for(std::vector<unsigned int>::iterator it=bvals.begin(); it != bvals.end(); it++)
                bvalFile << *it << " ";

            bvalFile.close();

            std::cout << "done." << std::endl;
        }
        else // !dwiConversion
        {
            //
            // Read input image (Nrrd file format)
            //

            std::cout << "Reading file \"" << inputFileName << "\"..." << std::flush;

            ImageFileReader::Pointer reader = ImageFileReader::New();
            reader->SetFileName(inputFileName.c_str());
            reader->Update();

            std::cout << "done." << std::endl;


            //
            // Convert file
            //

            std::cout << "Converting file..." << std::flush;

            // Create and initialize a new image with the same properties as the original
            Image::RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
            Image::Pointer newImage  = Image::New();
            newImage->SetRegions(region.GetSize());
            newImage->SetOrigin(reader->GetOutput()->GetOrigin());
            newImage->SetSpacing(reader->GetOutput()->GetSpacing());
            newImage->SetDirection(reader->GetOutput()->GetDirection());
            newImage->Allocate();
            newImage->FillBuffer(0);

            // Conversion loop
            ImageIterator inputIt(reader->GetOutput(), region);
            ImageIterator outputIt(newImage, newImage->GetLargestPossibleRegion());

            for(inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd() && !outputIt.IsAtEnd(); ++inputIt, ++outputIt)
            {
                outputIt.Set(inputIt.Get());
            } // for each voxel

            std::cout << "done." << std::endl;


            //
            // Write output image (Nifti file format)
            //

            std::cout << "Writing file \"" << outputFileName << ".nii.gz\"..." << std::flush;

            ImageFileWriter::Pointer writer = ImageFileWriter::New();
            writer->SetFileName(outputFileName+".nii.gz");
            writer->SetInput(reader->GetOutput());
            writer->Update();

            std::cout << "done." << std::endl;
        } // else !dwiConversion
    }
    catch(itk::ExceptionObject &err)
    {
        std::cout << "Error: " << err << std::endl;
    }
    catch(std::string &message)
    {
        std::cout << "Error: " << message << std::endl;
    }

    return EXIT_SUCCESS;
}
