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
#include "CmdLine.h"

// STL includes
#include "string"
#include "cstdlib"
#include "fstream"

// ITK includes
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkMatrix.h"


typedef double PixelType;
const unsigned int InDimension  = 3;
const unsigned int OutDimension = 4;

typedef itk::VectorImage<PixelType,InDimension> InSequence;
typedef itk::ImageFileReader<InSequence> InSequenceFileReader;
typedef itk::ImageFileWriter<InSequence> InSequenceFileWriter;
typedef itk::ImageRegionIteratorWithIndex<InSequence> InSequenceIterator;

typedef itk::Image<PixelType,OutDimension> OutSequence;
typedef itk::ImageFileReader<OutSequence> OutSequenceFileReader;
typedef itk::ImageFileWriter<OutSequence> OutSequenceFileWriter;

const unsigned int Dimension = 3;
typedef itk::Image<PixelType,Dimension> Image;
typedef itk::ImageFileReader<Image> ImageFileReader;
typedef itk::ImageFileWriter<Image> ImageFileWriter;
typedef itk::ImageRegionIterator<Image> ImageIterator;


int main(int argc, char *argv[])
{
    //
    // Parse program's arguments
    //

    // Define command line variables
    std::string inFileName;
    std::string outFileName;
    std::string vecFileName;
    bool dwi;

    // Define command line parser
    TCLAP::CmdLine cmd("Nrrd to Nifti Converter", ' ', "0.1");

    // Define command line arguments
    TCLAP::ValueArg<std::string> inArg("i", "input", "Input image", true, "", "string", cmd);
    TCLAP::ValueArg<std::string> outArg("o", "output", "Output image", false, "data.nii.gz", "string", cmd);
    TCLAP::ValueArg<std::string> vecArg("v", "vector", "Output vector", false, "data.bvec", "string", cmd);

    TCLAP::SwitchArg dwiArg("d", "dwi", "DWI image conversion", cmd, false);

    // Parse arguments
    cmd.parse(argc, argv);

    // Get back arguments' values
    inFileName   = inArg.getValue();
    outFileName  = outArg.getValue();
    vecFileName  = vecArg.getValue();
    dwi          = dwiArg.getValue();


    if(dwi)
    {
        //
        // Read input image (Nrrd file format)
        //

        std::cout << "Reading file..." << std::flush;

        InSequenceFileReader::Pointer reader = InSequenceFileReader::New();
        reader->SetFileName(inFileName.c_str());
        reader->Update();


        std::fstream headFile(inFileName.c_str(), std::fstream::in);

        std::string s;
        bool stop = false;

        while(!stop && (headFile >> s))
        {
            if(s == "DWMRI_gradient_0000:=")
                stop = true;
        }

        std::vector<double> vx;
        std::vector<double> vy;
        std::vector<double> vz;

        if(stop == true)
        {
            double x,y,z;

            while((headFile >> x) && (headFile >> y) && (headFile >> z))
            {
                vx.push_back(x);
                vy.push_back(y);
                vz.push_back(z);

                headFile >> s;
            }

        }

        headFile.close();

        std::cout << "done." << std::endl;


        //
        // Convert file
        //

        std::cout << "Converting file..." << std::endl;

        // Image properties
        OutSequence::SizeType size;
        size[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(0);
        size[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(1);
        size[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2);
        size[3] = reader->GetOutput()->GetVectorLength();

        std::cout << "\tsize: " << size[0] << "x" << size[1] << "x" << size[2] << "x" << size[3] << std::endl;

        OutSequence::PointType origin;
        origin[0] = reader->GetOutput()->GetOrigin()[0];
        origin[1] = reader->GetOutput()->GetOrigin()[1];
        origin[2] = reader->GetOutput()->GetOrigin()[2];
        origin[3] = 0;

        std::cout << "\torigin = (" << origin[0] << "," << origin[1] << "," << origin[2] << ")" << std::endl;

        OutSequence::SpacingType spacing;
        spacing[0] = reader->GetOutput()->GetSpacing()[0];
        spacing[1] = reader->GetOutput()->GetSpacing()[1];
        spacing[2] = reader->GetOutput()->GetSpacing()[2];
        spacing[3] = 1;

        std::cout << "\tspacing: " << spacing[0] << "x" << spacing[1] << "x" << spacing[2] << "x" << spacing[3] << std::endl;

        itk::Matrix<double,4,4> dirMat;
        itk::Matrix<double,3,3> iniMat = reader->GetOutput()->GetDirection();

        dirMat(0,0) = (iniMat(0,0) != 0 ? -iniMat(0,0) : iniMat(0,0));
        dirMat(0,1) = (iniMat(0,1) != 0 ? -iniMat(0,1) : iniMat(0,1));
        dirMat(0,2) = (iniMat(0,2) != 0 ? -iniMat(0,2) : iniMat(0,2));
        dirMat(0,3) = 0;

        dirMat(1,0) = (iniMat(1,0) != 0 ? -iniMat(1,0) : iniMat(1,0));
        dirMat(1,1) = (iniMat(1,1) != 0 ? -iniMat(1,1) : iniMat(1,1));
        dirMat(1,2) = (iniMat(1,2) != 0 ? -iniMat(1,2) : iniMat(1,2));
        dirMat(1,3) = 0;

        dirMat(2,0) = iniMat(2,0);
        dirMat(2,1) = iniMat(2,1);
        dirMat(2,2) = iniMat(2,2);
        dirMat(2,3) = 0;

        dirMat(3,0) = 0;
        dirMat(3,1) = 0;
        dirMat(3,2) = 0;
        dirMat(3,3) = 1;

        std::cout << "\torientation: (" << dirMat(0,0) << "," << dirMat(0,1) << "," << dirMat(0,2) << "), (" << dirMat(1,0) << "," << dirMat(1,1) << "," << dirMat(1,2) << "), (" << dirMat(2,0) << "," << dirMat(2,1) << "," << dirMat(2,2) << ")" << std::endl;


        // Temporary image file
        OutSequence::Pointer image = OutSequence::New();
        image->SetRegions(size);
        image->SetOrigin(origin);
        image->SetSpacing(spacing);
        image->SetDirection(dirMat);
        image->Allocate();
        image->FillBuffer(0);

        // Iterators
        InSequenceIterator in(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

        // Conversion loop
        for(in.GoToBegin(); !in.IsAtEnd(); ++in)
        {
            for(unsigned int k=0; k<reader->GetOutput()->GetVectorLength(); k++)
            {
                InSequence::IndexType  inIndex = in.GetIndex();

                OutSequence::IndexType outIndex;
                outIndex[0] = inIndex[0]; outIndex[1] = inIndex[1]; outIndex[2] = inIndex[2]; outIndex[3] = k;

                image->SetPixel(outIndex, in.Get()[k]);
            } // for each scalar
        } // for each vector

        std::cout << "done." << std::endl;


        //
        // Write output image (Nifti file format)
        //

        std::cout << "Writing file..." << std::flush;

        OutSequenceFileWriter::Pointer writer = OutSequenceFileWriter::New();
        writer->SetFileName(outFileName);
        writer->SetInput(image);

        try
        {
            writer->Update();

            std::fstream bvecFile(vecFileName.c_str(), std::fstream::out);

            for(std::vector<double>::iterator it=vx.begin(); it != vx.end(); it++)
                bvecFile << *it << " ";

            bvecFile << std::endl;

            for(std::vector<double>::iterator it=vy.begin(); it != vy.end(); it++)
                bvecFile << *it << " ";

            bvecFile << std::endl;

            for(std::vector<double>::iterator it=vy.begin(); it != vy.end(); it++)
                bvecFile << *it << " ";

            bvecFile << std::endl;

            bvecFile.close();

            std::cout << "done." << std::endl;
        }
        catch(itk::ExceptionObject &err)
        {
            std::cout << "Error: " << err << std::endl;
        }
    }
    else // not dwi image
    {
        //
        // Read input image (Nrrd file format)
        //

        std::cout << "Reading file..." << std::flush;

        ImageFileReader::Pointer reader = ImageFileReader::New();
        reader->SetFileName(inFileName.c_str());
        reader->Update();

        std::cout << "done." << std::endl;


        //
        // Convert file
        //

        std::cout << "Converting file..." << std::endl;

        // Image properties
        Image::SizeType size;
        size[0] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(0);
        size[1] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(1);
        size[2] = reader->GetOutput()->GetLargestPossibleRegion().GetSize(2);

        std::cout << "\tsize: " << size[0] << "x" << size[1] << "x" << size[2] << std::endl;

        Image::PointType origin;
        origin[0] = reader->GetOutput()->GetOrigin()[0];
        origin[1] = reader->GetOutput()->GetOrigin()[1];
        origin[2] = reader->GetOutput()->GetOrigin()[2];

        std::cout << "\torigin = (" << origin[0] << "," << origin[1] << "," << origin[2] << ")" << std::endl;

        Image::SpacingType spacing;
        spacing[0] = reader->GetOutput()->GetSpacing()[0];
        spacing[1] = reader->GetOutput()->GetSpacing()[1];
        spacing[2] = reader->GetOutput()->GetSpacing()[2];

        std::cout << "\tspacing: " << spacing[0] << "x" << spacing[1] << "x" << spacing[2] << std::endl;

        itk::Matrix<double,4,4> dirMat;
        itk::Matrix<double,3,3> iniMat = reader->GetOutput()->GetDirection();

        dirMat(0,0) = iniMat(0,0);
        dirMat(0,1) = iniMat(0,1);
        dirMat(0,2) = iniMat(0,2);
        dirMat(0,3) = 0;

        dirMat(1,0) = iniMat(1,0);
        dirMat(1,1) = iniMat(1,1);
        dirMat(1,2) = iniMat(1,2);
        dirMat(1,3) = 0;

        dirMat(2,0) = iniMat(2,0);
        dirMat(2,1) = iniMat(2,1);
        dirMat(2,2) = iniMat(2,2);
        dirMat(2,3) = 0;

        dirMat(3,0) = 0;
        dirMat(3,1) = 0;
        dirMat(3,2) = 0;
        dirMat(3,3) = 1;

        std::cout << "\torientation: (" << dirMat(0,0) << "," << dirMat(0,1) << "," << dirMat(0,2) << "), (" << dirMat(1,0) << "," << dirMat(1,1) << "," << dirMat(1,2) << "), (" << dirMat(2,0) << "," << dirMat(2,1) << "," << dirMat(2,2) << ")" << std::endl;

        // Temporary image file
        Image::Pointer image = Image::New();
        image->SetRegions(size);
        image->SetOrigin(origin);
        image->SetSpacing(spacing);
        image->Allocate();
        image->FillBuffer(0);

        // Iterators
        ImageIterator in(reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());
        ImageIterator out(image, image->GetLargestPossibleRegion());

        // Conversion loop
        for(in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
        {
            out.Set(in.Get());
        } // for each vector

        std::cout << "done." << std::endl;


        //
        // Write output image (Nifti file format)
        //

        std::cout << "Writing file..." << std::flush;

        ImageFileWriter::Pointer writer = ImageFileWriter::New();
        writer->SetFileName(outFileName);
        writer->SetInput(reader->GetOutput());

        try
        {
            writer->Update();

            std::cout << "done." << std::endl;
        }
        catch(itk::ExceptionObject &err)
        {
            std::cout << "Error: " << err << std::endl;
        }
    }



    return EXIT_SUCCESS;
}
