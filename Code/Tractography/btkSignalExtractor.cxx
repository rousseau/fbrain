/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

12 april 2010
< pontabry at unistra dot fr >

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty and the software's author, the holder of the
economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and, more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#include "btkSignalExtractor.h"


// STL includes
#include "iostream"
#include "cmath"

#ifndef NDEBUG
    #include "cassert"
#endif // NDEBUG

// ITK includes
#include "itkConstNeighborhoodIterator.h"


namespace btk
{

typedef itk::ConstNeighborhoodIterator<Sequence> SequenceNeighborhoodIterator;

SignalExtractor::SignalExtractor(const std::string &vectorsFileName, const std::string &dataFileName, const std::string &maskFileName, char displayMode)
{
    m_directions = 0;
    m_refIm      = 0;
    m_sigmas     = 0;

    m_displayMode = displayMode;

    this->readFiles(vectorsFileName, dataFileName, maskFileName);
}

SignalExtractor::~SignalExtractor()
{
    delete m_refIm;
}

void SignalExtractor::readFiles(const std::string &vectorsFileName, const std::string &dataFileName, const std::string &maskFileName)
{
    assert(!vectorsFileName.empty());
    assert(!dataFileName.empty());


    Display1(m_displayMode, std::cout << "Reading data files..." << std::endl);


    //
    // Reading dwi sequence
    //

        Display2(m_displayMode, std::cout << "\tReading dwi sequence..." << std::flush);

        SequenceReader::Pointer reader = SequenceReader::New();
        reader->SetFileName(dataFileName);
        reader->Update();               // reading file
        m_data = reader->GetOutput();   // getting data

        Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Reading image mask
    //

        Display2(m_displayMode, std::cout << "\tReading image mask..." << std::flush);

        MaskReader::Pointer mreader = MaskReader::New();
        mreader->SetFileName(maskFileName);
        mreader->Update();
        m_mask = mreader->GetOutput();

        Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Reading vectors file
    //

        Display2(m_displayMode, std::cout << "\tReading vectors..." << std::flush);

        std::vector<Real> values;

        // Open file
        std::fstream vectorsFile(vectorsFileName.c_str(), std::fstream::in);

        if(!vectorsFile.is_open())  // Unable to open file
        {
            std::cout << "Error: unable to read gradient directions file !" << std::endl;
            exit(EXIT_FAILURE);
        }
        else    // file opened
        {
            Real v = 0;

            while((vectorsFile >> v))
                values.push_back(v);

            m_directions = new std::vector<Direction>;
            m_refIm   = new std::vector<unsigned int>;
            unsigned int nbVectors  = values.size()/3;
            unsigned int nbVectors2 = nbVectors + nbVectors;

            for(unsigned int i=0; i<nbVectors; i++)
            {
                Vector tmp(values.at(i), values.at(i+nbVectors), values.at(i+nbVectors2));

                if(tmp.isNull())
                    m_refIm->push_back(i);
                else
                    m_directions->push_back(tmp.toDirection());
            }
        }

        // Close file
        vectorsFile.close();

        Display2(m_displayMode, std::cout << "done." << std::endl);


    Display1(m_displayMode, std::cout << "done." << std::endl);
}

void SignalExtractor::extract()
{
    SequenceRegion dataRegion = m_data->GetLargestPossibleRegion();
    MaskRegion maskRegion     = m_mask->GetLargestPossibleRegion();

    unsigned int nbOfImages = m_directions->size() + m_refIm->size();


    Display1(m_displayMode, std::cout << "Preparing data..." << std::endl);

    //
    // Verify sizes
    //

        Display2(m_displayMode, std::cout << "\tVerifying data..." << std::endl);

        if(dataRegion.GetSize(3) != nbOfImages)
        {
            std::cerr << "Error: bad data !" << std::endl;
            std::cerr << "There are " << dataRegion.GetSize(3) << " MRI images and ";
            std::cerr << nbOfImages << " vectors." << std::endl;
            exit(EXIT_FAILURE);
        }

        if(dataRegion.GetSize(0) != maskRegion.GetSize(0) &&
           dataRegion.GetSize(1) != maskRegion.GetSize(1) &&
           dataRegion.GetSize(2) != maskRegion.GetSize(2))
        {
            std::cerr << "Error: bad data !" << std::endl;
            std::cerr << "Data image is of size ";
            std::cerr << dataRegion.GetSize(0) << "x" << dataRegion.GetSize(1) << "x" << dataRegion.GetSize(2);
            std::cerr << " and image mask is of size ";
            std::cerr << maskRegion.GetSize(0) << "x" << maskRegion.GetSize(1) << "x" << maskRegion.GetSize(2);
            std::cerr << "!" << std::endl;
            exit(EXIT_FAILURE);
        }

        Display2(m_displayMode, std::cout << "\t\tThere are " << m_refIm->size() << " reference images and ");
        Display2(m_displayMode, std::cout << m_directions->size() << " gradient images of size ");
        Display2(m_displayMode, std::cout << dataRegion.GetSize(0) << "x" << dataRegion.GetSize(1) << "x" << dataRegion.GetSize(2) << "." << std::endl);

        Display2(m_displayMode, std::cout << "\tdone." << std::endl);


    //
    // Compute B0 reference image
    //

        Display2(m_displayMode, std::cout << "\tCompute reference image..." << std::flush);

        // Preparating mask iterator
        MaskIterator maskIt(m_mask, maskRegion);

        // Allocate space memory for B0 reference image
        Image::Pointer B0 = Image::New();

        Image::SizeType B0Size;
        B0Size[0] = dataRegion.GetSize(0); B0Size[1] = dataRegion.GetSize(1); B0Size[2] = dataRegion.GetSize(2);
        B0->SetRegions(B0Size);

        Image::SpacingType b0spacing;
        b0spacing[0] = m_data->GetSpacing()[0];
        b0spacing[1] = m_data->GetSpacing()[1];
        b0spacing[2] = m_data->GetSpacing()[2];
        B0->SetSpacing(b0spacing);

        Image::PointType b0origin;
        b0origin[0] = m_data->GetOrigin()[0];
        b0origin[1] = m_data->GetOrigin()[1];
        b0origin[2] = m_data->GetOrigin()[2];
        B0->SetOrigin(b0origin);

        itk::Matrix<Real,3,3> dirMat;
        itk::Matrix<Real,4,4> iniMat = m_data->GetDirection();

        dirMat(0,0) = iniMat(0,0);
        dirMat(0,1) = iniMat(0,1);
        dirMat(0,2) = iniMat(0,2);

        dirMat(1,0) = iniMat(1,0);
        dirMat(1,1) = iniMat(1,1);
        dirMat(1,2) = iniMat(1,2);

        dirMat(2,0) = iniMat(2,0);
        dirMat(2,1) = iniMat(2,1);
        dirMat(2,2) = iniMat(2,2);

        B0->SetDirection(dirMat);

        B0->Allocate();
        B0->FillBuffer(0);

        // Create B0 reference image
        ImageRegion B0Region = B0->GetLargestPossibleRegion();
        ImageIterator B0It(B0, B0Region);

        Sequence::IndexType sIndex;
        Sequence::SizeType sSize;

        std::vector<unsigned int>::iterator refImIt;

        for(refImIt = m_refIm->begin(); refImIt != m_refIm->end(); refImIt++)
        {
            // Set correct region in sequence
            sIndex[0] = 0; sIndex[1] = 0; sIndex[2] = 0;
            sIndex[3] = (*refImIt);
            dataRegion.SetIndex(sIndex);

            sSize[0] = B0Size[0]; sSize[1] = B0Size[1]; sSize[2] = B0Size[2];
            sSize[3] = 1;
            dataRegion.SetSize(sSize);

            SequenceIterator dataIt(m_data, dataRegion);

            // Fill B0 reference image (moy of all reference images)
            for(B0It.GoToBegin(), maskIt.GoToBegin(), dataIt.GoToBegin(); !B0It.IsAtEnd() && !maskIt.IsAtEnd() && !dataIt.IsAtEnd(); ++B0It, ++maskIt, ++dataIt)
                B0It.Set(B0It.Get() + dataIt.Get()/m_refIm->size());
        } // for refImIt

/*
        // Write B0 image
        ImageWriter::Pointer writer = ImageWriter::New();

        writer->SetFileName("B0.nii.gz");
        writer->SetInput(B0);

        try
        {
            writer->Update();
        }
        catch(itk::ImageFileWriterException &err)
        {
            std::cout << "Error: " << std::endl;
            std::cout << err << std::endl;
        }
*/
        Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Compute signal
    //

        Display2(m_displayMode, std::cout << "\tCompute signal..." << std::flush);

        // Allocation of space memory for signal images
        m_signal = Sequence::New();

        Sequence::SizeType signalSize;
        signalSize[0] = dataRegion.GetSize(0); signalSize[1] = dataRegion.GetSize(1);
        signalSize[2] = dataRegion.GetSize(2); signalSize[3] = m_directions->size();
        m_signal->SetRegions(signalSize);

        Sequence::SpacingType spacing = m_data->GetSpacing();
        m_signal->SetSpacing(spacing);

        Sequence::PointType origin = m_data->GetOrigin();
        m_signal->SetOrigin(origin);

        m_signal->SetDirection(m_data->GetDirection());

        m_signal->Allocate();


        // Normalize signal
        unsigned int count = 0;
        SequenceRegion signalRegion;

        for(unsigned int n = 0; n<nbOfImages; n++)
        {
            if(*find(m_refIm->begin(), m_refIm->end()-1, n) != n) // the nth image is not a reference image
            {
                // Define region for data sequence
                sIndex[0] = 0; sIndex[1] = 0; sIndex[2] = 0;
                sIndex[3] = n;
                dataRegion.SetIndex(sIndex);

                sSize[0] = B0Size[0]; sSize[1] = B0Size[1]; sSize[2] = B0Size[2];
                sSize[3] = 1;
                dataRegion.SetSize(sSize);

                SequenceIterator dataIt(m_data, dataRegion);

                // Define region for signal sequence
                sIndex[0] = 0; sIndex[1] = 0; sIndex[2] = 0;
                sIndex[3] = n - count;
                signalRegion.SetIndex(sIndex);

                sSize[0] = B0Size[0]; sSize[1] = B0Size[1]; sSize[2] = B0Size[2];
                sSize[3] = 1;
                signalRegion.SetSize(sSize);

                SequenceIterator signalIt(m_signal, signalRegion);

                // Fill current image (normalization : divide signal by reference image's one)
                for(B0It.GoToBegin(), dataIt.GoToBegin(), signalIt.GoToBegin(), maskIt.GoToBegin();
                        !B0It.IsAtEnd() && !dataIt.IsAtEnd() && !signalIt.IsAtEnd() && !maskIt.IsAtEnd();
                            ++B0It, ++dataIt, ++signalIt, ++maskIt)
                {
                    if(maskIt.Get() != 0)    // is in mask
                        signalIt.Set((B0It.Get() != 0. ? dataIt.Get() / B0It.Get() : 0.));
                    else    // is out of mask
                        signalIt.Set(0);
                }
            }
            else // the nth image is a reference image
            {
                count++;
            }
        } // for n

        Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Compute standard deviation (sigma)
    //

        Display2(m_displayMode, std::cout << "\tCompute standard deviation of signal images..." << std::flush);

        this->computeSigmas();

        Display2(m_displayMode, std::cout << "done." << std::endl);

    Display1(m_displayMode, std::cout << "done." << std::endl);
}

void SignalExtractor::save()
{
    assert(m_directions);
    assert(m_sigmas);


    Display1(m_displayMode, std::cout << "Saving data..." << std::endl);

    //
    // Saving signal
    //

        Display2(m_displayMode, std::cout << "\tsignal into file \"signal.nii.gz\"..." << std::flush);

        SequenceWriter::Pointer writer = SequenceWriter::New();

        writer->SetFileName("signal.nii.gz");
        writer->SetInput(m_signal);

        try
        {
            writer->Update();
            Display2(m_displayMode, std::cout << "done." << std::endl);
        }
        catch(itk::ImageFileWriterException &err)
        {
            std::cout << "Error: " << std::endl;
            std::cout << err << std::endl;
        }


    //
    // Saving directions
    //

        Display2(m_displayMode, std::cout << "\tvectors into file \"directions.txt\"..." << std::flush);

        std::fstream directionsFile("directions.txt", std::fstream::out);

        if(!directionsFile.is_open()) // Unable to open file
        {
            std::cerr << "Error: unable to write in file !" << std::endl;
            exit(EXIT_FAILURE);
        }
        else                       // File is open
        {
            for(std::vector<Direction>::iterator it = m_directions->begin(); it != m_directions->end(); it++)
                directionsFile << (*it).theta() << " " << (*it).phi() << " ";
        }

        directionsFile.close();

        std::cout << "done." << std::endl;


    //
    // Saving standard deviation
    //

        Display2(m_displayMode, std::cout << "\tstandard deviations into file \"sigmas.txt\"..." << std::flush);

        std::fstream sigmasFile("sigmas.txt", std::fstream::out);

        if(!sigmasFile.is_open()) // Unable to open file
        {
            std::cerr << "Error: unable to write in file !" << std::endl;
            exit(EXIT_FAILURE);
        }
        else                       // File is open
        {
            std::vector<Direction>::iterator dIt;

            for(std::vector<Real>::iterator it = m_sigmas->begin(); it != m_sigmas->end(); it++)
                sigmasFile << (*it) << " ";
        }

        sigmasFile.close();

        Display2(m_displayMode, std::cout << "done." << std::endl);

    Display1(m_displayMode, std::cout << "done." << std::endl);


    delete m_sigmas;
    delete m_directions;
}

std::vector<Direction> *SignalExtractor::GetDirections()
{
    return m_directions;
}

Sequence::Pointer SignalExtractor::GetSignal()
{
    return m_signal;
}

std::vector<Real> *SignalExtractor::GetSigmas()
{
    return m_sigmas;
}

Mask::Pointer SignalExtractor::GetMask()
{
    return m_mask;
}

void SignalExtractor::computeSigmas()
{
    // An Optimized Blockwise Non Local Means
    // Denoising Filter for 3D Magnetic Resonance
    // Images by Coupé et al.


    assert(m_signal);
    assert(!m_sigmas);


    m_sigmas = new std::vector<Real>;


    Real _sqrt6on7 = std::sqrt(6./7.);
    Real _1on6     = 1./6.;

    unsigned int xMax = m_signal->GetLargestPossibleRegion().GetSize(0);
    unsigned int yMax = m_signal->GetLargestPossibleRegion().GetSize(1);
    unsigned int zMax = m_signal->GetLargestPossibleRegion().GetSize(2);
    unsigned int kMax = m_signal->GetLargestPossibleRegion().GetSize(3);


    // Define loop region
    SequenceRegion loopRegion;

    Sequence::SizeType loopSize;
    loopSize[0] = xMax; loopSize[1] = yMax;
    loopSize[2] = zMax; loopSize[3] = 1;
    loopRegion.SetSize(loopSize);


    //
    // Compute coefficients
    //

        // Allocate space memory for coefficients image
        Sequence::Pointer coeffs = Sequence::New();
        coeffs->SetRegions(m_signal->GetLargestPossibleRegion());
        coeffs->Allocate();

        // For each image
        for(unsigned int k=0; k<kMax; k++)
        {
            // Set the correct index
            Sequence::IndexType loopIndex;
            loopIndex[0] = 0; loopIndex[1] = 0;
            loopIndex[2] = 0; loopIndex[3] = k;
            loopRegion.SetIndex(loopIndex);

            // Define signal iterator
            SequenceNeighborhoodIterator::RadiusType radius;
            radius[0] = 1; radius[1] = 1; radius[2] = 1; radius[3] = 0;

            SequenceNeighborhoodIterator signalIt(radius, m_signal, loopRegion);

            SequenceNeighborhoodIterator::OffsetType offset1 = { {-1, 0, 0, 0} };
            SequenceNeighborhoodIterator::OffsetType offset2 = { { 1, 0, 0, 0} };
            SequenceNeighborhoodIterator::OffsetType offset3 = { { 0,-1, 0, 0} };
            SequenceNeighborhoodIterator::OffsetType offset4 = { { 0, 1, 0, 0} };
            SequenceNeighborhoodIterator::OffsetType offset5 = { { 0, 0,-1, 0} };
            SequenceNeighborhoodIterator::OffsetType offset6 = { { 0, 0, 1, 0} };

            // Define coefficients images iterator
            SequenceIterator coeffsIt(coeffs, loopRegion);

            // Compute coefficients for each positions in current image
            for(signalIt.GoToBegin(), coeffsIt.GoToBegin(); !signalIt.IsAtEnd() && !coeffsIt.IsAtEnd(); ++signalIt, ++coeffsIt)
            {
                Real sum = signalIt.GetPixel(offset1);
                sum     += signalIt.GetPixel(offset2);
                sum     += signalIt.GetPixel(offset3);
                sum     += signalIt.GetPixel(offset4);
                sum     += signalIt.GetPixel(offset5);
                sum     += signalIt.GetPixel(offset6);

                coeffsIt.Set(_sqrt6on7 * (signalIt.GetCenterPixel() - _1on6*sum));
            } // for signalIt && coeffsIt
        } // for k


    //
    // Compute sigmas
    //

        Real norm = 1. / (xMax * yMax * zMax);

        // For each image
        for(unsigned int k=0; k<kMax; k++)
        {
            // Set the correct index
            Sequence::IndexType loopIndex;
            loopIndex[0] = 0; loopIndex[1] = 0;
            loopIndex[2] = 0; loopIndex[3] = k;
            loopRegion.SetIndex(loopIndex);

            // Define coefficients images iterator
            SequenceIterator coeffsIt(coeffs, loopRegion);

            Real sum = 0;

            // Compute standard deviation (sigma)
            for(coeffsIt.GoToBegin(); !coeffsIt.IsAtEnd(); ++coeffsIt)
                sum += coeffsIt.Get() * coeffsIt.Get();

            m_sigmas->push_back(std::sqrt(norm * sum));
        } // for k
}

} // namespace btk

