/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

24 februar 2010
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


#include "btkSHModelEstimator.h"

#ifndef NDEBUG
    #include "cassert"
#endif // NDEBUG

// STL includes
#include "cmath"

// Local includes
#include "btkSphericalHarmonics.h"


namespace btk
{

SHModelEstimator::SHModelEstimator(const std::string &signalFileName, const std::string &directionsFileName, const std::string &maskFileName, unsigned int order, Real lambda)
{
    m_directions = 0;
    m_Y          = 0;
    m_L          = 0;
    m_Omega      = 0;

    m_order  = order;
    m_R      = 0.5 * (m_order+1) * (m_order+2);
    m_lambda = lambda;

    m_2PI = 2. * M_PI;

    this->readFiles(signalFileName, directionsFileName, maskFileName);
}

SHModelEstimator::SHModelEstimator(Sequence::Pointer signal, std::vector<Direction> *directions, Mask::Pointer mask, unsigned int order, Real lambda, char displayMode)
{
    m_directions = directions;
    m_signal     = signal;
    m_mask       = mask;
    m_Y          = 0;
    m_L          = 0;
    m_Omega      = 0;

    m_order  = order;
    m_R      = 0.5 * (m_order+1) * (m_order+2);
    m_lambda = lambda;

    m_2PI = 2. * M_PI;

    m_displayMode = displayMode;
}

SHModelEstimator::~SHModelEstimator()
{
    delete m_Y;
    delete m_L;
    delete m_Omega;
}

void SHModelEstimator::estimate()
{
    SequenceRegion signalRegion = m_signal->GetLargestPossibleRegion();
    MaskRegion maskRegion       = m_mask->GetLargestPossibleRegion();


    Display1(m_displayMode, std::cout << "Estimating model..." << std::endl);


    //
    // Verify data
    //

        Display2(m_displayMode, std::cout << "\tVerifying data..." << std::endl);

        if(m_directions->size() != signalRegion.GetSize(3)) // there is a problem
        {
            std::cerr << "Error: ";
            std::cerr << "There are " << m_directions->size() << " gradient directions and ";
            std::cerr << signalRegion.GetSize(3) << " gradient images !" << std::endl;
            exit(EXIT_FAILURE);
        }

        if(signalRegion.GetSize(0) != maskRegion.GetSize(0) &&
           signalRegion.GetSize(1) != maskRegion.GetSize(1) &&
           signalRegion.GetSize(2) != maskRegion.GetSize(2))
        {
            std::cerr << "Error: bad data !" << std::endl;
            std::cerr << "Data image is of size ";
            std::cerr << signalRegion.GetSize(0) << "x" << signalRegion.GetSize(1) << "x" << signalRegion.GetSize(2);
            std::cerr << " and image mask is of size ";
            std::cerr << maskRegion.GetSize(0) << "x" << maskRegion.GetSize(1) << "x" << maskRegion.GetSize(2);
            std::cerr << "!" << std::endl;
            exit(EXIT_FAILURE);
        }

        if(m_R <= 0) // bad order
        {
            std::cerr << "Error: ";
            std::cerr << "Bad order (1/2 * (order+1) * (order+2) = " << m_R << " <= 0) !" << std::endl;
            exit(EXIT_FAILURE);
        }

        Display2(m_displayMode, std::cout << "\t\tThere are " << m_directions->size() << " gradient directions and images." << std::endl);
        Display2(m_displayMode, std::cout << "\tdone." << std::endl);


    //
    // Allocate memory for image
    //

        Display2(m_displayMode, std::cout << "\tAllocating memory space for model image..." << std::flush);

        m_model = Sequence::New();

        Sequence::SizeType modelSize;
        modelSize[0] = signalRegion.GetSize(0); modelSize[1] = signalRegion.GetSize(1);
        modelSize[2] = signalRegion.GetSize(2); modelSize[3] = m_R;
        m_model->SetRegions(modelSize);

        Sequence::SpacingType spacing = m_signal->GetSpacing();
        m_model->SetSpacing(spacing);

        Sequence::PointType origin = m_signal->GetOrigin();
        m_model->SetOrigin(origin);

        m_model->SetDirection(m_signal->GetDirection());

        m_model->Allocate();
        m_model->FillBuffer(0);

        Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    //  Compute needed matricies
    //

        Display2(m_displayMode, std::cout << "\tComputing needed matrices..." << std::endl);

        Display2(m_displayMode, std::cout << "\t\tSpherical harmonics basis matrix..." << std::flush);
        this->computeSHBasisMatrix();
        Display2(m_displayMode, std::cout << "done." << std::endl);

        Display2(m_displayMode, std::cout << "\t\tLaplace-Beltrami smoothing matrix..." << std::flush);
        this->computeLaplaceBeltramiMatrix();
        Display2(m_displayMode, std::cout << "done." << std::endl);

        Display2(m_displayMode, std::cout << "\t\tTransition matrix..." << std::flush);
        this->computeOmegaMatrix();
        Display2(m_displayMode, std::cout << "done." << std::endl);

        Display2(m_displayMode, std::cout << "\tdone." << std::endl);


    //
    //  Compute model image
    //

        Display2(m_displayMode, std::cout << "\tComputing model..." << std::flush);

        unsigned int xMax = signalRegion.GetSize(0);
        unsigned int yMax = signalRegion.GetSize(1);
        unsigned int zMax = signalRegion.GetSize(2);

        Sequence::SizeType signalSize;
        Sequence::IndexType signalIndex;
        Sequence::IndexType modelIndex;

        SequenceRegion modelRegion;

        Matrix &Omega = *m_Omega;

        for(unsigned int z=0; z<zMax; z++)
        {
            for(unsigned int x=0; x<xMax; x++)
            {
                for(unsigned int y=0; y<yMax; y++)
                {
                    Mask::IndexType mindex;
                    mindex[0] = x; mindex[1] = y; mindex[2] = z;

                    if(m_mask->GetPixel(mindex) != 0) // is in mask
                    {
                        // Define region for signal
                        signalIndex[0] = x; signalIndex[1] = y;
                        signalIndex[2] = z; signalIndex[3] = 0;
                        signalRegion.SetIndex(signalIndex);

                        signalSize[0] = 1; signalSize[1] = 1;
                        signalSize[2] = 1; signalSize[3] = m_directions->size();
                        signalRegion.SetSize(signalSize);

                        SequenceIterator signalIt(m_signal, signalRegion);


                        // Define region for model
                        modelIndex[0] = x; modelIndex[1] = y;
                        modelIndex[2] = z; modelIndex[3] = 0;
                        modelRegion.SetIndex(modelIndex);

                        modelSize[0] = 1; modelSize[1] = 1;
                        modelSize[2] = 1; modelSize[3] = m_R;
                        modelRegion.SetSize(modelSize);

                        SequenceIterator modelIt(m_model, modelRegion);


                        // Get S_p matrix
                        Matrix Sp(m_directions->size(), 1);
                        unsigned int i = 0;

                        for(signalIt.GoToBegin(); !signalIt.IsAtEnd(); ++signalIt)
                            Sp(i++,0) = signalIt.Get();


                        // Compute C_p matrix
                        Matrix Cp(m_R, 1);

                        Cp = Omega * Sp;


                        // Set model image
                        i = 0;

                        for(modelIt.GoToBegin(); !modelIt.IsAtEnd(); ++modelIt)
                            modelIt.Set(Cp(i++,0));
                    } // if in mask
                } // for y
            } // for x
        } // for z

        Display2(m_displayMode, std::cout << "done." << std::endl);


    Display1(m_displayMode, std::cout << "done." << std::endl);
}

void SHModelEstimator::save()
{
    std::cout << "Saving data..." << std::endl;

    //
    // Model file
    //

        std::cout << "\tmodel into file \"model.nii.gz\"..." << std::flush;

        SequenceWriter::Pointer writer = SequenceWriter::New();

        writer->SetFileName("model.nii.gz");
        writer->SetInput(m_model);

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


    std::cout << "done." << std::endl;

    delete m_directions;
}

Sequence::Pointer SHModelEstimator::GetModel()
{
    return m_model;
}

void SHModelEstimator::readFiles(const std::string &signalFileName, const std::string &directionsFileName, const std::string &maskFileName)
{
    assert(!directionsFileName.empty());
    assert(!signalFileName.empty());


    Display1(m_displayMode, std::cout << "Reading data files..." << std::endl);


    //
    // Reading signal
    //

        Display2(m_displayMode, std::cout << "\tReading signal image..." << std::flush);

        SequenceReader::Pointer reader = SequenceReader::New();
        reader->SetFileName(signalFileName);
        reader->Update();               // reading file
        m_signal = reader->GetOutput(); // getting data

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
    // Reading gradients directions file
    //

        Display2(m_displayMode, std::cout << "\tReading gradient directions..." << std::flush);

        // Open file
        std::fstream directionsFile(directionsFileName.c_str(), std::fstream::in);

        if(!directionsFile.is_open())  // Unable to open file
        {
            std::cout << "Error: unable to read gradient directions file !" << std::endl;
            exit(EXIT_FAILURE);
        }
        else                           // File is open
        {
            m_directions = new std::vector<Direction>;

            Real theta = 0;
            Real phi   = 0;

            while((directionsFile >> theta) && (directionsFile >> phi))
                m_directions->push_back(Direction(theta,phi));
        }

        // Close file
        directionsFile.close();

        Display2(m_displayMode, std::cout << "done." << std::endl);


    Display1(m_displayMode, std::cout << "done." << std::endl);
}

void SHModelEstimator::computeSHBasisMatrix()
{
    assert(m_directions);
    assert(!m_Y);

    // Allocate memory space for Y matrix
    m_Y = new Matrix(m_directions->size(), m_R);

    Matrix &Y = *m_Y;

    // Compute Y matrix
    Real _4PI  = 4.0 * M_PI;


    for(unsigned int u=0; u<m_directions->size(); u++)
    {
        unsigned int j = 0;

        for(unsigned int l=0; l<=m_order; l+=2)
        {
            Real coef = std::sqrt((2.0*(Real)l + 1.0) / _4PI);

            for(int m=-(int)l; m<=(int)l; m++)
                Y(u,j++) = coef * SphericalHarmonics::computeBasis(m_directions->at(u), l, m);
        } // for l
    } // for u

    #ifndef NDEBUG
        std::cerr << "Y (" << Y.Rows() << "x" << Y.Cols() << ")" << std::endl;
        Y.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

void SHModelEstimator::computeLaplaceBeltramiMatrix()
{
    assert(!m_L);

    m_L = new Matrix(m_R, m_R);

    Matrix &L = *m_L;

    L.Fill(0);

    unsigned int i = 0;

    for(unsigned int l=0; l<=m_order; l+=2)
    {
        Real lp1   = l+1;
        Real value = l*l * lp1*lp1;

        for(int m=-(int)l; m<=(int)l; m++)
        {
            L(i,i) = value;
            i++;
        } // for m
    } // for l

    #ifndef NDEBUG
        std::cerr << "L (" << L.Rows() << "x" << L.Cols() << ")" << std::endl;
        L.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

void SHModelEstimator::computeOmegaMatrix()
{
    assert(m_directions);
    assert(m_Y);
    assert(m_L);
    assert(!m_Omega);

    // Allocate memory space for Omega matrix
    m_Omega = new Matrix(m_R, m_directions->size());

    Matrix &Omega = *m_Omega;

    // Get needed matrices
    vnl_matrix<Real> tmp;
    Matrix &Y = *m_Y;
    Matrix &L = *m_L;
    Matrix Yt;
    Matrix YtY;
    Matrix YtYplL;
    Matrix YtYplLInverted;

    // Get Y transpose
    tmp = Y.GetTranspose();
    Yt  = tmp;

    // Get YtY
    YtY = Yt * Y;

    // Get YtY + lL
    L *= m_lambda;
    YtYplL = YtY + L;

    // Get (YtYplL)^-1
    tmp         = YtYplL.GetInverse();
    YtYplLInverted = tmp;

    #ifndef NDEBUG
        std::cerr << "Yt (" << Yt.Rows() << "x" << Yt.Cols() << ")" << std::endl;
        Yt.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;

        std::cerr << "YtY (" << YtY.Rows() << "x" << YtY.Cols() << ")" << std::endl;
        YtY.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;

        std::cerr << "YtYplL (" << YtYplL.Rows() << "x" << YtYplL.Cols() << ")" << std::endl;
        YtYplL.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;

        std::cerr << "YtYplLInverted (" << YtYplLInverted.Rows() << "x" << YtYplLInverted.Cols() << ")" << std::endl;
        YtYplLInverted.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG

    // Compute Omega
    // Omega = (Yt * Y + l*L)^-1 * Yt
    Omega = YtYplLInverted * Yt;

    #ifndef NDEBUG
        std::cerr << "Omega (" << Omega.Rows() << "x" << Omega.Cols() << ")" << std::endl;
        Omega.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

} // namespace btk

