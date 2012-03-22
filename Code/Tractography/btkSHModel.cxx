/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

12 februar 2010
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


#include "btkSHModel.h"


// STL includes
#include "cmath"

// Local includes
#include "btkSphericalHarmonics.h"


#define SPH_RESOLUTION 120
//#define MAXIMA_THRESHOLD 0.45


namespace btk
{

SHModel::SHModel(const std::string &filename, const std::string &directionsfilename, char displayMode)
{
    m_model      = 0;
    m_interp     = 0;
    m_directions = 0;

    m_displayMode = displayMode;

    m_Y     = 0;
    m_Yori  = 0;
    m_P     = 0;
    m_Sharp = 0;

    m_reso = SPH_RESOLUTION;

    m_maxTheta = M_PI;
    m_pasTheta = m_maxTheta / (Real)m_reso;

    m_maxPhi = 2.0 * M_PI;
//    m_maxPhi = M_PI;
    m_pasPhi = m_maxPhi / (Real)m_reso;

    m_4PI = 4. * M_PI;
    m_2PI = 2. * M_PI;


    // Reading files
    Sequence::Pointer model = this->readFiles(filename, directionsfilename);

    // Getting images sizes
    SequenceRegion smodelRegion = model->GetLargestPossibleRegion();

    unsigned int xMax = smodelRegion.GetSize(0);
    unsigned int yMax = smodelRegion.GetSize(1);
    unsigned int zMax = smodelRegion.GetSize(2);
    unsigned int kMax = smodelRegion.GetSize(3);

//    std::cout << "\tThere are " << kMax << " images of size ";
//    std::cout << xMax << "x" << yMax << "x" << zMax << "." << std::endl;


    // Compute order
    Real discr = 2.25 - 2.*(1. - kMax);
    m_order    = -1.5 + std::sqrt(discr);
    m_R        = kMax;


//    std::cout << "\tModel is at order " << m_order << "." << std::endl;


    // Generate directions following resolution
//    std::cout << "\tGenerating directions (sampling resolution " << m_reso << "x" << m_reso << ")..." << std::flush;
    m_directions = new std::vector<Direction>;

    for(unsigned int i=0; i<m_reso+1; i++)
    {
        for(unsigned int j=0; j<m_reso; j++)
            m_directions->push_back(Direction(i*m_pasTheta,j*m_pasPhi));
    }
//    for(unsigned int i=0; i<m_reso; i++)
//    {
//        for(unsigned int j=0; j<m_reso-1; j++)
//            m_directions->push_back(Direction(i*m_pasTheta,j*m_pasPhi));
//    }
//    std::cout << "done." << std::endl;


    // SH basis matrix
//    std::cout << "\tSpherical harmonics basis matrix..." << std::flush;
    this->computeSHBasisMatrix();
    this->computeSHBasisOriMatrix();
//    std::cout << "done." << std::endl;

    // Compute Legendre matrix
//    std::cout << "\tComputing Legendre matrix..." << std::flush;
    this->computeLegendreMatrix();
//    std::cout << "done." << std::endl;

    // Build sharp coefficients
//    std::cout << "\tBuilding sharper ODF matrix..." << std::flush;
    this->buildSharperODFMatrix();
//    std::cout << "done." << std::endl;


    // Allocate space memory for the array of images
//    std::cout << "\tAllocating space memory..." << std::flush;
    m_model  = new Image::Pointer[kMax];
    m_interp = new ImageInterpolator::Pointer[kMax];
//    std::cout << "done." << std::endl;


//    std::cout << "\tPreparing and interpolating data..." << std::flush;

    // Define images region
    ImageRegion iRegion;

    Image::IndexType iIndex;
    iIndex[0] = 0; iIndex[1] = 0; iIndex[2] = 0;
    iRegion.SetIndex(iIndex);

    Image::SizeType iSize;
    iSize[0] = xMax; iSize[1] = yMax; iSize[2] = zMax;
    iRegion.SetSize(iSize);


    // Define sequence region
    SequenceRegion sRegion;

    Sequence::SizeType sSize;
    sSize[0] = xMax; sSize[1] = yMax;
    sSize[2] = zMax; sSize[3] = 1;
    sRegion.SetSize(sSize);


    // Spacing and origin
    Sequence::SpacingType spacing = model->GetSpacing();
    Sequence::PointType origin    = model->GetOrigin();

    Image::SpacingType ispacing;
    ispacing[0] = spacing[0]; ispacing[1] = spacing[1]; ispacing[2] = spacing[2];
    Image::PointType iorigin;
    iorigin[0] = origin[0]; iorigin[1] = origin[1]; iorigin[2] = origin[2];

    itk::Matrix<Real,3,3> dirMat;
    itk::Matrix<Real,4,4> iniMat = model->GetDirection();

    dirMat(0,0) = iniMat(0,0);
    dirMat(0,1) = iniMat(0,1);
    dirMat(0,2) = iniMat(0,2);

    dirMat(1,0) = iniMat(1,0);
    dirMat(1,1) = iniMat(1,1);
    dirMat(1,2) = iniMat(1,2);

    dirMat(2,0) = iniMat(2,0);
    dirMat(2,1) = iniMat(2,1);
    dirMat(2,2) = iniMat(2,2);


    // Put sequence in array of images
    for(unsigned int k=0; k<kMax; k++)
    {
        // Allocate space memory
        m_model[k] = Image::New();
        m_model[k]->SetRegions(iSize);
        m_model[k]->SetSpacing(ispacing);
        m_model[k]->SetOrigin(iorigin);
        m_model[k]->SetDirection(dirMat);
        m_model[k]->Allocate();

        // Define correct sequence region
        Sequence::IndexType sIndex;
        sIndex[0] = 0; sIndex[1] = 0;
        sIndex[2] = 0; sIndex[3] = k;
        sRegion.SetIndex(sIndex);

        // Define iterators
        SequenceIterator sIt(model, sRegion);
        ImageIterator iIt(m_model[k], iRegion);

        // Copy data
        for(sIt.GoToBegin(), iIt.GoToBegin(); !sIt.IsAtEnd() && !iIt.IsAtEnd(); ++sIt, ++iIt)
            iIt.Set(sIt.Get());

        // Define interpolator
        m_interp[k] = ImageInterpolator::New();
        m_interp[k]->SetInputImage(m_model[k]);
    } // for k

//    std::cout << "done." << std::endl;


//    std::cout << "done." << std::endl;
}

SHModel::SHModel(Sequence::Pointer model, std::vector<Direction> *originalDirections, char displayMode)
{
    m_model      = 0;
    m_interp     = 0;
    m_directions = 0;

    m_Y     = 0;
    m_Yori  = 0;
    m_P     = 0;
    m_Sharp = 0;

    m_displayMode = displayMode;

    m_originalDirections = originalDirections;


    Display1(m_displayMode, std::cout << "Loading model..." << std::endl);

    m_reso = SPH_RESOLUTION;

    m_maxTheta = M_PI;
    m_pasTheta = m_maxTheta / (Real)m_reso;

    m_maxPhi = 2.0 * M_PI;
//    m_maxPhi = M_PI;
    m_pasPhi = m_maxPhi / (Real)m_reso;

    m_4PI = 4. * M_PI;
    m_2PI = 2. * M_PI;

    // Getting images sizes
    SequenceRegion smodelRegion = model->GetLargestPossibleRegion();

    unsigned int xMax = smodelRegion.GetSize(0);
    unsigned int yMax = smodelRegion.GetSize(1);
    unsigned int zMax = smodelRegion.GetSize(2);
    unsigned int kMax = smodelRegion.GetSize(3);

    Display2(m_displayMode, std::cout << "\tThere are " << kMax << " images of size ");
    Display2(m_displayMode, std::cout << xMax << "x" << yMax << "x" << zMax << "." << std::endl);


    // Compute order
    Real discr = 2.25 - 2.*(1. - kMax);
    m_order    = -1.5 + std::sqrt(discr);
    m_R        = kMax;


    Display2(m_displayMode, std::cout << "\tModel is at order " << m_order << "." << std::endl);


    // Generate directions following resolution
    Display2(m_displayMode, std::cout << "\tGenerating directions (sampling resolution " << m_reso << "x" << m_reso << ")..." << std::flush);
    m_directions = new std::vector<Direction>;

    for(unsigned int i=0; i<m_reso+1; i++)
    {
        for(unsigned int j=0; j<m_reso; j++)
            m_directions->push_back(Direction(i*m_pasTheta,j*m_pasPhi));
    }
//    for(unsigned int i=0; i<m_reso; i++)
//    {
//        for(unsigned int j=0; j<m_reso-1; j++)
//            m_directions->push_back(Direction(i*m_pasTheta,j*m_pasPhi));
//    }
    Display2(m_displayMode, std::cout << "done." << std::endl);


    // SH basis matrix
    Display2(m_displayMode, std::cout << "\tSpherical harmonics basis matrix..." << std::flush);
    this->computeSHBasisMatrix();
    this->computeSHBasisOriMatrix();
    Display2(m_displayMode, std::cout << "done." << std::endl);

    // Compute Legendre matrix
    Display2(m_displayMode, std::cout << "\tComputing Legendre matrix..." << std::flush);
    this->computeLegendreMatrix();
    Display2(m_displayMode, std::cout << "done." << std::endl);

    // Build sharp coefficients
    Display2(m_displayMode, std::cout << "\tBuilding sharper ODF matrix..." << std::flush);
    this->buildSharperODFMatrix();
    Display2(m_displayMode, std::cout << "done." << std::endl);


    // Allocate space memory for the array of images
    Display2(m_displayMode, std::cout << "\tAllocating space memory..." << std::flush);
    m_model  = new Image::Pointer[kMax];
    m_interp = new ImageInterpolator::Pointer[kMax];
    Display2(m_displayMode, std::cout << "done." << std::endl);


    Display2(m_displayMode, std::cout << "\tPreparing and interpolating data..." << std::flush);

    // Define images region
    ImageRegion iRegion;

    Image::IndexType iIndex;
    iIndex[0] = 0; iIndex[1] = 0; iIndex[2] = 0;
    iRegion.SetIndex(iIndex);

    Image::SizeType iSize;
    iSize[0] = xMax; iSize[1] = yMax; iSize[2] = zMax;
    iRegion.SetSize(iSize);


    // Define sequence region
    SequenceRegion sRegion;

    Sequence::SizeType sSize;
    sSize[0] = xMax; sSize[1] = yMax;
    sSize[2] = zMax; sSize[3] = 1;
    sRegion.SetSize(sSize);


    // Spacing and origin
    Sequence::SpacingType spacing = model->GetSpacing();
    Sequence::PointType origin    = model->GetOrigin();

    Image::SpacingType ispacing;
    ispacing[0] = spacing[0]; ispacing[1] = spacing[1]; ispacing[2] = spacing[2];
    Image::PointType iorigin;
    iorigin[0] = origin[0]; iorigin[1] = origin[1]; iorigin[2] = origin[2];

    itk::Matrix<Real,3,3> dirMat;
    itk::Matrix<Real,4,4> iniMat = model->GetDirection();

    dirMat(0,0) = iniMat(0,0);
    dirMat(0,1) = iniMat(0,1);
    dirMat(0,2) = iniMat(0,2);

    dirMat(1,0) = iniMat(1,0);
    dirMat(1,1) = iniMat(1,1);
    dirMat(1,2) = iniMat(1,2);

    dirMat(2,0) = iniMat(2,0);
    dirMat(2,1) = iniMat(2,1);
    dirMat(2,2) = iniMat(2,2);


    // Put sequence in array of images
    for(unsigned int k=0; k<kMax; k++)
    {
        // Allocate space memory
        m_model[k] = Image::New();
        m_model[k]->SetRegions(iSize);
        m_model[k]->SetSpacing(ispacing);
        m_model[k]->SetOrigin(iorigin);
        m_model[k]->SetDirection(dirMat);
        m_model[k]->Allocate();

        // Define correct sequence region
        Sequence::IndexType sIndex;
        sIndex[0] = 0; sIndex[1] = 0;
        sIndex[2] = 0; sIndex[3] = k;
        sRegion.SetIndex(sIndex);

        // Define iterators
        SequenceIterator sIt(model, sRegion);
        ImageIterator iIt(m_model[k], iRegion);

        // Copy data
        for(sIt.GoToBegin(), iIt.GoToBegin(); !sIt.IsAtEnd() && !iIt.IsAtEnd(); ++sIt, ++iIt)
            iIt.Set(sIt.Get());

        // Define interpolator
        m_interp[k] = ImageInterpolator::New();
        m_interp[k]->SetInputImage(m_model[k]);
    } // for k

    Display2(m_displayMode, std::cout << "done." << std::endl);


    Display1(m_displayMode, std::cout << "done." << std::endl);
}

Sequence::Pointer SHModel::readFiles(const std::string &filename, const std::string &directionsfilename)
{
//    std::cout << "Loading model file \"" << filename << "\"..." << std::endl;

    //
    // Read model file
    //

//        std::cout << "\tReading model's file..." << std::flush;
        SequenceReader::Pointer reader = SequenceReader::New();
        reader->SetFileName(filename);
        reader->Update();
//        std::cout << "done." << std::endl;

    //
    // Read direction's file
    //

//        std::cout << "\tReading gradient directions' file..." << std::flush;

        // Open file
        std::fstream dirFile(directionsfilename.c_str(), std::fstream::in);

        if(!dirFile.is_open())  // Unable to open file
        {
            std::cout << "Error: unable to read file !" << std::endl;
            exit(EXIT_FAILURE);
        }
        else                           // File is open
        {
            m_originalDirections = new std::vector<Direction>;

            Real theta = 0, phi = 0;

            while((dirFile >> theta) && (dirFile >> phi))
                m_originalDirections->push_back(Direction(theta,phi));
        }

        // Close file
        dirFile.close();

//        std::cout << "done." << std::endl;



    return reader->GetOutput();
}

SHModel::~SHModel()
{
    delete[] m_interp;
    delete[] m_model;

    delete m_directions;

    delete m_Y;
    delete m_Yori;
    delete m_P;
    delete m_Sharp;
}

Real SHModel::signalAt(Direction u, Point p)
{
    // Compute spherical harmonics basis matrix
    Matrix Y(1, m_R);

    unsigned int j = 0;

    for(unsigned int l=0; l<=m_order; l+=2)
    {
        Real coeff = std::sqrt((2.*l + 1.) / m_4PI);

        for(int m=-(int)l; m<=(int)l; m++)
            Y(0,j++) = coeff * SphericalHarmonics::computeBasis(u,l,m);
    } // for l


    // Evaluate signal at u
    Matrix C(m_R, 1);

    ImageInterpolator::ContinuousIndexType index;
    index[0] = p.x(); index[1] = p.y(); index[2] = p.z();

    for(unsigned int i=0; i<m_R; i++)
        C(i,0) = m_interp[i]->EvaluateAtContinuousIndex(index);


    return (Y * C)(0,0);
}

Matrix SHModel::signalAt(Point p)
{
    assert(m_P);
    assert(m_Yori);


    Matrix &Y     = *m_Yori;


    // Evaluate Orientation Diffusion Function at each gradient direction
    Matrix C(m_R, 1);

    ImageInterpolator::ContinuousIndexType index;
    index[0] = p.x(); index[1] = p.y(); index[2] = p.z();

    for(unsigned int i=0; i<m_R; i++)
        C(i,0) = m_interp[i]->EvaluateAtContinuousIndex(index);


    return (Y * C);
}

Matrix SHModel::signalAt(Point p, std::vector<Direction> *g)
{
    // Compute spherical harmonics basis matrix
    Matrix Y(g->size(), m_R);

    // Compute Y matrix
    Real _4PI  = 4.0 * M_PI;

    unsigned int j = 0;

    for(unsigned int u=0; u<g->size(); u++)
    {
        for(unsigned int l=0; l<=m_order; l+=2)
        {
            Real coef = std::sqrt((2.0*(Real)l + 1.0) / _4PI);

            for(int m=-(int)l; m<=(int)l; m++)
                Y(u,j++) = coef * SphericalHarmonics::computeBasis(g->at(u), l, m);
        } // for l

        j = 0;
    } // for u


    // Evaluate Orientation Diffusion Function at each gradient direction
    Matrix C(m_R, 1);

    ImageInterpolator::ContinuousIndexType index;
    index[0] = p.x(); index[1] = p.y(); index[2] = p.z();

    for(unsigned int i=0; i<m_R; i++)
        C(i,0) = m_interp[i]->EvaluateAtContinuousIndex(index);


    return (Y * C);
}

Real SHModel::odfAt(Direction u, Point p)
{
    assert(m_P);
    assert(m_Sharp);


    Matrix &P     = *m_P;
    Matrix &Sharp = *m_Sharp;


    // Compute spherical harmonics basis matrix
    Matrix Y(1, m_R);

    unsigned int j = 0;

    for(unsigned int l=0; l<=m_order; l+=2)
    {
        Real coeff = std::sqrt((2.*l + 1.) / m_4PI);

        for(int m=-(int)l; m<=(int)l; m++)
            Y(0,j++) = coeff * SphericalHarmonics::computeBasis(u,l,m);
    } // for l


    // Evaluate Orientation Diffusion Function at u
    Matrix C(m_R, 1);

    ImageInterpolator::ContinuousIndexType index;
    index[0] = p.x(); index[1] = p.y(); index[2] = p.z();

    for(unsigned int i=0; i<m_R; i++)
        C(i,0) = m_interp[i]->EvaluateAtContinuousIndex(index);


    return (Y * (P * (Sharp * C)))(0,0);
//    return (Y * (P * C))(0,0);
}

Matrix SHModel::odfAt(Point p)
{
    assert(m_P);
    assert(m_Y);


    Matrix &Y     = *m_Y;
    Matrix &P     = *m_P;
    Matrix &Sharp = *m_Sharp;


    // Evaluate Orientation Diffusion Function at each gradient direction
    Matrix C(m_R, 1);

    ImageInterpolator::ContinuousIndexType index;
    index[0] = p.x(); index[1] = p.y(); index[2] = p.z();

    for(unsigned int i=0; i<m_R; i++)
        C(i,0) = m_interp[i]->EvaluateAtContinuousIndex(index);


    return (Y * (P * (Sharp * C)));
//    return (Y * (P * C));
}

Direction SHModel::getMaxDirectionAt(Point p)
{
    Matrix Psi       = this->odfAt(p);
    Real max         = Psi(0,0);
    Direction maxDir = m_directions->at(0);

    for(unsigned int i=1; i<Psi.Rows(); i++)
    {
        if(Psi(i,0) > max)
        {
            max    = Psi(i,0);
            maxDir = m_directions->at(i);
        } // Psi(i,0) <= max
    } // for each sampled directions

    return maxDir;
}

std::vector<Direction> SHModel::getMaxDirectionsAt(Point p)
{
    std::vector<Direction> maxima;

    // Get data (max and min values at this position)
    Matrix Psi = this->odfAt(p);

    Real min = Psi(0,0), max = Psi(0,0);
    for(unsigned int i=1; i<Psi.Rows(); i++)
    {
        if(min > Psi(i,0))
            min = Psi(i,0);

        if(max < Psi(i,0))
            max = Psi(i,0);
    }

    unsigned int thetaRes = m_reso+1;
    unsigned int phiRes   = m_reso;


    for(unsigned int i=1; i<thetaRes-1; i++)
    {
        for(unsigned int j=0; j<phiRes; j++)
        {
            Real odf1, odf2, odf3, odf4, odf5, odf6, odf7, odf8;
            Real odf = Psi(phiRes*i + j, 0);

            if(j == 0)
            {
                odf1 = Psi(phiRes*(i-1)+(phiRes-1), 0);
                odf2 = Psi(phiRes*(i-1)+j, 0);
                odf3 = Psi(phiRes*(i-1)+(j+1), 0);
                odf4 = Psi(phiRes*i+(phiRes-1), 0);
                odf5 = Psi(phiRes*i+(j+1), 0);
                odf6 = Psi(phiRes*(i+1)+(phiRes-1), 0);
                odf7 = Psi(phiRes*(i+1)+j, 0);
                odf8 = Psi(phiRes*(i+1)+(j+1), 0);
            }
            else if(j == phiRes-1)
            {
                odf1 = Psi(phiRes*(i-1)+(j-1), 0);
                odf2 = Psi(phiRes*(i-1)+j, 0);
                odf3 = Psi(phiRes*(i-1)+(0), 0);
                odf4 = Psi(phiRes*i+(j-1), 0);
                odf5 = Psi(phiRes*i+(0), 0);
                odf6 = Psi(phiRes*(i+1)+(j-1), 0);
                odf7 = Psi(phiRes*(i+1)+j, 0);
                odf8 = Psi(phiRes*(i+1)+(0), 0);
            }
            else // j != 0 && j != phiRes-1
            {
                odf1 = Psi(phiRes*(i-1)+(j-1), 0);
                odf2 = Psi(phiRes*(i-1)+j, 0);
                odf3 = Psi(phiRes*(i-1)+(j+1), 0);
                odf4 = Psi(phiRes*i+(j-1), 0);
                odf5 = Psi(phiRes*i+(j+1), 0);
                odf6 = Psi(phiRes*(i+1)+(j-1), 0);
                odf7 = Psi(phiRes*(i+1)+j, 0);
                odf8 = Psi(phiRes*(i+1)+(j+1), 0);
            }

            if(odf > odf1 && odf > odf2 && odf > odf3 && odf > odf4 && odf > odf5 && odf > odf6 && odf > odf7 && odf > odf8)
            {
                if((Psi(phiRes*i+j,0)-min)/(max-min) > 0.9)
                    maxima.push_back(Direction(i*m_pasTheta,j*m_pasPhi));
            }
        } // for each phi
    } // for each theta

//    unsigned int phiRes = m_reso-1;
//    bool isGreater;
//    unsigned int j;
//
//    // i=0
//    isGreater = true;
//    j = 0;
//
//    do
//    {
//        isGreater = ( Psi(phiRes + j,0) < Psi(0,0) );
//        j++;
//    } while(isGreater && j<phiRes);
//
//    if(isGreater)
//    {
//        if((Psi(0,0)-min)/(max-min) > MAXIMA_THRESHOLD)
//        {
//            maxima.push_back(Direction(0,0));
//            maxima.push_back(Direction(M_PI,0));
//        }
//    }
//
////    // i=|theta|
////    isGreater = true;
////    j = 0;
////
////    do
////    {
////        isGreater = ( Psi(phiRes*phiRes,0) < Psi(m_reso-1,0) );
////        j++;
////    } while(isGreater && j<phiRes);
////
////    if(isGreater)
////    {
////        if((Psi(m_reso-1,0)-min)/(max-min) > MAXIMA_THRESHOLD)
////        {
////            maxima.push_back(Direction(M_PI,0));
////            maxima.push_back(Direction(0,0));
////        }
////    }
//
//    // j=0 & general case
//    for(unsigned int i=1; i<m_reso-1; i++)
//    {
//        for(unsigned int j=0; j<phiRes; j++)
//        {
//            Real odf, odf1, odf2, odf3, odf4, odf5, odf6, odf7, odf8;
//
//            if(j == 0)
//            {
//                odf = Psi(phiRes*i,0);
//
//                odf1 = Psi(phiRes*(i-1),0);
//                odf2 = Psi(phiRes*(i+1),0);
//                odf3 = Psi(phiRes*(i-1) + (phiRes-1),0);
//                odf4 = Psi(phiRes*i + (phiRes-1),0);
//                odf5 = Psi(phiRes*(i+1) + (phiRes-1),0);
//                odf6 = Psi(phiRes*(i-1) + 1,0);
//                odf7 = Psi(phiRes*i + 1,0);
//                odf8 = Psi(phiRes*(i+1) + 1,0);
//            }
//            else if(j == phiRes-1)
//            {
//                odf = Psi(phiRes*i,phiRes-1);
//
//                odf1 = Psi(phiRes*(i-1) + (phiRes-1),0);
//                odf2 = Psi(phiRes*(i+1) + (phiRes-1),0);
//                odf3 = Psi(phiRes*(i-1) + (phiRes-2),0);
//                odf4 = Psi(phiRes*i + (phiRes-2),0);
//                odf5 = Psi(phiRes*(i+1) + (phiRes-2),0);
//                odf6 = Psi(phiRes*(i-1),0);
//                odf7 = Psi(phiRes*i,0);
//                odf8 = Psi(phiRes*(i+1),0);
//            }
//            else // general case
//            {
//                odf = Psi(phiRes*i + j,0);
//
//                odf1 = Psi(phiRes*(i-1) + (j+1),0);
//                odf2 = Psi(phiRes*i + (j+1),0);
//                odf3 = Psi(phiRes*(i+1) + (j+1),0);
//                odf4 = Psi(phiRes*(i-1) + j,0);
//                odf5 = Psi(phiRes*(i+1) + j,0);
//                odf6 = Psi(phiRes*(i-1) + (j-1),0);
//                odf7 = Psi(phiRes*i + (j-1),0);
//                odf8 = Psi(phiRes*(i+1) + (j-1),0);
//            }
//
//            if(odf > odf1 && odf > odf2 && odf > odf3 && odf > odf4 && odf > odf5 && odf > odf6 && odf > odf7 && odf > odf8)
//            {
//                if((Psi(phiRes*i,0)-min)/(max-min) > MAXIMA_THRESHOLD)
//                {
//                    Direction u(i*m_pasTheta,j*m_pasPhi);
//                    Vector v = u.toVector();
//
//                    maxima.push_back(u);
//                    maxima.push_back(Vector(-v.x(), -v.y(), -v.z()).toDirection());
//                }
//            }
//        }
//    }


    return maxima;
}

void SHModel::computeLegendreMatrix()
{
    assert(!m_P);

    m_P = new Matrix(m_R, m_R);

    Matrix &P = *m_P;

    P.Fill(0);


    unsigned int j = 0;

    for(unsigned int l=0; l<=m_order; l+=2)
    {
        for(int m=-(int)l; m<=(int)l; m++)
        {
            P(j,j) = m_2PI * SphericalHarmonics::legendrePolynomialInZero(l);
            j++;
        } // for m
    } // for l

    #ifndef NDEBUG
        std::cerr << "P (" << P.Rows() << "x" << P.Cols() << ")" << std::endl;
        P.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

void SHModel::computeSHBasisMatrix()
{
    assert(m_directions);
    assert(!m_Y);


    // Allocate memory space for Y matrix
    m_Y = new Matrix(m_directions->size(), m_R);

    Matrix &Y = *m_Y;

    // Compute Y matrix
    Real _4PI  = 4.0 * M_PI;

    unsigned int j = 0;

    for(unsigned int u=0; u<m_directions->size(); u++)
    {
        for(unsigned int l=0; l<=m_order; l+=2)
        {
            Real coef = std::sqrt((2.0*(Real)l + 1.0) / _4PI);

            for(int m=-(int)l; m<=(int)l; m++)
                Y(u,j++) = coef * SphericalHarmonics::computeBasis(m_directions->at(u), l, m);
        } // for l

        j = 0;
    } // for u

    #ifndef NDEBUG
        std::cerr << "Y (" << Y.Rows() << "x" << Y.Cols() << ")" << std::endl;
        Y.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

void SHModel::computeSHBasisOriMatrix()
{
    assert(m_originalDirections);
    assert(!m_Yori);


    // Allocate memory space for Y matrix
    m_Yori = new Matrix(m_originalDirections->size(), m_R);

    Matrix &Y = *m_Yori;

    // Compute Y matrix
    Real _4PI  = 4.0 * M_PI;

    unsigned int j = 0;

    for(unsigned int u=0; u<m_originalDirections->size(); u++)
    {
        for(unsigned int l=0; l<=m_order; l+=2)
        {
            Real coef = std::sqrt((2.0*(Real)l + 1.0) / _4PI);

            for(int m=-(int)l; m<=(int)l; m++)
                Y(u,j++) = coef * SphericalHarmonics::computeBasis(m_originalDirections->at(u), l, m);
        } // for l

        j = 0;
    } // for u

    #ifndef NDEBUG
        std::cerr << "Y (" << Y.Rows() << "x" << Y.Cols() << ")" << std::endl;
        Y.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

void SHModel::buildSharperODFMatrix()
{
    assert(!m_Sharp);


    m_Sharp = new Matrix(m_R,m_R);

    Matrix &Sharp = *m_Sharp;

    Sharp.Fill(0);


    // Parameters
    Real l1 = 1.0;
    Real l2 = 0.1;
    Real a  = 0.9; // 1 - l2 / l1
    Real b  = 1500.;

    // Pre-computed values
    Real sqrta = std::sqrt(a);
    Real pre1 = std::asin(sqrta);
    Real pre2 = std::sqrt(0.1) * sqrta; // sqrt(1-a) * sqrt(a)
    Real coeff = 4. * b * std::sqrt(l1*l2);

    // Coefficients
    Real A0, A2, A4, A6, A8;

    if(m_order >=0)
        A0 = 1./(std::pow(a,0.5)) * 2. * pre1;

    if(m_order >= 2)
        A2 = -1./(2.*std::pow(a,1.5)) * ( (-3.+2.*a)*pre1 + 3.*pre2 );

    if(m_order >= 4)
        A4 = 1./(32.*std::pow(a,2.5)) * ( (105.-120.*a+24.*std::pow(a,2.))*pre1 + (-105.+50.*a)*pre2 );

    if(m_order >= 6)
        A6 = -1./(128.*std::pow(a,3.5)) * ( (-1155.+1890.*a-840.*std::pow(a,2.)+80.*std::pow(a,3.))*pre1 + (1155.-1120.*a+196.*std::pow(a,2.))*pre2 );

    if(m_order >= 8)
        A8 = 1./(8192.* std::pow(a,4.5)) *
            ( (225225.-480480.*a+332640.* std::pow(a,2.)-80640.* std::pow(a,3.)+4480.* std::pow(a,4.))*pre1 +
            (-225225.+330330.*a-132440.* std::pow(a,2.)+12176.* std::pow(a,3.))*pre2 );


    // Compute sharper coefficients in a diagonal matrix
    if(m_R >= 1)
        Sharp(0,0) = coeff/A0;

    if(m_R >= 6)
    {
        for(unsigned int i=1; i<6; i++)
            Sharp(i,i) = coeff/A2;
    }

    if(m_R >= 15)
    {
        for(unsigned int i=6; i<15; i++)
            Sharp(i,i) = coeff/A4;
    }

    if(m_R >= 28)
    {
        for(unsigned int i=15; i<28; i++)
            Sharp(i,i) = coeff/A6;
    }

    if(m_R >= 45)
    {
        for(unsigned int i=28; i<45; i++)
            Sharp(i,i) = coeff/A8;
    }

    #ifndef NDEBUG
        std::cerr << "Sharp (" << Sharp.Rows() << "x" << Sharp.Cols() << ")" << std::endl;
        Sharp.GetVnlMatrix().print(std::cerr);
        std::cerr << std::endl;
    #endif // NDEBUG
}

} // namespace btk

