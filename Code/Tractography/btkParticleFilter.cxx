/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

31 march 2010
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

#include "btkParticleFilter.h"


// STL includes
#include "iostream"
#include "fstream"
#include "cstdlib"
#include "ctime"
#include "sstream"
#include "cmath"
#include "algorithm"
#include "utility"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkColorTransferFunction.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkAppendPolyData.h"

// ITK includes
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPoint.h"
#include "itkContinuousIndex.h"


#define Ind(i,j) (m_M*(i) + (j))

#define KAPPA 180


namespace btk
{

ParticleFilter::ParticleFilter(Signal *signal, SHModel *model, APrioriDensity aPriori, LikelihoodDensity likelihood, ImportanceDensity importance,
                               const std::string maskFileName, Image::SizeType size, Image::PointType origin, Image::SpacingType spacing,
                               unsigned int M, Point x0, Real epsilon, Real stepSize, unsigned int maxLength) :
        m_aPriori(aPriori), m_likelihood(likelihood), m_importance(importance)
{
    std::cout << "\tInitializing filter..." << std::endl;

    m_k         = 0;
    m_epsilon   = epsilon;
    m_M         = M;
    m_stepSize  = stepSize;
    m_maxLength = maxLength;
    m_origin    = origin;
    m_spacing   = spacing;
    m_model     = model;
    m_signal    = signal;
    m_dirNum    = 1;
    m_kx        = m_stepSize/m_spacing[0];
    m_ky        = m_stepSize/m_spacing[1];
    m_kz        = m_stepSize/m_spacing[2];
    m_lps       = false;

    m_vStepSize.push_back(m_kx);
    m_vStepSize.push_back(m_ky);
    m_vStepSize.push_back(m_kz);

    std::cout << "\t\tFilter will use " << M << " particles." << std::endl;
    std::cout << "\t\tResampling treshold is set as " << epsilon << "." << std::endl;
    std::cout << "\t\tThe moving step is fixed at " << m_stepSize << " mm." << std::endl;
    std::cout << "\t\tThere will be " << m_maxLength << " steps." << std::endl;


    // Initialize random number generator
    std::srand(std::time(NULL));


    // Create density map
    m_map = Image::New();
    m_map->SetRegions(size);
    m_map->SetOrigin(origin);
    m_map->SetSpacing(spacing);
    m_map->SetDirection(m_model->GetDirection());
    m_map->Allocate();
    m_map->FillBuffer(0);

    itk::Point<Real,3> worldPoint;
    worldPoint[0] = x0.x(); worldPoint[1] = x0.y(); worldPoint[2] = x0.z();

    itk::ContinuousIndex<Real,3> continuousIndex;
    if(!m_map->TransformPhysicalPointToContinuousIndex(worldPoint, continuousIndex))
    {
        std::cout << "Error: continuous index not in image !" << std::endl;
        std::cerr << "(" << worldPoint[0] << "," << worldPoint[1] << "," << worldPoint[2] << ") --> (" << continuousIndex[0] << "," << continuousIndex[1] << "," << continuousIndex[2] << ")" << std::endl;
        exit(EXIT_FAILURE);
    }

    m_x0 = Point(continuousIndex[0], continuousIndex[1], continuousIndex[2]);


    // Read mask image
    MaskReader::Pointer reader = MaskReader::New();
    reader->SetFileName(maskFileName);
    reader->Update();
    m_mask = reader->GetOutput();


    std::cout << "\tdone." << std::endl;
}

ParticleFilter::ParticleFilter(Signal *signal, SHModel *model, APrioriDensity aPriori, LikelihoodDensity likelihood, ImportanceDensity importance,
                               Mask::Pointer mask, Image::SizeType size, Image::PointType origin, Image::SpacingType spacing,
                               unsigned int M, Point x0, Real epsilon, Real stepSize, char displaMode) :
        m_aPriori(aPriori), m_likelihood(likelihood), m_importance(importance)
{
    m_displayMode = displaMode;


    Display2(m_displayMode, std::cout << "\tInitializing filter..." << std::endl);

    m_k         = 0;
    m_epsilon   = epsilon;
    m_M         = M;
    m_stepSize  = stepSize;
    m_origin    = origin;
    m_spacing   = spacing;
    m_model     = model;
    m_signal    = signal;
    m_dirNum    = 1;
    m_kx        = m_stepSize/m_spacing[0];
    m_ky        = m_stepSize/m_spacing[1];
    m_kz        = m_stepSize/m_spacing[2];
    m_lps       = false;

    m_vStepSize.push_back(m_kx);
    m_vStepSize.push_back(m_ky);
    m_vStepSize.push_back(m_kz);


    Display2(m_displayMode, std::cout << "\t\tFilter will use " << M << " particles." << std::endl);
    Display2(m_displayMode, std::cout << "\t\tResampling treshold is set as " << epsilon << "." << std::endl);
    Display2(m_displayMode, std::cout << "\t\tThe moving step is fixed at " << m_stepSize << " mm." << std::endl);


    // Initialize random number generator
    std::srand(std::time(NULL));


    // Create density map
    m_map = Image::New();
    m_map->SetRegions(size);
    m_map->SetOrigin(origin);
    m_map->SetSpacing(spacing);
    m_map->SetDirection(m_model->GetDirection());
    m_map->Allocate();
    m_map->FillBuffer(0);

    itk::Point<Real,3> worldPoint;
    worldPoint[0] = x0.x(); worldPoint[1] = x0.y(); worldPoint[2] = x0.z();

    itk::ContinuousIndex<Real,3> continuousIndex;
    if(!m_map->TransformPhysicalPointToContinuousIndex(worldPoint, continuousIndex))
    {
        std::cout << "Error: continuous index not in image !" << std::endl;
        std::cerr << "(" << worldPoint[0] << "," << worldPoint[1] << "," << worldPoint[2] << ") --> (" << continuousIndex[0] << "," << continuousIndex[1] << "," << continuousIndex[2] << ")" << std::endl;
        exit(EXIT_FAILURE);
    }

    m_x0 = Point(continuousIndex[0], continuousIndex[1], continuousIndex[2]);

    m_mask = mask;


    Display2(m_displayMode, std::cout << "\tdone." << std::endl);
}

void ParticleFilter::run(int label)
{
    //
    // Initialization
    //

    // Get max direction à starting point
    std::vector<Direction> initDir = m_model->getMaxDirectionsAt(m_x0);

    // Find maximal direction
    Direction maxDir;
    Real maxValue = 0;

    for(std::vector<Direction>::iterator it=initDir.begin(); it != initDir.end(); it++)
    {
        Real tmp = m_model->odfAt(*it, m_x0);

        if(maxValue < tmp)
        {
            maxValue = tmp;
            maxDir   = *it;
        }
    }

    // Compute Symetric direction
    Direction symDir(
        M_PI - maxDir.theta(),
        (maxDir.phi() < M_PI) ? maxDir.phi() + M_PI : maxDir.phi() - M_PI
    );


    //
    // Filtering
    //

    this->run(label, maxDir);
    this->ComputeMap();
    Particle map1 = this->GetMAP();

    this->run(label, symDir);
    this->ComputeMap();
    Particle map2 = this->GetMAP();

    if(map1.length() + map2.length() > 2)
        this->ComputeFiber(map1,map2);
}

void ParticleFilter::run(int label, Direction dir)
{
    m_k     = 0;
    m_cloud = std::vector<Particle>(m_M, Particle(m_x0));


    //
    // Initial sampling (vMF in mean direction dir and a priori kappa)
    //

    Display2(m_displayMode, std::cout << "\tBegin initial sampling..." << std::flush);


    Real *weights    = new Real[m_M];
    unsigned int nbInsPart = m_M;


    #pragma omp parallel for
    for(unsigned int m=0; m<m_M; m++)
    {
        // Simulate a direction using initial density
        Direction u0 = m_importance.simulate(dir,m_importance.computeConcentration(dir,m_x0));

        // Move particle
        bool isInside = m_cloud[m].addToPath(u0.toVector()*m_stepSize, m_mask);

        // Set the initial weight of the particle
        if(isInside == 0 || isInside == 2) // outside mask or in exclusion zone
        {
            nbInsPart--;
            m_cloud[m].setActive(false);
        }
    } // for i

    Real w    = 1.0/(Real)m_M;

    #pragma omp parallel for
    for(unsigned int m=0; m<m_M; m++)
    {
        weights[m] = w;
        m_cloud[m].setWeight(w);
    }

//    this->saveCloudInVTK(label, m_k, m_x0);

    m_k++;

    Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Sequential sampling
    //


    while(nbInsPart > 0)
    {
        Display2(m_displayMode, std::cout << "\tBegin sampling " << m_k << "..." << std::flush);

        #pragma omp parallel for
        for(unsigned int m=0; m<m_M; m++)
        {
            if(m_cloud[m].isActive()) // m is active
            {
                // Informations about current particle
                Direction ukm1 = m_cloud[m].lastVector().toDirection();
                Point xk       = m_cloud[m].lastPoint();

                // Simulating next direction according importance density
                Direction mu = m_importance.computeMeanDirection(xk, ukm1);
                Real kappa   = mu.isNull() ? m_aPriori.getConcentration() : m_importance.computeConcentration(mu,xk);
                Direction uk = mu.isNull() ? m_importance.simulate(ukm1, kappa) : m_importance.simulate(mu, kappa);

                // Move particle
                m_cloud[m].addToPath(uk.toVector()*m_stepSize, m_mask);

                // Compute particle's weight
                if(mu.isNull())
                {
                    Real apriori    = m_aPriori.compute(uk, ukm1);
                    Real importance = m_importance.compute(uk, ukm1, kappa);
                    m_cloud[m].addLikelihood(0);

                    weights[m] = std::exp( std::log(m_cloud[m].weight()) + apriori - importance );
                }
                else // mu is not null
                {
                    Real likelihood = m_likelihood.compute(uk, xk, mu);
                    Real apriori    = m_aPriori.compute(uk, ukm1);
                    Real importance = m_importance.compute(uk, mu, kappa);
                    m_cloud[m].addLikelihood(likelihood);

                    weights[m] = std::exp( std::log(m_cloud[m].weight()) + likelihood + apriori - importance );
                }

                if(!std::isfinite(weights[m]) || std::isnan(weights[m]))
                    weights[m] = 0;
            }
        } // for i particles

        // Compute the sum of particles' weight
        Real sum   = 0;

        for(unsigned int m=0; m<m_M; m++)
        {
            if(m_cloud[m].isActive()) // m is active
            {
                sum += weights[m];
            }
        } // for each particle

        if(sum <= 0.0)
            sum = 1;
        assert(sum > 0.0);

        // Normalize particles' weights
        #pragma omp parallel for
        for(unsigned int m=0; m<m_M; m++)
        {
            if(m_cloud[m].isActive()) // m is active
            {
                assert(weights[m]/sum >= 0.0);
                m_cloud[m].setWeight(weights[m]/sum);
            }
        }


        // Compute sum square of particles' weights
        Real sumSquare = 0;
        for(unsigned int m=0; m<m_M; m++)
        {
            if(m_cloud[m].isActive()) // m is active
            {
                sumSquare += m_cloud[m].weight() * m_cloud[m].weight();
            }
        }

        // Compute resampling threshold
        Real ESS = 1.0 / sumSquare;

        if(!std::isfinite(ESS))
            nbInsPart = 0;
        else
        {
            // Resample if ESS is below a threshold
            // keeping proportionnality of weights
            if(ESS < m_epsilon/*(nbInsPart)*/)
            {
                Display2(m_displayMode, std::cout << " (Resampling, ESS = " << (ESS/nbInsPart)*100.0 << "%) " << std::flush);
                this->ResampleCloud(nbInsPart);
            }
            else // ESS >= m_epsilon*(nbInsPart)
                Display2(m_displayMode, std::cout << " (ESS = " << (ESS/nbInsPart)*100.0 << "%) " << std::flush);


            for(unsigned int m=0; m<m_M; m++)
            {
                if(m_cloud[m].isOutside() && m_cloud[m].isActive())
                {
                    m_cloud[m].setActive(false);
                    nbInsPart--;
                }
            }

        }

//        this->saveCloudInVTK(label, m_k, m_x0);

        m_k++;

        Display2(m_displayMode, std::cout << "done." << std::endl);
    } // while there are active particles or for k steps

    delete[] weights;

//    this->saveCloudInVTK(label, m_k-1, m_x0);

    m_dirNum++;
}

void ParticleFilter::ResampleCloud(unsigned int nbInsPart)
{
    assert(nbInsPart <= m_M);
    assert(nbInsPart > 0);

    Real cumul          = 0;
    Real *intervals     = new Real[nbInsPart];
    unsigned int *index = new unsigned int[nbInsPart];


    // Create proportionnal intervals
    // between 0 and 1 for inside particles
    unsigned int count = 0;
    for(unsigned int m=0; m<m_M; m++)
    {
        if(m_cloud[m].isActive()) // m is active
        {
            cumul += m_cloud[m].weight();
            intervals[count] = cumul;
            index[count] = m;
            count++;
        }
    } // for i


    Real w    = 1.0/(Real)nbInsPart;

    assert(w >= 0.0);

    // Get M particles from the cloud
    // keeping proportionnality
    // (multinomial resampling)
    for(unsigned int m=0; m<m_M; m++)
    {
        if(m_cloud[m].isActive()) // m is active
        {
            // Simulate x ~ U(0,1)
            Real x = (Real)rand() / (Real)RAND_MAX;

            bool found = false;
            unsigned int i = 0;

            do
            {
                if(x < intervals[i])
                    found = true;
                else
                    i++;
            } while(!found && i < nbInsPart);

            m_cloud[m].SetLastPoint(m_cloud[index[i]].lastPoint());
            m_cloud[m].SetLastVector(m_cloud[index[i]].lastVector());
            m_cloud[m].SetLastWeight(w);
        }
    } // for m


    // Clean space memory
    delete[] index;
    delete[] intervals;
}


void ParticleFilter::saveCloudInVTK(int label, unsigned int step, Point begin)
{
    //
    // Build VTK PolyData from particle's cloud
    //

    vtkSmartPointer<vtkPoints>       points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>     lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> weights = vtkSmartPointer<vtkDoubleArray>::New();

    vtkIdType pid;
    double valueWeight = 0;

    weights->SetNumberOfComponents(1);

    // create all points and lines
    for(std::vector<Particle>::iterator pIt = m_cloud.begin(); pIt != m_cloud.end(); pIt++)
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        ImageContinuousIndex cix;
        Image::PointType wx;
        unsigned int id = 0;

        Point x0 = pIt->getPoint(0);

        cix[0] = x0.x(); cix[1] = x0.y(); cix[2] = x0.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

        if(m_lps)
            pid = points->InsertNextPoint(wx[0], wx[1], wx[2]);
        else // ras
            pid = points->InsertNextPoint(-wx[0], -wx[1], wx[2]);

        valueWeight = pIt->getWeight(id++);
        weights->InsertNextTupleValue(&valueWeight);

        for(unsigned int i=1; i<=pIt->length(); i++)
        {
            Point xk = pIt->getPoint(i);

            Point xkm1 = pIt->getPoint(i-1);
            Vector vk  = pIt->getVector(i-1);
            Point p    = xk + vk*(-1);

            if(!(p == xkm1)) // points are note linked
            {
                cix[0] = p.x(); cix[1] = p.y(); cix[2] = p.z();
                m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

                if(m_lps)
                    pid = points->InsertNextPoint(wx[0], wx[1], wx[2]);
                else // ras
                    pid = points->InsertNextPoint(-wx[0], -wx[1], wx[2]);

                valueWeight = 1.0/m_M;
                weights->InsertNextTupleValue(&valueWeight);
            }

            cix[0] = xk.x(); cix[1] = xk.y(); cix[2] = xk.z();
            m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

            if(m_lps)
            pid = points->InsertNextPoint(wx[0], wx[1], wx[2]);
            else // ras
            pid = points->InsertNextPoint(-wx[0], -wx[1], wx[2]);

            valueWeight = pIt->getWeight(id++);
            weights->InsertNextTupleValue(&valueWeight);

            line->GetPointIds()->SetId(0, pid-1);
            line->GetPointIds()->SetId(1, pid);
            lines->InsertNextCell(line);
        } // for i
    } // for each particle

    // create vtk polydata object
    vtkSmartPointer<vtkPolyData> cloud = vtkSmartPointer<vtkPolyData>::New();
    cloud->SetPoints(points);
    cloud->SetLines(lines);
    cloud->GetPointData()->SetScalars(weights);


    //
    // Save VTK PolyData object into VTK file
    //

    std::stringstream filename;
    filename << "cloud-" << label << "-" << begin << "-" << m_dirNum << "-" << step << ".vtk";

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInput(cloud);
    writer->SetFileName(filename.str().c_str());
    writer->Write();
}

Particle ParticleFilter::GetMAP()
{
    Particle map(m_x0);

    if(m_k > 2)
    {
        unsigned int t = m_k-1;
        Real *delta = new Real[t*m_M];
        Real *psi   = new Real[(t-1)*m_M];


        // Compute delta and psi by dynamic programming

        #pragma omp parallel for
        for(unsigned int m=0; m<m_M; m++)
        {
            if(m_cloud[m].length() > 1)
            {
                assert(m_cloud[m].length() > 1);
                delta[Ind(0,m)] = m_aPriori.compute(m_cloud[m].getVector(1).toDirection(), m_cloud[m].getVector(0).toDirection()) + m_cloud[m].getLikelihood(0);
            }
            else
                delta[Ind(0,m)] = std::log(1.0/m_M);
        }


        for(unsigned int k=1; k<t; k++)
        {
            #pragma omp parallel for
            for(unsigned int m=0; m<m_M; m++)
            {
                if(k+1 < m_cloud[m].length())
                {
                    Real max          = MIN_REAL;
                    unsigned int imax = 0;

                    assert(k+1 <= m_cloud[m].length());
                    Point xkp1m = m_cloud[m].getPoint(k+1);

                    for(unsigned int i=0; i<m_M; i++)
                    {
                        if(k < m_cloud[i].length())
                        {
                            assert(k+1 <= m_cloud[i].length());
                            Point xkp1i   = m_cloud[i].getPoint(k+1);
                            Real distance = std::sqrt( (xkp1i.x()-xkp1m.x())*(xkp1i.x()-xkp1m.x()) + (xkp1i.y()-xkp1m.y())*(xkp1i.y()-xkp1m.y()) + (xkp1i.z()-xkp1m.z())*(xkp1i.z()-xkp1m.z()));

                            if(distance <= 0.01)
                            {
                                assert(k+1 < m_cloud[m].length());
                                assert(k < m_cloud[i].length());
                                Real tmp = delta[Ind(k-1,i)] + m_aPriori.compute(m_cloud[m].getVector(k+1).toDirection(), m_cloud[i].getVector(k).toDirection());

                                if(tmp > max)
                                {
                                    max  = tmp;
                                    imax = i;
                                }
                            }
                        }
                    }

                    if(EQUAL(max,MIN_REAL))
                    {
                        delta[Ind(k,m)] = m_cloud[m].getLikelihood(k);
                        psi[Ind(k-1,m)] = m;
                    }
                    else // max <> MIN_REAL
                    {
                        delta[Ind(k,m)] = m_cloud[m].getLikelihood(k) + max;
                        psi[Ind(k-1,m)] = imax;
                    }
                }
                else // k+1 >= m_cloud[m].length()
                {
                    if(k == 1)
                        delta[Ind(1,m)] = delta[Ind(0,m)];
                    else // k > 1
                        delta[Ind(k,m)] = delta[Ind(k-1,m)] + std::abs(delta[Ind(k-1,m)] - delta[Ind(k-2,m)]);

                    psi[Ind(k-1,m)] = m;
                }
            }
        }


        // Backtracking

        std::vector<Point> points;

        Real max          = MIN_REAL;
        unsigned int mmax = 0;

        for(unsigned int m=0; m<m_M; m++)
        {
            Real tmp = delta[Ind(t-1,m)];

            if(max < tmp)
            {
                max  = tmp;
                mmax = m;
            }
        }

        int len = m_cloud[mmax].length();
        points.push_back(m_cloud[mmax].getPoint(len));

        for(int k=t-2; k>=0; k--)
        {
            mmax = psi[Ind(k,mmax)];

            int len = m_cloud[mmax].length();

            if(k+1 <= len)
                points.push_back(m_cloud[mmax].getPoint(k+1));
            else // k+1 > len
                points.push_back(m_cloud[mmax].getPoint(len));
        }


        delete[] psi;
        delete[] delta;


        // Compute MAP estimate

        for(int k=points.size()-1; k>=0; k--)
        {
            Point p1 = map.lastPoint();
            Point p2 = points[k];
            Vector v(p2.x()-p1.x(), p2.y()-p1.y(), p2.z()-p1.z());
            map.addToPath(v*m_stepSize, m_mask);
            map.setWeight(1);
        }
    }


    return map;
}

double ParticleFilter::ComputeGFA(Point p)
{
    Matrix S = m_signal->signalAt(p);

    double average = 0;
    for(unsigned int i=0; i<S.Rows(); i++)
        average += S(i,0);
    average /= S.Rows();

    double std = 0;
    for(unsigned int i=0; i<S.Rows(); i++)
        std += (S(i,0)-average) * (S(i,0)-average);
    std /= S.Rows();
    std = std::sqrt(std);

    double rms = 0;
    for(unsigned int i=0; i<S.Rows(); i++)
        rms += S(i,0) * S(i,0);
    rms /= S.Rows();
    rms = std::sqrt(rms);

    return std / rms;
}

void ParticleFilter::ComputeFiber(Particle map1, Particle map2)
{
    // VTK structures
    m_fiber = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyLine>   line = vtkSmartPointer<vtkPolyLine>::New();
    vtkSmartPointer<vtkDoubleArray> gfa = vtkSmartPointer<vtkDoubleArray>::New();

    // Data
    vtkIdType pid;
    unsigned int id = 0;
    double valueGFA = 0;

    ImageContinuousIndex ci;
    Image::PointType wp;


    // Build fiber with the MAP estimate in world coordinates
    gfa->SetNumberOfComponents(1);
    line->GetPointIds()->SetNumberOfIds(map1.length() + map2.length()/* - 2*/);

    // First point
    Point p = map1.getPoint(map1.length());
    ci[0] = p.x(); ci[1] = p.y(); ci[2] = p.z();
    m_map->TransformContinuousIndexToPhysicalPoint(ci,wp);

    if(m_lps)
        pid = points->InsertNextPoint(wp[0], wp[1], wp[2]);
    else // ras
        pid = points->InsertNextPoint(-wp[0], -wp[1], wp[2]);

    valueGFA = this->ComputeGFA(p);
    gfa->InsertNextTuple(&valueGFA);
    line->GetPointIds()->SetId(id++,pid);

    // First part
    for(int k=map1.length()-2; k>=0; k--)
    {
        assert(k+2 <= map1.length());
        assert(k+1 <= map1.length());
        assert(k   <= map1.length());

        Point p1 = map1.getPoint(k+2);
        Point p2 = map1.getPoint(k+1);
        Point p3 = map1.getPoint(k);

        Image::PointType wp1, wp2, wp3;
        ImageContinuousIndex cip;

        cip[0] = p1.x(); cip[1] = p1.y(); cip[2] = p1.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp1);
        cip[0] = p2.x(); cip[1] = p2.y(); cip[2] = p2.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp2);
        cip[0] = p3.x(); cip[1] = p3.y(); cip[2] = p3.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp3);

        if(m_lps)
            pid = points->InsertNextPoint( (wp1[0]+wp2[0]+wp3[0])/3.0, (wp1[1]+wp2[1]+wp3[1])/3.0, (wp1[2]+wp2[2]+wp3[2])/3.0 );
        else // ras
            pid = points->InsertNextPoint( -(wp1[0]+wp2[0]+wp3[0])/3.0, -(wp1[1]+wp2[1]+wp3[1])/3.0, (wp1[2]+wp2[2]+wp3[2])/3.0 );

        valueGFA = this->ComputeGFA(Point( (p1.x()+p2.x()+p3.x())/3.0, (p1.y()+p2.y()+p3.y())/3.0, (p1.z()+p2.z()+p3.z())/3.0 ));
        gfa->InsertNextTuple(&valueGFA);
        line->GetPointIds()->SetId(id++,pid);
    }

    // Middle point
    if(map1.length() >= 1 && map2.length() >= 1)
    {
        assert(1 <= map1.length());
        assert(1 <= map2.length());

        Point p1 = map1.getPoint(1);
        Point p2 = map2.getPoint(0);
        Point p3 = map2.getPoint(1);

        Image::PointType wp1, wp2, wp3;
        ImageContinuousIndex cip;

        cip[0] = p1.x(); cip[1] = p1.y(); cip[2] = p1.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp1);
        cip[0] = p2.x(); cip[1] = p2.y(); cip[2] = p2.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp2);
        cip[0] = p3.x(); cip[1] = p3.y(); cip[2] = p3.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp3);

        if(m_lps)
            pid = points->InsertNextPoint( (wp1[0]+wp2[0]+wp3[0])/3.0, (wp1[1]+wp2[1]+wp3[1])/3.0, (wp1[2]+wp2[2]+wp3[2])/3.0 );
        else // ras
            pid = points->InsertNextPoint( -(wp1[0]+wp2[0]+wp3[0])/3.0, -(wp1[1]+wp2[1]+wp3[1])/3.0, (wp1[2]+wp2[2]+wp3[2])/3.0 );

        valueGFA = this->ComputeGFA(Point( (p1.x()+p2.x()+p3.x())/3.0, (p1.y()+p2.y()+p3.y())/3.0, (p1.z()+p2.z()+p3.z())/3.0 ));
        gfa->InsertNextTuple(&valueGFA);
        line->GetPointIds()->SetId(id++,pid);
    }


    // Last part
    for(unsigned int k=2; k<=map2.length(); k++)
    {
        assert(k-2 <= map2.length());
        assert(k-1 <= map2.length());
        assert(k   <= map2.length());

        Point p1 = map2.getPoint(k-2);
        Point p2 = map2.getPoint(k-1);
        Point p3 = map2.getPoint(k);

        Image::PointType wp1, wp2, wp3;
        ImageContinuousIndex cip;

        cip[0] = p1.x(); cip[1] = p1.y(); cip[2] = p1.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp1);
        cip[0] = p2.x(); cip[1] = p2.y(); cip[2] = p2.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp2);
        cip[0] = p3.x(); cip[1] = p3.y(); cip[2] = p3.z();
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp3);

        if(m_lps)
            pid = points->InsertNextPoint( (wp1[0]+wp2[0]+wp3[0])/3.0, (wp1[1]+wp2[1]+wp3[1])/3.0, (wp1[2]+wp2[2]+wp3[2])/3.0 );
        else // ras
            pid = points->InsertNextPoint( -(wp1[0]+wp2[0]+wp3[0])/3.0, -(wp1[1]+wp2[1]+wp3[1])/3.0, (wp1[2]+wp2[2]+wp3[2])/3.0 );

        valueGFA = this->ComputeGFA(Point( (p1.x()+p2.x()+p3.x())/3.0, (p1.y()+p2.y()+p3.y())/3.0, (p1.z()+p2.z()+p3.z())/3.0 ));
        gfa->InsertNextTuple(&valueGFA);
        line->GetPointIds()->SetId(id++,pid);
    }

    lines->InsertNextCell(line);
    m_fiber->SetPoints(points);
    m_fiber->SetLines(lines);
    m_fiber->GetPointData()->SetScalars(gfa);
}

void ParticleFilter::saveFiber(int label/*, unsigned int step, Point begin*/)
{
    std::stringstream filename;
    filename << "fiber-" << label << ".txt";

    std::fstream file(filename.str().c_str(), std::fstream::out);

    double p[3]; p[0] = -10000; p[1] = -10000; p[2] = -10000;
    vtkSmartPointer<vtkPoints> points = m_fiber->GetPoints();
    for(vtkIdType id=0; id<points->GetNumberOfPoints(); id++)
    {
        double tmp[3];
        points->GetPoint(id, tmp);

        if(!(EQUAL(p[0], tmp[0]) && EQUAL(p[1], tmp[1]) && EQUAL(p[2], tmp[2])))
            file << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;

        p[0] = tmp[0]; p[1] = tmp[1]; p[2] = tmp[2];
    }

    file.close();
}

vtkSmartPointer<vtkPolyData> ParticleFilter::GetFiber()
{
    return m_fiber;
}

void ParticleFilter::ComputeMap()
{
    for(std::vector<Particle>::iterator particle = m_cloud.begin(); particle != m_cloud.end(); particle++)
    {

        for(unsigned int k=0; k<=particle->length(); k++)
        {
            // Get data
            assert(k <= particle->length());
            Point p = particle->getPoint(k);
//            Real  w = (k == 0) ? 1.0/(Real)m_M : particle->getWeight(k-1);
            Real w  = 1;

            // Tri-linear interpolation
            short x = (short)std::floor(p.x());
            short y = (short)std::floor(p.y());
            short z = (short)std::floor(p.z());

            Real wx = (Real)(p.x() - x);
            Real wy = (Real)(p.y() - y);
            Real wz = (Real)(p.z() - z);

            // Set particle's weight at particle's position
            Image::IndexType index;

            index[0] = x; index[1] = y; index[2] = z;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + (1-wz)*(1-wx)*(1-wy) * w);

            index[0] = x; index[1] = y+1; index[2] = z;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + (1-wz)*(1-wx)*wy * w);

            index[0] = x+1; index[1] = y; index[2] = z;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + (1-wz)*wx*(1-wy) * w);

            index[0] = x+1; index[1] = y+1; index[2] = z;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + (1-wz)*wx*wy * w);

            index[0] = x; index[1] = y; index[2] = z+1;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + wz*(1-wx)*(1-wy) * w);

            index[0] = x; index[1] = y+1; index[2] = z+1;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + wz*(1-wx)*wy * w);

            index[0] = x+1; index[1] = y; index[2] = z+1;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + wz*wx*(1-wy) * w);

            index[0] = x+1; index[1] = y+1; index[2] = z+1;
            if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < m_map->GetLargestPossibleRegion().GetSize(0) &&
               0 <= (unsigned int)index[1] && (unsigned int)index[1] < m_map->GetLargestPossibleRegion().GetSize(1) &&
               0 <= (unsigned int)index[2] && (unsigned int)index[2] < m_map->GetLargestPossibleRegion().GetSize(2))
                m_map->SetPixel(index, m_map->GetPixel(index) + wz*wx*wy * w);
        }
    }
}

void ParticleFilter::saveConnectionMap(int label, Point begin)
{
    // Normalize image
    ImageIterator it(m_map, m_map->GetLargestPossibleRegion());

    Real max = 0;
    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
        if(max < it.Get())
            max = it.Get();
    }

    for(it.GoToBegin(); !it.IsAtEnd(); ++it)
        it.Set( it.Get() / max );


    //
    // Write data
    //

    // Write connectivity map
    ImageWriter::Pointer writer = ImageWriter::New();

    std::stringstream filename;
    filename << "map-" << label << "-" << begin << "-" << m_dirNum << "-" << m_k << ".nii.gz";
    writer->SetFileName(filename.str().c_str());
    writer->SetInput(m_map);

    try
    {
        writer->Update();
    }
    catch(itk::ImageFileWriterException &err)
    {
        std::cout << "Error: " << std::endl;
        std::cout << err << std::endl;
    }
}

Image::Pointer ParticleFilter::GetConnectionMap()
{
    return m_map;
}

void ParticleFilter::SetLPSOn()
{
    m_lps = true;
}

} // namespace btk

