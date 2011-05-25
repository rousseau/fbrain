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

#include "btkDTFPParticleFilter.h"


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

// OpenMP includes
#include "omp.h"


namespace btk
{

DTFPParticleFilter::DTFPParticleFilter(DTFPSignal *model, DTFPAPrioriDensity aPriori, DTFPLikelihoodDensity likelihood, DTFPImportanceDensity importance,
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

DTFPParticleFilter::DTFPParticleFilter(DTFPSignal *model, DTFPAPrioriDensity aPriori, DTFPLikelihoodDensity likelihood, DTFPImportanceDensity importance,
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
    m_dirNum    = 1;
    m_kx        = m_stepSize/m_spacing[0];
    m_ky        = m_stepSize/m_spacing[1];
    m_kz        = m_stepSize/m_spacing[2];
    m_lps       = false;

    m_vStepSize.push_back(m_kx);
    m_vStepSize.push_back(m_ky);
    m_vStepSize.push_back(m_kz);

    m_mask = mask;


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


    Display2(m_displayMode, std::cout << "\tdone." << std::endl);
}

void DTFPParticleFilter::run(int label)
{
    //
    // Initialization
    //

    // Get max direction à starting point
    itk::DiffusionTensor3D<Real> tensor = m_model->DiffusionTensorAt(m_x0);
    itk::FixedArray<Real,3> eigenValues;
    itk::Matrix<Real,3,3> eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues,eigenVectors);
//    Pr(eigenValues[0]); Pr(eigenValues[1]); Pr(eigenValues[2]);
//    Pr(eigenVectors(0,0)); Pr(eigenVectors(0,1)); Pr(eigenVectors(0,2));
//    Pr(eigenVectors(1,0)); Pr(eigenVectors(1,1)); Pr(eigenVectors(1,2));
//    Pr(eigenVectors(2,0)); Pr(eigenVectors(2,1)); Pr(eigenVectors(2,2));
//    Pr(tensor.GetFractionalAnisotropy());

    Direction maxDir = Vector(eigenVectors(2,0),eigenVectors(2,1),eigenVectors(2,2)).toDirection();

    //
    // Filtering
    //

    // Compute Symetric direction
    Direction symDir(
        M_PI - maxDir.theta(),
        (maxDir.phi() < M_PI) ? maxDir.phi() + M_PI : maxDir.phi() - M_PI
    );

    this->run(label, maxDir);
    this->ComputeMap();
    DTFPParticle map1 = this->GetMAP();

    this->run(label, symDir);
    this->ComputeMap();
    DTFPParticle map2 = this->GetMAP();

    this->ComputeFiber(map1,map2);
    this->saveFiber(label,m_k,m_x0);
}


void DTFPParticleFilter::run(int label, Direction dir)
{
    m_k     = 0;
    m_cloud = std::vector<DTFPParticle>(m_M, DTFPParticle(m_x0));


    //
    // Initial sampling (vMF in mean direction dir and a priori kappa)
    //

    Display2(m_displayMode, std::cout << "\tBegin initial sampling..." << std::flush);


    unsigned int nbInsPart = m_M;
    Real *weights = new Real[m_M];


    #pragma omp parallel for
    for(unsigned int m=0; m<m_M; m++)
    {
        // Simulate a direction using initial density
        Real kappa   = m_importance.computeConcentration(dir,m_x0);
        Direction u0 = m_importance.simulate(dir,kappa);

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
    Real logw = std::log(w);

    #pragma omp parallel for
    for(unsigned int m=0; m<m_M; m++)
    {
        weights[m] = logw;
        m_cloud[m].setWeight(w);
    }

//    this->saveCloudInVTK(label, m_k, m_x0);


    m_k++;

    Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Sequential sampling
    //


    while(nbInsPart > 0/*m_k < 61*/)
    {
        Display2(m_displayMode, std::cout << "\tBegin sampling " << m_k << "..." << std::flush);

        #pragma omp parallel for
        for(unsigned int m=0; m<m_M; m++)
        {
            if(m_cloud[m].isActive() && !m_cloud[m].isOutside()) // m is active
            {
                // Informations about current particle
                Direction ukm1 = m_cloud[m].lastVector().toDirection();
                Point xk       = m_cloud[m].lastPoint();

                // Simulating next direction according importance density
                Direction mu = m_importance.computeMeanDirection(xk, ukm1);
                Real kappa   = m_importance.computeConcentration(mu,xk);
                Direction uk = m_importance.simulate(mu, kappa);

                // Move particle
                char isInside = m_cloud[m].addToPath(uk.toVector()*m_stepSize, m_mask);
//if(isInside != 0)
//{
                // Compute particle's weight
                Real likelihood = m_likelihood.compute(uk, xk, mu);
                Real apriori    = m_aPriori.compute(uk, ukm1);
                Real importance = m_importance.compute(uk, mu, kappa);

                weights[m] = m_cloud[m].weight() * std::exp(likelihood + apriori - importance);
//}
//else
//    weights[m] = 0;

                assert(weights[m] >= 0.0);
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
                sumSquare += m_cloud[m].weight() * m_cloud[m].weight();
        }

        // Compute resampling threshold
        Real ESS = 1.0 / sumSquare;

        if(!std::isfinite(ESS))
            nbInsPart = 0;
        else
        {
            // Resample if ESS is below a threshold
            // keeping proportionnality of weights
            if(ESS < m_epsilon*(nbInsPart))
            {
                Display2(m_displayMode, std::cout << " (Resampling, ESS = " << (ESS/nbInsPart)*100.0 << "%) " << std::flush);
                this->ResampleCloud(nbInsPart,weights);
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
//Pr("ok");
//        this->saveCloudInVTK(label, m_k, m_x0);
//        this->ComputeFiber(this->GetMAP(),this->GetMAP());
//        this->ComputeMap();
//        this->saveFiber(label, m_k-1, m_x0);
//        this->saveConnectionMap(label,m_x0);

        m_k++;

        Display2(m_displayMode, std::cout << "done." << std::endl);
    } // while there are active particles or for k steps


    delete[] weights;

//    this->saveCloudInVTK(label, m_k-1, m_x0);

    m_dirNum++;
}

void DTFPParticleFilter::ResampleCloud(unsigned int nbInsPart, Real *weights)
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


    std::vector<DTFPParticle> cloud;
    Real w    = 1.0/(Real)nbInsPart;

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

            cloud.push_back(m_cloud[index[i]]);
            cloud.back().SetLastWeight(w);
        }
        else    // not in mask
            cloud.push_back(m_cloud[m]);
    } // for m

    m_cloud = cloud;


    // Clean space memory
    delete[] index;
    delete[] intervals;
}


void DTFPParticleFilter::saveCloudInVTK(int label, unsigned int step, Point begin)
{
    //
    // Build VTK PolyData from particle's cloud
    //

    vtkSmartPointer<vtkPoints>      points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>    lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();

    vtkIdType pid[1];
    colors->SetNumberOfComponents(1);

    // create all points and lines
    for(std::vector<DTFPParticle>::iterator pIt = m_cloud.begin(); pIt != m_cloud.end(); pIt++)
    {
//        if(pIt->isActive())
//        {
            ImageContinuousIndex cix;
            Image::PointType wx;

            Point x0 = pIt->getPoint(0);

            cix[0] = x0.x(); cix[1] = x0.y(); cix[2] = x0.z();
            m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

            if(m_lps)
                points->InsertNextPoint(wx[0], wx[1], wx[2]);
            else // ras
                points->InsertNextPoint(-wx[0], -wx[1], wx[2]);

            double color[1]; color[0] = pIt->weight();
            colors->InsertNextTupleValue(color);

            for(unsigned int i=1; i<=pIt->length(); i++)
            {
                vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
                Point xk = pIt->getPoint(i);

                cix[0] = xk.x(); cix[1] = xk.y(); cix[2] = xk.z();
                m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

                if(m_lps)
                    pid[0] = points->InsertNextPoint(wx[0], wx[1], wx[2]);
                else // ras
                    pid[0] = points->InsertNextPoint(-wx[0], -wx[1], wx[2]);

                Point xkm1 = pIt->getPoint(i-1);
                Vector vk  = pIt->getVector(i-1);
                Point p    = xk + vk*(-1);

                if(!(p == xkm1)) // points are note linked
                {
                    cix[0] = p.x(); cix[1] = p.y(); cix[2] = p.z();
                    m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

                    if(m_lps)
                        pid[0] = points->InsertNextPoint(wx[0], wx[1], wx[2]);
                    else // ras
                        pid[0] = points->InsertNextPoint(-wx[0], -wx[1], wx[2]);

                    cix[0] = xk.x(); cix[1] = xk.y(); cix[2] = xk.z();
                    m_map->TransformContinuousIndexToPhysicalPoint(cix, wx);

                    if(m_lps)
                        pid[0] = points->InsertNextPoint(wx[0], wx[1], wx[2]);
                    else // ras
                        pid[0] = points->InsertNextPoint(-wx[0], -wx[1], wx[2]);
                }

                line->GetPointIds()->SetId(0, pid[0]-1);
                line->GetPointIds()->SetId(1, pid[0]);
                lines->InsertNextCell(line);

                colors->InsertNextTupleValue(color);
            } // for i
//        }
    } // for each particle

    // create vtk polydata object
    vtkSmartPointer<vtkPolyData> cloud = vtkSmartPointer<vtkPolyData>::New();
    cloud->SetPoints(points);
    cloud->SetLines(lines);
    cloud->GetPointData()->SetScalars(colors);


    //
    // Save VTK PolyData object into VTK file
    //

    std::stringstream filename;
    filename << "cloud-" << label << "-" << begin << "-" << m_dirNum << "-" << step << ".vtk";

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInput(cloud);
    writer->SetFileName(filename.str().c_str());
    writer->SetFileTypeToBinary();
    writer->Write();
}

DTFPParticle DTFPParticleFilter::GetMAP()
{
    DTFPParticle map(m_x0);
    Real max = m_cloud[0].weight();

    for(std::vector<DTFPParticle>::iterator it = m_cloud.begin(); it != m_cloud.end(); it++)
    {
        Real weight = it->weight();

        if(max < weight)
        {
            max = weight;
            map = *it;
        }
    }
//Pr(max);
//Pr(map.length());
    return map;
}

void DTFPParticleFilter::ComputeFiber(DTFPParticle map1, DTFPParticle map2)
{
    // VTK structures
    vtkSmartPointer<vtkAppendPolyData> append = vtkSmartPointer<vtkAppendPolyData>::New();

    vtkSmartPointer<vtkPolyData> fiber1       = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points1        = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines1      = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> fiber2       = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points2        = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines2      = vtkSmartPointer<vtkCellArray>::New();

    vtkIdType pid[1];


    // Build fiber with the MAP estimate in world coordinates

    Point x0 = map1.getPoint(0);
    ImageContinuousIndex cix0;
    cix0[0] = x0.x(); cix0[1] = x0.y(); cix0[2] = x0.z();
    Image::PointType wx0;
    m_map->TransformContinuousIndexToPhysicalPoint(cix0, wx0);

    if(m_lps)
        points1->InsertNextPoint(wx0[0], wx0[1], wx0[2]);
    else // ras
        points1->InsertNextPoint(-wx0[0], -wx0[1], wx0[2]);

    for(unsigned int k=1; k<map1.length(); k++)
    {
        Point p = map1.getPoint(k);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

        ImageContinuousIndex cip;
        cip[0] = p.x(); cip[1] = p.y(); cip[2] = p.z();
        Image::PointType wp;
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp);

        if(m_lps)
            pid[0] = points1->InsertNextPoint(wp[0], wp[1], wp[2]);
        else // ras
            pid[0] = points1->InsertNextPoint(-wp[0], -wp[1], wp[2]);

        line->GetPointIds()->SetId(0, pid[0]-1);
        line->GetPointIds()->SetId(1, pid[0]);
        lines1->InsertNextCell(line);
    }

    fiber1->SetPoints(points1);
    fiber1->SetLines(lines1);


    x0 = map2.getPoint(0);
    cix0[0] = x0.x(); cix0[1] = x0.y(); cix0[2] = x0.z();
    m_map->TransformContinuousIndexToPhysicalPoint(cix0, wx0);

    if(m_lps)
        points2->InsertNextPoint(wx0[0], wx0[1], wx0[2]);
    else // ras
        points2->InsertNextPoint(-wx0[0], -wx0[1], wx0[2]);

    for(unsigned int k=1; k<map2.length(); k++)
    {
        Point p = map2.getPoint(k);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

        ImageContinuousIndex cip;
        cip[0] = p.x(); cip[1] = p.y(); cip[2] = p.z();
        Image::PointType wp;
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp);

        if(m_lps)
            pid[0] = points2->InsertNextPoint(wp[0], wp[1], wp[2]);
        else // ras
            pid[0] = points2->InsertNextPoint(-wp[0], -wp[1], wp[2]);

        line->GetPointIds()->SetId(0, pid[0]-1);
        line->GetPointIds()->SetId(1, pid[0]);
        lines2->InsertNextCell(line);
    }

    fiber2->SetPoints(points2);
    fiber2->SetLines(lines2);


    // create vtk polydata object
    append->AddInput(fiber1);
    append->AddInput(fiber2);
    append->Update();
    m_fiber = append->GetOutput();
}

void DTFPParticleFilter::saveFiber(int label, unsigned int step, Point begin)
{
    // Save VTK PolyData object into VTK file
    std::stringstream filename;
    filename << "fiber-" << label << "-" << begin << "-" << m_dirNum << "-" << step << ".vtk";

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInput(m_fiber);
    writer->SetFileName(filename.str().c_str());
    writer->SetFileTypeToBinary();
    writer->Write();
}

vtkSmartPointer<vtkPolyData> DTFPParticleFilter::GetFiber()
{
    return m_fiber;
}

void DTFPParticleFilter::ComputeMap()
{
    for(std::vector<DTFPParticle>::iterator particle = m_cloud.begin(); particle != m_cloud.end(); particle++)
    {

        for(unsigned int k=0; k<=particle->length(); k++)
        {
            // Get data
            Point p = particle->getPoint(k);
//            Real  w = (k == 0) ? 1.0/(Real)m_M : particle->getWeight(k-1);
            Real w = 1;

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

void DTFPParticleFilter::saveConnectionMap(int label, Point begin)
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

Image::Pointer DTFPParticleFilter::GetConnectionMap()
{
    return m_map;
}

void DTFPParticleFilter::SetLPSOn()
{
    m_lps = true;
}

} // namespace btk
