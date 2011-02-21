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
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkColorTransferFunction.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkAppendPolyData.h"

// ITK includes
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPoint.h"
#include "itkContinuousIndex.h"

// OpenMP includes
#include "omp.h"

//#define NOPARA
#define NODISPLAY


namespace btk
{

ParticleFilter::ParticleFilter(SHModel *model, APrioriDensity aPriori, LikelihoodDensity likelihood, ImportanceDensity importance,
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

    std::cout << "\t\tFilter will use " << M << " particles." << std::endl;
    std::cout << "\t\tResampling treshold is set as " << epsilon << "." << std::endl;
    std::cout << "\t\tThe moving step is fixed at " << m_stepSize << " of voxel." << std::endl;
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

ParticleFilter::ParticleFilter(SHModel *model, APrioriDensity aPriori, LikelihoodDensity likelihood, ImportanceDensity importance,
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


    Display2(m_displayMode, std::cout << "\t\tFilter will use " << M << " particles." << std::endl);
    Display2(m_displayMode, std::cout << "\t\tResampling treshold is set as " << epsilon << "." << std::endl);
    Display2(m_displayMode, std::cout << "\t\tThe moving step is fixed at " << m_stepSize << " of voxel." << std::endl);


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

    #ifndef NOPARA
    #pragma omp parallel for
    #endif // NOPARA
    for(unsigned int i=0; i<m_M; i++)
    {
        // Simulate a direction using initial density
        Direction u0 = m_aPriori.simulate(dir);

        // Move particle
        Vector v0     = u0.toVector() * m_stepSize;
        bool isInside = m_cloud[i].addToPath(v0, m_mask);

        // Set the initial weight of the particle
        if(isInside)
            m_cloud[i].setWeight(1.0/(Real)m_M);
        else
            m_cloud[i].setWeight(0);
    } // for i

    m_k++;

    Display2(m_displayMode, std::cout << "done." << std::endl);


    //
    // Next sampling
    //

    unsigned int nbOfActiveParticles = m_M;
    Real *weights = new Real[m_M];

    for(unsigned int m=0; m<m_M; m++)
        weights[m] = 0;


    while(nbOfActiveParticles > 0)
    {
        Display2(m_displayMode, std::cout << "\tBegin sampling " << m_k << "..." << std::flush);

        #ifndef NODISPLAY
            std::cerr << std::endl << std::endl << "step " << m_k << std::endl;
        #endif // NODISPLAY

        #ifndef NOPARA
        #pragma omp parallel for
        #endif // NOPARA
        for(unsigned int i=0; i<m_M; i++)
        {
            if(m_cloud[i].isActive())
            {
                // Informations about current particle
                Direction ukm1 = m_cloud[i].lastVector().toDirection();
                Point xk       = m_cloud[i].lastPoint();

                // Simulating next direction according importance density
                Direction mu = m_importance.computeMeanDirection(xk, ukm1);
                Real kappa   = m_importance.computeConcentration(mu, xk);
                Direction uk = m_importance.simulate(mu, kappa);

                // Move particle
                Vector vk = uk.toVector() * m_stepSize;

                bool isInside = m_cloud[i].addToPath(vk, m_mask);

                if(isInside)
                {
                    // Compute particle's weight
                    Real likelihood = m_likelihood.compute(uk, xk, mu);
                    Real apriori    = m_aPriori.compute(uk, ukm1);
                    Real importance = m_importance.compute(uk, mu, kappa);

                    weights[i] = std::log(m_cloud[i].weight()) + likelihood + apriori - importance;
                }
                else
                {
                    nbOfActiveParticles--;
                    weights[i] = 0;
                }
            } // particle not active
            else
                weights[i] = 0;

            #ifndef NODISPLAY
                std::cerr << "_w[" << i << "] = " << weights[i] << std::endl;
            #endif // NODISPLAY
        } // for i particles

        if(nbOfActiveParticles > 0)
        {
            // Search minimal weight
            Real min = 0;

            for(unsigned int i=0; i<m_M; i++)
            {
                if(m_cloud[i].isActive() && std::isfinite(weights[i]) && min > weights[i])
                    min = weights[i];
            } // for each particle

            #ifndef NODISPLAY
                std::cerr << std::endl << "min = " << min << std::endl;
            #endif // NODISPLAY

            // Compute the sum of particles' weight (after an interval shift)
            Real sum   = 0;
            Real shift = (min < 0) ? -min+1 : 0;

            #ifndef NODISPLAY
                std::cerr << "shift = " << shift << std::endl << std::endl;
            #endif // NODISPLAY

            for(unsigned int i=0; i<m_M; i++)
            {
                if(m_cloud[i].isActive())
                {
                    if(std::isfinite(weights[i]))
                        weights[i] += shift;
                    else // infinite number
                        weights[i] = 0;

                    #ifndef NODISPLAY
                        std::cerr << "w[" << i << "] = " << weights[i] << std::endl;
                    #endif // NODISPLAY

                    sum += weights[i];
                }
            } // for each particle

            #ifndef NODISPLAY
                std::cerr << std::endl << "sum = " << sum << std::endl << std::endl;
            #endif // NODISPLAY

            // Normalize particles' weights
            #ifndef NOPARA
            #pragma omp parallel for
            #endif // NOPARA
            for(unsigned int i=0; i<m_M; i++)
            {
                m_cloud[i].setWeight(weights[i]/sum);

                #ifndef NODISPLAY
                    std::cerr << "~w[" << i << "] = " << m_cloud[i].weight() << std::endl;
                #endif // NODISPLAY
            }


            // Compute sum square of particles' weights
            Real sumSquare = 0;
            for(unsigned int i=0; i<m_M; i++)
            {
                sumSquare += m_cloud[i].weight() * m_cloud[i].weight();
            }



            // Compute resampling threshold
            Real ESS = 1. / sumSquare;

            if(!std::isfinite(ESS))
                nbOfActiveParticles = 0;
            else
            {

                #ifndef NODISPLAY
                    std::cerr << std::endl << "actives particles remaining = " << nbOfActiveParticles << std::endl << std::endl;
                #endif // NODISPLAY

                // Resample if ESS is below a threshold
                // keeping proportionnality of weights
                if(ESS < m_epsilon*m_M)
                {
                    Display2(m_displayMode, std::cout << " (Resampling, ESS = " << ESS << ") " << std::flush);

                    Real cumul = 0;

                    Real *intervals = new Real[m_M-1];


                    // Create proportionnal intervals
                    // between 0 and 1
                    for(unsigned int i=0; i<m_M-1; i++)
                    {
                        cumul += m_cloud[i].weight();
                        intervals[i] = cumul;
                    } // for i


                    nbOfActiveParticles = m_M;


                    // Get M particles from the cloud
                    // keeping proportionnality
                    // (multinomial resampling)
                    for(unsigned int m=0; m<m_M; m++)
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
                        } while(!found && i < m_M-1);

                        m_cloud[m].SetLastPoint(m_cloud[i].lastPoint());
                        m_cloud[m].SetLastVector(m_cloud[i].lastVector());
                        m_cloud[m].setWeight(1.0/(Real)m_M);

                        if(m_cloud[i].isActive())
                            m_cloud[m].SetActive();
                        else
                        {
                            m_cloud[m].SetInactive();
                            nbOfActiveParticles--;
                        }
                    } // for m


                    // Clean space memory
                    delete[] intervals;

                } // ESS < m_epsilon
                else
                    Display2(m_displayMode, std::cout << " (ESS = " << ESS << ") " << std::flush);
            }

            m_k++;

        } // endif nbOfActiveParticles > 0

        Display2(m_displayMode, std::cout << "done." << std::endl);
    } // for k steps

    delete[] weights;

    m_dirNum++;
}

void ParticleFilter::saveCloudInVTK(int label, unsigned int step, Point begin)
{
    //
    // Build VTK PolyData from particle's cloud
    //

    vtkSmartPointer<vtkPoints> points      = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> lines    = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> colors = vtkSmartPointer<vtkDoubleArray>::New();

    vtkIdType pid[1];

    // color map (search min and max values)
    Real max = m_cloud[0].weight();
    Real min = m_cloud[0].weight();

    std::vector<Particle>::iterator it;
    for(it = m_cloud.begin(); it != m_cloud.end(); it++)
    {
        if(it->weight() > max)
            max = it->weight();

        if(it->weight() < min)
            min = it->weight();
    } // for it

    vtkSmartPointer<vtkColorTransferFunction> map = vtkSmartPointer<vtkColorTransferFunction>::New();
    map->AddRGBPoint(min, 1, 0, 0);
    map->AddRGBPoint(max, 0, 0, 1);


    // create all points and lines
    std::vector<Particle>::iterator pIt;
    for(pIt = m_cloud.begin(); pIt != m_cloud.end(); pIt++)
    {
        Point x0 = pIt->getPoint(0);
        points->InsertNextPoint(x0.x()*m_spacing[0]+ m_origin[0], x0.y()*m_spacing[1] + m_origin[1], x0.z()*m_spacing[2] + m_origin[2]);

        double color[3];
        map->GetColor(pIt->weight(), color);
        colors->InsertNextTupleValue(color);

        for(unsigned int i=1; i<pIt->length()+1; i++)
        {
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            Point p = pIt->getPoint(i);

            pid[0] = points->InsertNextPoint(p.x()*m_spacing[0] + m_origin[0], p.y()*m_spacing[1] + m_origin[1], p.z()*m_spacing[2] + m_origin[2]);
            line->GetPointIds()->SetId(0, pid[0]-1);
            line->GetPointIds()->SetId(1, pid[0]);
            lines->InsertNextCell(line);

            map->GetColor(pIt->weight(), color);
            colors->InsertNextTupleValue(color);
        } // for i
    } // for pIt

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

Particle ParticleFilter::GetMAP()
{
//*
    // MAP estimate of probablity law
    // Search the particle with maximal importance weight
    Particle map = m_cloud[0];
    Real max = m_cloud[0].weight();

    for(unsigned int m=1; m<m_cloud.size(); m++)
    {
        if(max < m_cloud[m].weight())
        {
            max = m_cloud[m].weight();
            map = m_cloud[m];
        }
    } // for each particle in cloud

    return map;
//*/

/*
	std::vector<Point> map;

	for(unsigned int k=0; k<m_k; k++)
	{
		unsigned int imap = 0;
		Real         max  = 0;
		unsigned int num  = 0;

		for(unsigned int m=0; m<m_M; m++)
		{
			if(m_cloud[m].length() <= k)
			{
				Real tmp = m_cloud[m].getWeight(k);

				if(max < tmp)
				{
					max  = tmp;
					imap = m;
				}

				num++;
			}
		} // for each particle m

		if(num > 0)
			map.push_back(m_cloud[imap].getPoint(k));
	} // for each step k

	return map;
*/
}

void ParticleFilter::ComputeFiber(Particle map1, Particle map2)
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
    points1->InsertNextPoint(wx0[0], wx0[1], wx0[2]);

    for(unsigned int k=1; k<map1.length(); k++)
    {
        Point p = map1.getPoint(k);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

        ImageContinuousIndex cip;
        cip[0] = p.x(); cip[1] = p.y(); cip[2] = p.z();
        Image::PointType wp;
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp);
        pid[0] = points1->InsertNextPoint(wp[0], wp[1], wp[2]);

        line->GetPointIds()->SetId(0, pid[0]-1);
        line->GetPointIds()->SetId(1, pid[0]);
        lines1->InsertNextCell(line);
    }

    fiber1->SetPoints(points1);
    fiber1->SetLines(lines1);


    x0 = map1.getPoint(0);
    cix0[0] = x0.x(); cix0[1] = x0.y(); cix0[2] = x0.z();
    m_map->TransformContinuousIndexToPhysicalPoint(cix0, wx0);
    points2->InsertNextPoint(wx0[0], wx0[1], wx0[2]);

    for(unsigned int k=1; k<map2.length(); k++)
    {
        Point p = map1.getPoint(k);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

        ImageContinuousIndex cip;
        cip[0] = p.x(); cip[1] = p.y(); cip[2] = p.z();
        Image::PointType wp;
        m_map->TransformContinuousIndexToPhysicalPoint(cip, wp);
        pid[0] = points2->InsertNextPoint(wp[0], wp[1], wp[2]);

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

void ParticleFilter::saveFiber(int label, unsigned int step, Point begin)
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

vtkSmartPointer<vtkPolyData> ParticleFilter::GetFiber()
{
    return m_fiber;
}

void ParticleFilter::ComputeMap()
{
    // Accumulator map
    Image::Pointer accumulator = Image::New();
    accumulator->SetRegions(m_map->GetLargestPossibleRegion().GetSize());
    accumulator->SetOrigin(m_map->GetOrigin());
    accumulator->SetSpacing(m_map->GetSpacing());
    accumulator->SetDirection(m_map->GetDirection());
    accumulator->Allocate();


    //
    // Data extraction
    //

    for(unsigned int k=0; k<=m_k; k++)
    {
        accumulator->FillBuffer(0);

        for(unsigned int m=0; m<m_M; m++)
        {
            if(k <= m_cloud[m].length())
            {
                // Get data
                Point p = m_cloud[m].getPoint(k);
                Real  w = m_cloud[m].getWeight(k);
//                Real  w = 1;

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
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + (1-wz)*(1-wx)*(1-wy) * w);

                index[0] = x; index[1] = y+1; index[2] = z;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + (1-wz)*(1-wx)*wy * w);

                index[0] = x+1; index[1] = y; index[2] = z;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + (1-wz)*wx*(1-wy) * w);

                index[0] = x+1; index[1] = y+1; index[2] = z;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + (1-wz)*wx*wy * w);

                index[0] = x; index[1] = y; index[2] = z+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + wz*(1-wx)*(1-wy) * w);

                index[0] = x; index[1] = y+1; index[2] = z+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + wz*(1-wx)*wy * w);

                index[0] = x+1; index[1] = y; index[2] = z+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + wz*wx*(1-wy) * w);

                index[0] = x+1; index[1] = y+1; index[2] = z+1;
                if(0 <= (unsigned int)index[0] && (unsigned int)index[0] < accumulator->GetLargestPossibleRegion().GetSize(0) &&
                   0 <= (unsigned int)index[1] && (unsigned int)index[1] < accumulator->GetLargestPossibleRegion().GetSize(1) &&
                   0 <= (unsigned int)index[2] && (unsigned int)index[2] < accumulator->GetLargestPossibleRegion().GetSize(2))
                    accumulator->SetPixel(index, accumulator->GetPixel(index) + wz*wx*wy * w);
            }
        } // for each particle

        // Add accumulator to connectivity map
        ImageIterator in(accumulator, accumulator->GetLargestPossibleRegion());
        ImageIterator out(m_map, m_map->GetLargestPossibleRegion());

        for(in.GoToBegin(), out.GoToBegin(); !in.IsAtEnd() && !out.IsAtEnd(); ++in, ++out)
            out.Set(out.Get() + in.Get());
    } // for each step
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

} // namespace btk

