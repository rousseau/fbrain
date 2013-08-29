/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 04/01/2013
  Author(s): Julien Pontabry (pontabry@unistra.fr)

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

==========================================================================*/


#include "btkParticleFilteringTractographyAlgorithm.h"


// VTK includes
#include "vtkPolyLine.h"
#include "vtkCellArray.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkDoubleArray.h"
#include "vtkLine.h"
#include "vtkPolyDataWriter.h"
#include "vtkPointData.h"

#include "sstream"
#include "cfloat"

// Local includes
#include "btkImageHelper.h"


// Definitions
//#define INDEX(i,j) (m_NumberOfParticles*(i) + (j))


namespace btk
{

ParticleFilteringTractographyAlgorithm::ParticleFilteringTractographyAlgorithm() : Superclass(), m_NumberOfParticles(200), m_ParticleStepSize(0.5), m_ResamplingThreshold(0.05*m_NumberOfParticles), m_CurveConstraint(30.0), m_InitialWeight(1.0/m_NumberOfParticles), m_ImportanceDensity(), m_LikelihoodDensity(), m_PriorDensity(m_CurveConstraint), m_ProbabilityMap(NULL)
{
    // ----
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::Initialize()
{
    // Probability densities
    m_ImportanceDensity = ImportanceDensity(m_DiffusionModel, m_DiffusionSignal, m_ThresholdAngle);
    m_LikelihoodDensity = LikelihoodDensity(m_DiffusionSignal, m_DiffusionModel);

    // Probability map
    m_ProbabilityMap = btk::ImageHelper< DiffusionSignal,ProbabilityMap >::CreateNewImageFromPhysicalSpaceOf(m_DiffusionSignal);

    // Set number of threads to 1 (multi-threading is supported by opem mp instead)
    //this->SetNumberOfThreads(1); // NOTE : now the seed propagation is not multi-threaded, so we thread on seeds instead...
}

//----------------------------------------------------------------------------------------

vtkSmartPointer< vtkPolyData > ParticleFilteringTractographyAlgorithm::PropagateSeed(Self::PhysicalPoint point)
{
    // Diffusion directions provided by the model at point
    std::vector< btk::GradientDirection > nextDirections = m_DiffusionModel->MeanDirectionsAt(point);

    // Graphical representation structures
    vtkSmartPointer< vtkPoints >        vpoints = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray >       lines = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkPolyData > currentFiber = vtkSmartPointer< vtkPolyData >::New();

    // Processing
    for(std::vector< btk::GradientDirection >::iterator it = nextDirections.begin(); it != nextDirections.end(); it++)
    {
        //
        // Estimate Fiber
        //

        std::vector< Self::PhysicalPoint > points;
        points.push_back(point);

        btk::GradientDirection nextDirection = *it;

        this->PropagateSeed(points, nextDirection);

        //
        // Build graphical fiber
        //

        if(points.size() > 1)
        {
            vtkSmartPointer< vtkPolyLine > line = vtkSmartPointer< vtkPolyLine >::New();

            line->GetPointIds()->SetNumberOfIds(points.size());

            unsigned int nbOfPreviousPoints = vpoints->GetNumberOfPoints();

            vpoints->InsertNextPoint(-points[0][0], -points[0][1], points[0][2]);
            line->GetPointIds()->SetId(0,nbOfPreviousPoints);

            for(unsigned int i = 1; i < points.size(); i++)
            {
                // Insert cells (points and lines)
                vpoints->InsertNextPoint(-points[i][0], -points[i][1], points[i][2]);
                line->GetPointIds()->SetId(i,i+nbOfPreviousPoints);
            }

            lines->InsertNextCell(line);
        }
    } // for each direction

    // Build the whole trajectories
    currentFiber->SetPoints(vpoints);
    currentFiber->SetLines(lines);

//    std::stringstream filename;
//    filename << "pmap-" << point << "-" << ".nii.gz";
//    btk::ImageHelper< ProbabilityMap >::WriteImage(m_ProbabilityMap, filename.str());

    return currentFiber;
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::PropagateSeed(std::vector< Self::PhysicalPoint > &points, GradientDirection nextDirection)
{
//    btkTicTocInit();
//    btkTic();

    // Initialize particle's cloud
    std::vector< Particle > cloud(m_NumberOfParticles, Particle(points.back()));

    unsigned int numberOfActiveParticles = m_NumberOfParticles;
    unsigned int      numberOfIterations = 0;

    // Get the first point
    PhysicalPoint x0 = points.back();

    // Initial sampling
//    #pragma omp parallel for default(shared) schedule(dynamic)
    for(unsigned int m = 0; m < m_NumberOfParticles; m++)
    {
        // Simulate a direction using initial density
        GradientDirection v0 = m_ImportanceDensity.Simulate(nextDirection, m_ImportanceDensity.EstimateConcentrationParameter(nextDirection, x0)) * m_ParticleStepSize;

        // Move particle (inside the mask)
        Self::PhysicalPoint x1 = x0 + v0;

        // Check if the physical point is in the mask
        Self::MaskImage::IndexType maskIndex;

        if(m_Mask->TransformPhysicalPointToIndex(x1, maskIndex) && m_Mask->GetPixel(maskIndex) > 0)
        {
            cloud[m].AddToPath(v0, x1, m_InitialWeight);
        }
        else // x1 is not inside the mask
        {
            numberOfActiveParticles--;
            cloud[m].Desactivate();
        }
    } // for m particles

    numberOfIterations++;


    // Sequential sampling
    while(numberOfActiveParticles > 0)
    {
        // Move particle and update weight
//        #pragma omp parallel for default(shared) schedule(dynamic)
        for(unsigned int m = 0; m < m_NumberOfParticles; m++)
        {
            if(cloud[m].IsActive())
            {
                // State of current particle
                GradientDirection   vkm1 = cloud[m].GetLastDirection();
                Self::PhysicalPoint   xk = cloud[m].GetLastPoint();

                // Simulate next direction
                GradientDirection mu = m_ImportanceDensity.GetMeanDirection(xk, vkm1);

                double         kappa;
                GradientDirection vk;

                if(mu.IsNull())
                {
                    kappa = m_PriorDensity.GetConcentration();
                    vk    = m_ImportanceDensity.Simulate(vkm1, kappa) * m_ParticleStepSize;
                }
                else // mean is not null
                {
                    kappa = m_ImportanceDensity.EstimateConcentrationParameter(mu, xk);
                    vk    = m_ImportanceDensity.Simulate(mu, kappa) * m_ParticleStepSize;
                }

                // Move particle and update weight (inside mask)
                Self::PhysicalPoint xkp1 = xk + vk;
                MaskImage::IndexType maskIndex;

                if(m_Mask->TransformPhysicalPointToIndex(xkp1, maskIndex) && m_Mask->GetPixel(maskIndex) > 0)
                {
                    double weight = cloud[m].GetLastWeight();
                    double likelihood = 1.0;

                    if(mu.IsNull())
                    {
                        double      prior = m_PriorDensity.Evaluate(vk, vkm1);
                        double importance = m_ImportanceDensity.Evaluate(vk, vkm1, kappa);

                        weight *= prior / importance;
                    }
                    else // mean is not null
                    {
                               likelihood = m_LikelihoodDensity.Evaluate(vk, xk, mu);
                        double      prior = m_PriorDensity.Evaluate(vk, vkm1);
                        double importance = m_ImportanceDensity.Evaluate(vk, mu, kappa);

                        weight *= std::exp( likelihood + prior - importance );
                    }

                    cloud[m].AddToPath(vk, xkp1, weight);
                    cloud[m].AddLikelihood(likelihood);
                }
                else // xkp1 is not inside the mask
                {
                    numberOfActiveParticles--;
                    cloud[m].Desactivate();
                }
            }
        } // for m particles

        // Normalize weights (compute the sum and normalize) and compute sum squared
        double         sum = 0.0;
        double sumSquarred = 0.0;

//        #pragma omp parallel for default(shared) schedule(dynamic)
        for(unsigned int m = 0; m < m_NumberOfParticles; m++)
        {
            if(cloud[m].IsActive())
            {
                double weight = cloud[m].GetLastWeight();
                sum         += weight;
                sumSquarred += weight * weight;
            }
        } // for m particles

//        #pragma omp parallel for default(shared) schedule(dynamic)
        for(unsigned int m = 0; m < m_NumberOfParticles; m++)
        {
            if(cloud[m].IsActive())
            {
                cloud[m].NormalizeLastWeightWith(sum);
            }
        }

        // Compute resampling threshold
        double ESS = 1.0 / sumSquarred;

        // Resample if ESS is below a threshold
        // keeping proportionnality of weights
        if(ESS < m_ResamplingThreshold)
        {
            this->ResampleParticlesCloud(cloud, numberOfActiveParticles);
        }

        numberOfIterations++;
    }

//    btkCoutMacro("Fin propagation.");
//    btkToc();

//    this->SaveCloud(cloud);

    // Build output
//    this->ComputeProbabilityMap(cloud);
//    btkTic();
    this->ComputeMaximumAPosteriori(points, cloud, numberOfIterations, nextDirection);
//    btkCoutMacro("Fin MAP.");
//    btkToc();
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::ResampleParticlesCloud(std::vector< Particle > &cloud, unsigned numberOfActiveParticles)
{
    // Initialize
    double cumulative = 0.0;
    std::vector< double > intervals(numberOfActiveParticles, 0.0);
    std::vector< unsigned int > indices(numberOfActiveParticles, 0);

    // Create proportionnal intervals between 0 and 1 for inside particles
    unsigned int count = 0;
    for(unsigned int m = 0; m < m_NumberOfParticles; m++)
    {
        if(cloud[m].IsActive())
        {
            cumulative      += cloud[m].GetLastWeight();
            intervals[count] = cumulative;
            indices[count]   = m;
            count++;
        }
    } // for m particles

    double initialWeight = 1.0 / numberOfActiveParticles;

    // Get particles from the cloud keeping proportionnality (multinomial resampling)
//    #pragma omp parallel for default(shared) schedule(dynamic)
    for(unsigned int m = 0; m < m_NumberOfParticles; m++)
    {
        if(cloud[m].IsActive())
        {
            // Simulate x ~ U(0,1)
            double x = static_cast< double >(rand()) / static_cast< double >(RAND_MAX);

            bool found = false;
            unsigned int i = 0;

            do
            {
                if(x < intervals[i])
                {
                    found = true;
                }
                else
                {
                    i++;
                }
            } while(!found && i < numberOfActiveParticles-1);

            cloud[m].Resample(cloud[indices[i]].GetLastPoint(), initialWeight);
        }
    } // for m particles
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::SaveCloud(std::vector< Particle > &cloud)
{
    //
    // Build VTK PolyData from particle's cloud
    //

    vtkSmartPointer< vtkPoints >       points = vtkSmartPointer< vtkPoints >::New();
    vtkSmartPointer< vtkCellArray >     lines = vtkSmartPointer< vtkCellArray >::New();
    vtkSmartPointer< vtkDoubleArray > weights = vtkSmartPointer< vtkDoubleArray >::New();

    weights->SetNumberOfComponents(1);

    // Create all points and lines
    for(unsigned int m = 0; m < cloud.size(); m++)
    {
        unsigned int pathLength = cloud[m].GetPathLength();

        if(pathLength > 1)
        {
            // Initial point
            Self::PhysicalPoint point = cloud[m].GetPointAtStep(0);
            points->InsertNextPoint(point[0], point[1], point[2]);

            double weight = cloud[m].GetWeightAtStep(0);
            weights->InsertNextTupleValue(&weight);

            // Next points
            for(unsigned int i = 1; i < pathLength-1; i++)
            {
                Self::PhysicalPoint estimatedPoint = cloud[m].GetPointAtStep(i-1) + cloud[m].GetVectorAtStep(i-1);
                point = cloud[m].GetPointAtStep(i);

                vtkIdType pid = points->InsertNextPoint(point[0], point[1], point[2]);
                weight = cloud[m].GetWeightAtStep(i);
                weights->InsertNextTupleValue(&weight);

                // If point is linked with previous, add a line between them
                if(point == estimatedPoint)
                {
                    vtkSmartPointer< vtkLine > line = vtkSmartPointer< vtkLine >::New();
                    line->GetPointIds()->SetId(0, pid-1);
                    line->GetPointIds()->SetId(1, pid);
                    lines->InsertNextCell(line);
                }
            } // for i

            // Last point
            Self::PhysicalPoint estimatedPoint = cloud[m].GetPointAtStep(pathLength-2) + cloud[m].GetVectorAtStep(pathLength-2);
            point = cloud[m].GetPointAtStep(pathLength-1);

            vtkIdType pid = points->InsertNextPoint(point[0], point[1], point[2]);
            weight = cloud[m].GetWeightAtStep(pathLength-2);
            weights->InsertNextTupleValue(&weight);

            // If point is linked with previous, add a line between them
            if(point == estimatedPoint)
            {
                vtkSmartPointer< vtkLine > line = vtkSmartPointer< vtkLine >::New();
                line->GetPointIds()->SetId(0, pid-1);
                line->GetPointIds()->SetId(1, pid);
                lines->InsertNextCell(line);
            }
        }
    } // for each particle

    // create vtk polydata object
    vtkSmartPointer< vtkPolyData > polydata = vtkSmartPointer< vtkPolyData >::New();
    polydata->SetPoints(points);
    polydata->SetLines(lines);
    polydata->GetPointData()->SetScalars(weights);


    //
    // Save VTK PolyData object into VTK file
    //

    std::stringstream filename;
    filename << "cloud-" << cloud[0].GetPointAtStep(0) << "-" << cloud[0].GetPathLength() << ".vtk";

    vtkSmartPointer< vtkPolyDataWriter > writer = vtkSmartPointer< vtkPolyDataWriter >::New();
    writer->SetInput(polydata);
    writer->SetFileName(filename.str().c_str());
    writer->Write();
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::ComputeProbabilityMap(std::vector< Particle > &cloud)
{
    ProbabilityMap::RegionType region = m_ProbabilityMap->GetLargestPossibleRegion();

    for(unsigned int m = 0; m < cloud.size(); m++)
    {
        bool stop = false;

        for(unsigned int k = 0; k < cloud[m].GetPathLength()-1 && !stop; k++)
        {
            PhysicalPoint pk = cloud[m].GetPointAtStep(k);
            double    weight = cloud[m].GetWeightAtStep(k);

            // Tri-linear interpolation
            DiffusionSignal::ContinuousIndex cindex;
            if(m_ProbabilityMap->TransformPhysicalPointToContinuousIndex(pk, cindex))
            {
                // Get current index coordinates
                unsigned short x = static_cast< unsigned short >(std::floor(cindex[0]));
                unsigned short y = static_cast< unsigned short >(std::floor(cindex[1]));
                unsigned short z = static_cast< unsigned short >(std::floor(cindex[2]));

                // Compute distances
                double dx = cindex[0] - x;
                double dy = cindex[1] - y;
                double dz = cindex[2] - z;

                // Set interpolated values
                ProbabilityMap::IndexType index;

                index[0] = x; index[1] = y; index[2] = z;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + (1-dz)*(1-dx)*(1-dy) * weight);
                }

                index[0] = x; index[1] = y+1; index[2] = z;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + (1-dz)*(1-dx)*dy * weight);
                }

                index[0] = x+1; index[1] = y; index[2] = z;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + (1-dz)*dx*(1-dy) * weight);
                }

                index[0] = x+1; index[1] = y+1; index[2] = z;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + (1-dz)*dx*dy * weight);
                }

                index[0] = x; index[1] = y; index[2] = z+1;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + dz*(1-dx)*(1-dy) * weight);
                }

                index[0] = x; index[1] = y+1; index[2] = z+1;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + dz*(1-dx)*dy * weight);
                }

                index[0] = x+1; index[1] = y; index[2] = z+1;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + dz*dx*(1-dy) * weight);
                }

                index[0] = x+1; index[1] = y+1; index[2] = z+1;
                if(region.IsInside(index))
                {
                    m_ProbabilityMap->SetPixel(index, m_ProbabilityMap->GetPixel(index) + dz*dx*dy * weight);
                }
            }
            else // pk is no longer in image
            {
                stop = true;
            }
        } // for each step of the particle's path
    } // for each particle
}

//----------------------------------------------------------------------------------------

void ParticleFilteringTractographyAlgorithm::ComputeMaximumAPosteriori(std::vector< Self::PhysicalPoint > &map, std::vector< Particle > &cloud, unsigned int numberOfIterations, btk::GradientDirection nextDirection)
{
    if(numberOfIterations > 2)
    {
        // Initialization
        unsigned int t = numberOfIterations-1;
//        std::vector< double > delta(t*numberOfIterations, 0.0);
        std::vector< double > delta(t*m_NumberOfParticles, 0.0);
//        std::vector< double >   psi((t-1)*numberOfIterations, 0.0);
        std::vector< double >   psi((t-1)*m_NumberOfParticles, 0.0);

        //
        // Compute delta and psi using dynamic programming
        //

        for(unsigned int m = 0; m < m_NumberOfParticles; m++)
        {
            if(cloud[m].GetPathLength()-1 > 1)
            {
                assert(INDEX(0,m) < delta.size());
                delta[INDEX(0,m)] = std::log(m_PriorDensity.Evaluate(cloud[m].GetVectorAtStep(1), cloud[m].GetVectorAtStep(0))) + cloud[m].GetLikelihoodAtStep(0);
            }
            else // cloud[m].GetPathLength-1 <= 1
            {
                assert(INDEX(0,m) < delta.size());
                delta[INDEX(0,m)] = std::log(0.1 / m_NumberOfParticles);
            }
        } // for each particle

        for(unsigned int k = 1; k < t; k++)
        {
            for(unsigned int m = 0; m < m_NumberOfParticles; m++)
            {
                if(k+1 < cloud[m].GetPathLength()-1)
                {
                    double max        = DBL_MIN;
                    unsigned int imax = 0;

                    PhysicalPoint xkp1_m = cloud[m].GetPointAtStep(k+1);

                    for(unsigned int i = 0; i < m_NumberOfParticles; i++)
                    {
                        if(k < cloud[i].GetPathLength()-1)
                        {
                            PhysicalPoint actual_xkp1_i = cloud[i].GetPointAtStep(k+1);
                            PhysicalPoint should_xkp1_i = cloud[i].GetPointAtStep(k) + cloud[i].GetVectorAtStep(k);
                            double    distance_particle = std::sqrt( (actual_xkp1_i[0]-xkp1_m[0])*(actual_xkp1_i[0]-xkp1_m[0]) + (actual_xkp1_i[1]-xkp1_m[1])*(actual_xkp1_i[1]-xkp1_m[1]) + (actual_xkp1_i[2]-xkp1_m[2])*(actual_xkp1_i[2]-xkp1_m[2]));
                            double      distance_points = std::sqrt( (should_xkp1_i[0]-actual_xkp1_i[0])*(should_xkp1_i[0]-actual_xkp1_i[0]) + (should_xkp1_i[1]-actual_xkp1_i[1])*(should_xkp1_i[1]-actual_xkp1_i[1]) + (should_xkp1_i[2]-actual_xkp1_i[2])*(should_xkp1_i[2]-actual_xkp1_i[2]));

                            if(distance_points < 1.0 && distance_particle < 1.0)
                            {
                                assert(INDEX(k-1,i) < delta.size());
                                double tmp = delta[INDEX(k-1,i)] + std::log(m_PriorDensity.Evaluate(cloud[m].GetVectorAtStep(k+1), cloud[i].GetVectorAtStep(k)));

                                if(tmp > max)
                                {
                                    max  = tmp;
                                    imax = i;
                                }
                            }
                        }
                    } // for each other particle

                    if(btkFloatingEqual(max, DBL_MIN)) // no maximum found
                    {
                        assert(INDEX(k,m) < delta.size());
                        delta[INDEX(k,m)] = cloud[m].GetLikelihoodAtStep(k);
                        psi[INDEX(k-1,m)] = m;
                    }
                    else // max != MIN_REAL
                    {
                        assert(INDEX(k,m) < delta.size());
                        delta[INDEX(k,m)] = cloud[m].GetLikelihoodAtStep(k) + max;
                        psi[INDEX(k-1,m)] = imax;
                    }
                }
                else // k+1 >= cloud[m].GetPathLength()-1
                {
                    if(k == 1)
                    {
                        assert(INDEX(1,m) < delta.size());
                        assert(INDEX(0,m) < delta.size());
                        delta[INDEX(1,m)] = delta[INDEX(0,m)];
                    }
                    else // k > 1
                    {
                        assert(INDEX(k,m) < delta.size());
                        assert(INDEX(k-1,m) < delta.size());
                        assert(INDEX(k-2,m) < delta.size());
                        delta[INDEX(k,m)] = delta[INDEX(k-1,m)] + std::abs(delta[INDEX(k-1,m)] - delta[INDEX(k-2,m)]);
                    }

                    psi[INDEX(k-1,m)] = m;
                }
            } // for each particle
        } // for each step


        //
        // Backtracking
        //

        std::vector< PhysicalPoint > backwardPoints;

        double max        = DBL_MIN;
        unsigned int mmax = 0;

        for(unsigned int m = 0; m < m_NumberOfParticles; m++)
        {
            assert(INDEX(t-1,m) < delta.size());
            double tmp = delta[INDEX(t-1,m)];

            if(max < tmp)
            {
                max  = tmp;
                mmax = m;
            }
        } // for each particle

        int len = cloud[mmax].GetPathLength()-1;
//        PhysicalPoint backwardPoint = cloud[mmax].GetPointAtStep(len);
//        backwardPoints.push_back(backwardPoint);
        backwardPoints.push_back(cloud[mmax].GetPointAtStep(len));

        for(int k = t-2; k >= 0; k--)
        {
            mmax = psi[INDEX(k,mmax)];

            int len = cloud[mmax].GetPathLength()-1;

            if(k+1 <= len)
            {
                backwardPoints.push_back(cloud[mmax].GetPointAtStep(k+1));
            }
            else // k+1 > len
            {
                backwardPoints.push_back(cloud[mmax].GetPointAtStep(len));
            }
        } // for each step


        //
        // Get forward MAP estimate
        //

        for(int k = backwardPoints.size()-1; k >= 0; k--)
        {
            map.push_back(backwardPoints[k]);
        }
    }
    else if(numberOfIterations > 1)
    {
//        double maxWeight = cloud[0].GetLastWeight();
        double maxWeight = DBL_MIN;
        unsigned int   i = 0;

        for(unsigned int m = 0; m < cloud.size(); m++)
        {
            if(cloud[m].GetPathLength() > 1)
            {
//                double weight = cloud[m].GetLastWeight();
                double weight = m_PriorDensity.Evaluate(cloud[m].GetLastDirection(), nextDirection);

                if(weight > maxWeight)
                {
                    weight = maxWeight;
                    i      = m;
                }
            }
        } // for each particle

        if(cloud[i].GetPathLength() > 1)
        {
            map.push_back(cloud[i].GetLastPoint());
        }
    }
//    else // numberOfIterations <= 1
//    {
//        // Do nothing
//    }


    //
    // Use mobile mean for smoothing
    //

    for(unsigned int step = 0; step < 10; step++)
    {
        for(unsigned int i = 2; i < map.size(); i++)
        {
            PhysicalPoint p;
            p[0] = (map[i-2][0] + map[i-1][0] + map[i][0]) / 3.0;
            p[1] = (map[i-2][1] + map[i-1][1] + map[i][1]) / 3.0;
            p[2] = (map[i-2][2] + map[i-1][2] + map[i][2]) / 3.0;

            map[i-1] = p;
        }
    }
}

} // namespace btk
