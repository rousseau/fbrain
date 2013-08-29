/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date:
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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


#include "btkWeightedEstimationFilter.h"
#include "iomanip"


static itk::SimpleMutexLock Mutex;
namespace btk
{

WeightedEstimationFilter::WeightedEstimationFilter() : Superclass()
{
    m_NbLoop  = 0;
    m_Percent = 0.0;
    m_Radius  = 0;

}

//----------------------------------------------------------------------------------------

WeightedEstimationFilter::~WeightedEstimationFilter()
{
    // ----
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    // ----
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::SetDiffusionDataset(DatasetPointer dataset)
{
    m_Dataset = dataset;
    Self::SetNumberOfRequiredInputs(0);
    //    Self::SetNumberOfThreads(1);
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::AllocateOutputs()
{
    ////////////////////////////////////////////////////////////////////////////
    //
    // Allocate memory for output image
    //
    DiffusionSequence::RegionType    sequenceRegion     = m_Dataset -> GetSequenceRegion();
    DiffusionSequence::SpacingType   sequenceSpacing    = m_Dataset -> GetSpacing();
    DiffusionSequence::PointType     sequenceOrigin     = m_Dataset -> GetOrigin();
    DiffusionSequence::DirectionType sequenceDirection  = m_Dataset -> GetDirection();

    m_Size4D = sequenceRegion.GetSize();

    m_OutputSequence = Self::GetOutput();

    m_OutputSequence->SetRegions(sequenceRegion);
    m_OutputSequence->SetSpacing(sequenceSpacing);
    m_OutputSequence->SetOrigin(sequenceOrigin);
    m_OutputSequence->SetDirection(sequenceDirection);
    m_OutputSequence->Allocate();
    m_OutputSequence->FillBuffer(0);

}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::BeforeThreadedGenerateData()
{
    // ----
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{

    ////////////////////////////////////////////////////////////////////////////
    // Set radius
    ////////////////////////////////////////////////////////////////////////////
    NeighborhoodIterator::RadiusType radius;
    radius[0]= m_Radius * m_Dataset->GetSpacing()[0];
    radius[1]= m_Radius * m_Dataset->GetSpacing()[1];
    radius[2]= 0;

    itk::ImageRegionIterator< Self::OutputImageType > outIt(m_OutputSequence, outputRegionForThread);
    for(outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
    {
        VectorType              distanceVector;
        VectorType              signalVector;
        GradientTableType       directionTable;

        DiffusionSequence::IndexType queryIndex4D = outIt.GetIndex();
        DiffusionSequence::PointType queryPoint4D;
        m_OutputSequence->TransformIndexToPhysicalPoint(queryIndex4D,queryPoint4D);

        Point3DType queryPoint3D;
        queryPoint3D[0] = queryPoint4D[0];
        queryPoint3D[1] = queryPoint4D[1];
        queryPoint3D[2] = queryPoint4D[2];

        for (DataIterator dataIt = m_Dataset ->begin() ; dataIt != m_Dataset -> end(); ++dataIt)
        {
            Index3DType index3D;
            (*dataIt)->TransformPhysicalPointToIndex(queryPoint3D,index3D);

            ////////////////////////////////////////////////////////////////////////////
            // Consider only slices which have a distance <= Radius form the queryIndex
            ////////////////////////////////////////////////////////////////////////////
            if(std::abs(index3D[2] ) <= std::ceil(m_Radius) &&  (*dataIt)->GetLargestPossibleRegion().IsInside(index3D))
            {
                index3D[2]=0;

                ////////////////////////////////////////////////////////////////////////////
                // Search the neighobr in the slice according to the radius
                ////////////////////////////////////////////////////////////////////////////
                NeighborhoodIterator NeighborIt(radius, *dataIt,(*dataIt)->GetRequestedRegion());
                NeighborIt.SetLocation(index3D);
                for (unsigned int i = 0; i < NeighborIt.Size(); ++i)
                {
                    ////////////////////////////////////////////////////////////////////////////
                    // Get value and distance of the neighbor
                    ////////////////////////////////////////////////////////////////////////////
                    float voxelValue = static_cast<float>(NeighborIt.GetPixel(i));

                    Point3DType neighborPoint;
                    (*dataIt)->TransformIndexToPhysicalPoint(index3D,neighborPoint);
                    float distance = neighborPoint.EuclideanDistanceTo(queryPoint3D);

                    ////////////////////////////////////////////////////////////////////////////
                    // Add neighbor parameters in containers
                    ////////////////////////////////////////////////////////////////////////////
                    distanceVector.push_back(distance);
                    signalVector.push_back(voxelValue);
                    directionTable.push_back((*dataIt)->GetGradientDirection());

                } // end for NeighborIt
            } // end slice selection NeighborIt
        } // end for through dataset

        if(distanceVector.size() !=0)
        {
            ////////////////////////////////////////////////////////////////////////////
            // Spherical harmonic decomposition and Model estimation
            ////////////////////////////////////////////////////////////////////////////
            WeightedEstimationPointer Estimation = WeightedEstimationType::New();

            // Weighted Spherical harmonic decomposition part
            Estimation -> SetNeighborsDistances(distanceVector);
            Estimation -> SetSignalValues(signalVector);
            Estimation -> SetGradientTable(directionTable);
            Estimation -> SetSphericalHarmonicsOrder(4);
            Estimation -> SetSphericalResolution(300);
            Estimation -> Initialize();
            Estimation -> Update();

            outIt.Set(Estimation->SignalAt(m_Dataset->GetGradientTable()[queryIndex4D[3]]));
        }
        else
        {
            outIt.Set(0);
        }

        Mutex.Lock();
        m_NbLoop+=1;
        m_Percent = double(m_NbLoop)*100.0/((m_Size4D[0]*m_Size4D[1]*m_Size4D[2]*m_Size4D[3]));
        std::cout<<"\r  -> Estimation ... "<<std::fixed<<std::setprecision(1)<<m_Percent<<" %";
        std::cout<<std::fixed<<std::setprecision(5);
        Mutex.Unlock();
    }
}



} // end namespace btk
