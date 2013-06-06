/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 04/06/2013
  Author(s):Frederic Champ (champ(at)unistra.fr)
  
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

#include <iomanip>

// Local inclues
#include "btkSphericalHarmonics.h"
#include "btkImageHelper.h"
#include "itkSimpleFastMutexLock.h"

static itk::SimpleMutexLock Mutex;
namespace btk
{

WeightedEstimationFilter::WeightedEstimationFilter() : Superclass()
{
    m_NbLoop = 0;
    m_Percent = 0.0;
    m_Sigma = 1.0;
    m_EllipsoidSize = 1.0;
}

//----------------------------------------------------------------------------------------

WeightedEstimationFilter::~WeightedEstimationFilter()
{
    // ----
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::SetInputSequence(SequenceConstPointer InputSequence)
{
    m_InputSequence = btk::ImageHelper<TSequence>::DeepCopy(InputSequence);
    m_GradientTable = InputSequence -> GetGradientTable();
    m_Size4D = m_InputSequence -> GetLargestPossibleRegion().GetSize();

    Self::SetNumberOfRequiredInputs(0);

    //Self::SetNumberOfThreads(1);
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::AllocateOutputs()
{
    ////////////////////////////////////////////////////////////////////////////
    //
    // Allocate memory for output image
    //
    DiffusionSequence::RegionType           sequenceSize = m_InputSequence->GetRequestedRegion();
    DiffusionSequence::SpacingType     sequenceSpacing = m_InputSequence->GetSpacing();
    DiffusionSequence::PointType        sequenceOrigin = m_InputSequence->GetOrigin();
    DiffusionSequence::DirectionType sequenceDirection = m_InputSequence->GetDirection();

    Self::OutputImageType::Pointer output = Self::GetOutput();
    output->SetRegions(sequenceSize);
    output->SetSpacing(sequenceSpacing);
    output->SetOrigin(sequenceOrigin);
    output->SetDirection(sequenceDirection);
    output->Allocate();
    output->FillBuffer(0);
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::BeforeThreadedGenerateData()
{
    ////////////////////////////////////////////////////////////////////////////
    //
    // Generate a gradient image
    //
    GradientFilterPointer gradientFilter = GradientFilterType::New();
    gradientFilter -> SetInput(m_InputSequence );
    gradientFilter -> SetUseImageSpacing(true);
    gradientFilter -> SetUseImageDirection(true);
    gradientFilter -> Update();

    m_GradientImage = gradientFilter -> GetOutput();
}

//----------------------------------------------------------------------------------------

void WeightedEstimationFilter::ThreadedGenerateData(const OutputImageRegionType &outputRegionForThread, itk::ThreadIdType threadId)
{

    float gradX      = 1.0;
    float gradY      = 1.0;
    float gradZ      = 1.0;
    float norm       = 0.0;
    int sizeSearch   = m_EllipsoidSize +1 ; // Search windows arount the query point

    GradientIterator It(m_GradientImage,outputRegionForThread);
    itk::ImageRegionIterator< Self::OutputImageType > outIt(this->GetOutput(), outputRegionForThread);

    for(It.GoToBegin(), outIt.GoToBegin(); !It.IsAtEnd(), !outIt.IsAtEnd(); ++It, ++outIt)
    {


        TGradient::IndexType    gradIndex = It.GetIndex();
        VectorType        IndexVector;
        VectorType        SignalVector;
        GradientTableType       DirectionTable;

        ////////////////////////////////////////////////////////////////////////////
        //
        // Get gradients values thru the 3 directions
        // and normalize them
        //
        float valX= m_GradientImage->GetPixel(gradIndex)[0];
        float valY= m_GradientImage->GetPixel(gradIndex)[1];
        float valZ= m_GradientImage->GetPixel(gradIndex)[2];

        if(valX>0.0)
        {
            gradX=valX;
        }
        if(valY>0.0)
        {
            gradY= valY;
        }
        if(valZ>0)
        {
            gradZ= valZ;
        }

        norm = std::sqrt(gradX*gradX+gradY*gradY+gradZ*gradZ);

        ////////////////////////////////////////////////////////////////////////////
        //
        // Calculate semi-principal axts of the ellipsoid as the inverse of
        // normalized gradient values multiplied by the "m_EllipsoidSize" factor.
        //
        if(norm !=0.0)
        {
            gradX = (1-std::abs(gradX/norm))*m_EllipsoidSize;
            gradY = (1-std::abs(gradY/norm))*m_EllipsoidSize;
            gradZ = (1-std::abs(gradZ/norm))*m_EllipsoidSize;
        }
        else // Find a good way
        {
            gradX = m_EllipsoidSize;
            gradY = m_EllipsoidSize;
            gradZ = m_EllipsoidSize;
        }

        ////////////////////////////////////////////////////////////////////////////
        //
        // Ellipsoid set up.
        //
        EllipsoidFunctionPointer EllipsoidFunction = EllipsoidFunctionType::New();

        EllipsoidFunctionVector center;
        center[0]=gradIndex[0];
        center[1]=gradIndex[1];
        center[2]=gradIndex[2];
        EllipsoidFunction -> SetCenter(center);

        ////////////////////////////////////////////////////////////////////////////
        //
        // To have principal axes gradX, gradY and gradZ are multiplied by 2.0
        //
        EllipsoidFunctionVector axesLength;
        axesLength[0] =  gradX*2.0;
        axesLength[1] =  gradY*2.0;
        axesLength[2] =  gradZ*2.0;
        EllipsoidFunction -> SetAxes(axesLength);

        ////////////////////////////////////////////////////////////////////////////
        //
        // Set orientation of ellipsoid axes
        //
        double data[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        vnl_matrix<double> orientations (data, 3, 3);

        EllipsoidFunction -> SetOrientations(orientations);

        ////////////////////////////////////////////////////////////////////////////
        //
        // Get all neighbors inside the ellipsoid and who are are not outliers
        //
        for(int x=-sizeSearch; x<=sizeSearch; x++)
        {
            for(int y=-sizeSearch; y<=sizeSearch; y++)
            {
                for(int z=-sizeSearch; z<=sizeSearch; z++)
                {
                    // Get neighbors in the search window
                    double NeighborIndexes[3];
                    NeighborIndexes[0]=gradIndex[0]+x;
                    NeighborIndexes[1]=gradIndex[1]+y;
                    NeighborIndexes[2]=gradIndex[2]+z;

                    // check if the index is not out of the image
                    if( NeighborIndexes[0] >= 0 && NeighborIndexes[0] < m_Size4D[0] &&
                            NeighborIndexes[1] >= 0 && NeighborIndexes[1] < m_Size4D[1] &&
                            NeighborIndexes[2] >= 0 && NeighborIndexes[2] < m_Size4D[2] )
                    {
                        // check if the neighbor is inside the ellipsoid
                        bool IsInsideEllipsoid = EllipsoidFunction -> Evaluate(NeighborIndexes);

                        if(IsInsideEllipsoid)
                        {
                            unsigned int SliceIndex = NeighborIndexes[2];
                            for(unsigned int i = 0; i<m_Size4D[3]; i++)
                            {

                                ////////////////////////////////////////////////////////////////////////////
                                //
                                // If the neighobrs is not an outliers, his distance from the query point,
                                // his value and his gradient direction are added to vectors who will be use
                                // to compute the weighted spherical harmonics decomposition
                                //
                                if(!m_BoolIndexes[SliceIndex][i])
                                {

                                    GradPoint QueryPoint;
                                    GradPoint NeighborPoint;
                                    SequenceIndexType NeighborSequenceIndexes;
                                    NeighborSequenceIndexes[0] = NeighborIndexes[0];
                                    NeighborSequenceIndexes[1] = NeighborIndexes[1];
                                    NeighborSequenceIndexes[2] = NeighborIndexes[2];
                                    NeighborSequenceIndexes[3] = i;

                                    float voxelValue = static_cast<float>(m_InputSequence -> GetPixel(NeighborSequenceIndexes));

                                    m_GradientImage -> TransformIndexToPhysicalPoint(gradIndex,QueryPoint);
                                    m_GradientImage -> TransformIndexToPhysicalPoint(NeighborSequenceIndexes,NeighborPoint) ;

                                    float norm = 0.0;
                                    for (unsigned int k=0;k<3;k++)
                                    {
                                        norm += (QueryPoint[k]-NeighborPoint[k])*(QueryPoint[k]-NeighborPoint[k]);
                                    }
                                    norm=std::sqrt(norm);

                                    IndexVector.push_back(norm);
                                    SignalVector.push_back(voxelValue);
                                    DirectionTable.push_back(m_GradientTable[i]);

                                } // end if is not outliers
                            } // end for all volumes
                        } // end if is inside ellipsoid
                    } // end if is inside buffered region
                } // end z search
            } // end y search
        } // end x search


        if(m_VerboseMod)
        {
            btkCoutMacro("  QueryIndex : "<<gradIndex<<" is outliers ? : "<<m_BoolIndexes[gradIndex[2]][gradIndex[3]]);
            btkCoutMacro("  Number of neighbors non outliers : "<<IndexVector.size())
        }

        if(IndexVector.size() !=0)
        {
            ////////////////////////////////////////////////////////////////////////////
            //
            // Spherical harmonic decomposition and Model estimation
            //
            WeightedEstimationPointer Estimation = WeightedEstimationType::New();

            // Weighted Spherical harmonic decomposition part
            Estimation -> SetNeighborsDistances(IndexVector);
            Estimation -> SetSignalValues(SignalVector);
            Estimation -> SetGradientTable(DirectionTable);
            Estimation -> SetSphericalHarmonicsOrder(4);

            if(!m_BoolIndexes[gradIndex[2]][gradIndex[3]])
                Estimation -> SetQueryValue(m_InputSequence -> GetPixel(gradIndex));

            Estimation -> SetSigma(m_Sigma);
            Estimation -> SetSphericalResolution(300);
            Estimation -> Initialize();
            Estimation -> Update();

            outIt.Set(Estimation->SignalAt(m_GradientTable[gradIndex[3]]));
        }

        if(!m_VerboseMod)
        {
            Mutex.Lock();
            m_NbLoop+=1;
            m_Percent = double(m_NbLoop)*100.0/((m_Size4D[0]*m_Size4D[1]*m_Size4D[2]*m_Size4D[3]));
            std::cout<<"\r  -> "<<std::fixed<<std::setprecision(1)<<m_Percent<<" %";
            std::cout<<std::fixed<<std::setprecision(5);
            Mutex.Unlock();
        }
    }



}

//----------------------------------------------------------------------------------------

}// end namespace btk
