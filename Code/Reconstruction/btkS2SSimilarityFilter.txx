/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date:
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

#ifndef __S2SSimilarityFilter_txx__
#define __S2SSimilarityFilter_txx__

#include "btkS2SSimilarityFilter.h"


namespace btk
{

template <typename TImage >
S2SSimilarityFilter<TImage>::S2SSimilarityFilter()
{
    m_VerboseMod = false;
    m_Delimiter  = ";";
    m_Method     = "Normalized";
    m_SliceNumber = 0;
    m_VolumeNumber = 0;
    m_NumberOfSlices = 0;
    m_NumberOfVolumes = 0;
}

//----------------------------------------------------------------------------------------

template <typename TImage >
S2SSimilarityFilter<TImage>::~S2SSimilarityFilter()
{
    //---
}


//----------------------------------------------------------------------------------------
template <typename TImage >
void
S2SSimilarityFilter<TImage>::SetSliceNumber(unsigned int SliceNumber)
{
    m_SliceNumber = SliceNumber;
    m_NumberOfSlices = 1;

}

//----------------------------------------------------------------------------------------

template <typename TImage >
void
S2SSimilarityFilter<TImage>::SetVolumeNumber(unsigned int VolumeNumber)
{
    if(VolumeNumber != 0)
    {
        m_VolumeNumber = VolumeNumber;
        m_NumberOfVolumes = 1;
    }
}

//----------------------------------------------------------------------------------------

template < typename TImage >
void S2SSimilarityFilter< TImage >::WriteData(std::string FileName)
{

    std::cout<<"  -> Writing similarities in "<<FileName<<"... ";

    std::ofstream fichierout(FileName.c_str(), std::ios::trunc);
    if(fichierout)
    {
        fichierout<<"slice"<<m_Delimiter;
        for (unsigned short i= m_VolumeNumber; i < m_VolumeNumber+m_NumberOfVolumes; i++)
        {
            fichierout<<"dwi"<<i<<m_Delimiter;
        }
        fichierout<<"mean"<<m_Delimiter<<"std"<<std::endl;

        for (unsigned short j=0; j < m_NumberOfSlices; j++)
        {
            fichierout <<j+m_SliceNumber<<m_Delimiter;
            for (unsigned short i=0; i < m_NumberOfVolumes; i++)
            {
                fichierout <<m_Similarity[j][i]<<m_Delimiter;
            }
            fichierout <<m_Mean[j]<<m_Delimiter<<m_StandardDeviation[j]<<std::endl;
        }
        fichierout.close();
        btkCoutMacro(" Done.");
    }
    else
        std::cerr << "Unreadable file !" << std::endl;
}
//----------------------------------------------------------------------------------------

template < typename TImage >
void S2SSimilarityFilter<TImage>::Initialize()
{
    SequenceRegionType region4D    = m_InputSequence -> GetLargestPossibleRegion();
    SequenceSizeType   input4Dsize = region4D.GetSize();

    m_Slice4DSize     = input4Dsize;
    m_Slice4DSize[2]  = 1;
    m_Slice4DSize[3]  = 0;

    if(m_NumberOfVolumes == 0)
    {
        m_VolumeNumber     = 1;
        m_NumberOfVolumes  = input4Dsize[3]-1;
    }
    if(m_NumberOfSlices == 0)
    {
        m_SliceNumber     = 0;
        m_NumberOfSlices  = input4Dsize[2];
    }

    m_Mean.resize(m_NumberOfSlices,0.0);
    m_StandardDeviation.resize(m_NumberOfSlices,0.0);
    m_Similarity.resize(m_NumberOfSlices);
}
//----------------------------------------------------------------------------------------
template < typename TImage >
void S2SSimilarityFilter<TImage>::Update()
{
    this->Initialize();
    for (unsigned int i = 0; i < m_NumberOfSlices; i++)
    {
        m_Similarity[i].resize(m_NumberOfVolumes);
    }
    unsigned int i,j;
    int nbloop=0;
    double percent;



    #pragma omp parallel for private(i,j) schedule(dynamic)
    for (j = m_SliceNumber; j <m_SliceNumber+m_NumberOfSlices; j++)
    {
        unsigned int VectorSliceIndex = j - m_SliceNumber;
        for (i = m_VolumeNumber; i < m_VolumeNumber + m_NumberOfVolumes; i++)
        {
            unsigned int VectorVolumeIndex =  i -m_VolumeNumber;

            SequenceIndexType  Slice4DIndex;
            Slice4DIndex[0]=0;
            Slice4DIndex[1]=0;
            Slice4DIndex[2] = j;
            Slice4DIndex[3] = i;

            unsigned int ImageIndex = i;
            unsigned int SliceIndex = j;


            ////////////////////////////////////////////////////////////////////////////
            //
            // Extraction of the reference slice
            //

            ImagePointer SliceReference;
            if(m_Reference)
            {
                ImageIndexType RefSliceIndex;
                RefSliceIndex[0] = 0;
                RefSliceIndex[1] = 0;
                RefSliceIndex[2] = SliceIndex;

                ImageSizeType  RefSliceSize;
                RefSliceSize[0] = m_Slice4DSize[0];
                RefSliceSize[1] = m_Slice4DSize[1];
                RefSliceSize[2] = m_Slice4DSize[2];

                ImageRegionType RefSliceRegion;
                RefSliceRegion.SetIndex( RefSliceIndex );
                RefSliceRegion.SetSize( RefSliceSize );


                ImageExtractorPointer ReferenceExtractor = ImageExtractor::New();
                ReferenceExtractor -> SetInput( m_Reference );
                ReferenceExtractor -> SetExtractionRegion( RefSliceRegion );
                ReferenceExtractor -> SetDirectionCollapseToSubmatrix( );
                ReferenceExtractor -> Update();

                SliceReference = ReferenceExtractor->GetOutput();
            }
            else
            {
                SequenceIndexType B0SliceIndex;
                B0SliceIndex[0] = 0;
                B0SliceIndex[1] = 0;
                B0SliceIndex[2] = SliceIndex;
                B0SliceIndex[3] = 0;

                SequenceRegionType B0SliceRegion;
                B0SliceRegion.SetIndex( B0SliceIndex );
                B0SliceRegion.SetSize( m_Slice4DSize );

                SequenceExtractorPointer B0Extractor = SequenceExtractor::New();
                B0Extractor -> SetInput( m_InputSequence );
                B0Extractor -> SetExtractionRegion( B0SliceRegion );
                B0Extractor -> SetDirectionCollapseToSubmatrix( );
                B0Extractor -> Update();

                SliceReference = B0Extractor->GetOutput();
            }

            ////////////////////////////////////////////////////////////////////////////
            //
            // Extraction of the DWI's slice
            //
            SequenceRegionType    Slice4DRegion;
            Slice4DRegion.SetIndex( Slice4DIndex );
            Slice4DRegion.SetSize( m_Slice4DSize );


            SequenceExtractorPointer ImageExtractor = SequenceExtractor::New();
            ImageExtractor -> SetExtractionRegion( Slice4DRegion );
            ImageExtractor -> SetInput( m_InputSequence );
            ImageExtractor -> SetDirectionCollapseToSubmatrix();
            ImageExtractor -> Update();

            ////////////////////////////////////////////////////////////////////////////
            //
            // Calculate Similarity between slices
            //
            MetricPointer Metric = MetricType::New();
            Metric -> SetFixedImage(SliceReference);
            Metric -> SetMovingImage( ImageExtractor->GetOutput() );
            Metric -> SetMethod(m_Method);
            Metric -> SetVerboseMod(m_VerboseMod);
            Metric -> Initialize();

            float value = std::fabs( Metric->GetValue( ) );
            m_Similarity[VectorSliceIndex][VectorVolumeIndex]= value;

            ////////////////////////////////////////////////////////////////////////////
            //
            // Display
            //
            if(m_VerboseMod)
                 #pragma omp critical
                std::cout<<"  Slice "<<SliceIndex<<" - Image "<<ImageIndex<<" - similarity = "<<value<<std::endl;
            else
            {
                percent = double(nbloop+1.0)*100.0/(m_NumberOfSlices*(m_NumberOfVolumes));
                #pragma omp critical
                std::cout<<"\r  -> Caculting similarities ... "<<std::fixed<<std::setprecision(1)<<percent<<" %";
                std::cout<<std::fixed<<std::setprecision(5);
            }

            ////////////////////////////////////////////////////////////////////////////

            nbloop+=1;
        }

        if(m_NumberOfVolumes > 1)
        {
            ////////////////////////////////////////////////////////////////////////////
            //
            // Calculate mean
            //
            for (unsigned short i=0; i < m_NumberOfVolumes; i++)
            {
                m_Mean[VectorSliceIndex]+=  m_Similarity[VectorSliceIndex][i]/(m_NumberOfVolumes);
            }

            if(m_VerboseMod)
                std::cout<<std::endl<<"    Slices "<<j<<" - mean = "<<m_Mean[VectorSliceIndex];

            ////////////////////////////////////////////////////////////////////////////
            //
            // Calculate standard deviation
            //
            double sum=0.0;
            for (unsigned short i=0; i < m_NumberOfVolumes; i++)
            {
                sum+=pow((m_Similarity[VectorSliceIndex][i]-m_Mean[VectorSliceIndex]),2);
            }

            m_StandardDeviation[VectorSliceIndex]=sqrt(sum/double(m_NumberOfVolumes));

            if(m_VerboseMod)
                btkCoutMacro(" - std = "<<m_StandardDeviation[VectorSliceIndex]);
        }

    }
    std::cout<<std::endl;



}



}//namespace btk

#endif

