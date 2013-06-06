/*==========================================================================

  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique

  Date: 04/06/2013
  Author(s): Frederic Champ (champ(at)unistra.fr)

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

#ifndef _btkOutlierCorrectionFilter_txx
#define _btkOutlierCorrectionFilter_txx

#include "btkOutlierCorrectionFilter.h"
#include "btkDiffusionSequenceHelper.h"

namespace btk
{

template <typename TImage >
OutlierCorrectionFilter<TImage>::OutlierCorrectionFilter()
{
    m_Factor = 3;
    m_NumberOfIterations = 10;
    m_VerboseMod = false;
    m_InputSequence = TSequence::New();
    m_Delimiter=";";
    m_Method ="SH";
    m_Rgra = 0.15;
    m_Sigma = 1.0;
    m_EllipsoidSize = 1.0;

}

//----------------------------------------------------------------------------------------

template <typename TImage >
OutlierCorrectionFilter<TImage>::~OutlierCorrectionFilter()
{
    //---
}

//----------------------------------------------------------------------------------------

template < typename TImage >
typename OutlierCorrectionFilter<TImage>::SequencePointer
OutlierCorrectionFilter<TImage>::ExtractOutlier(SequenceConstPointer InputSequence, std::vector<unsigned int> CorruptedVolumes, unsigned int SliceNumber)
{

    ////////////////////////////////////////////////////////////////////////////
    //
    // Initializations
    //
    SequenceRegionType InputRegion = InputSequence -> GetLargestPossibleRegion();
    SequenceSizeType   InputSize   = InputRegion.GetSize();
    SequenceIndexType  SliceIndex  = InputRegion.GetIndex();

    std::vector< GradientDirection > GradientTable = InputSequence -> GetGradientTable();
    std::vector< unsigned short >    BValues       = InputSequence -> GetBValues();
    std::vector<unsigned int>        VolumeIndexes;

    SliceIndex[2]=SliceNumber;
    unsigned int NumberOfVolumes = InputSize[3];

    ////////////////////////////////////////////////////////////////////////////
    //
    // Create a vector with volume indexes.
    //
    for(unsigned int i=0; i< NumberOfVolumes; i++)
    {
        VolumeIndexes.push_back(i);
    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //  Then, if there are outliers for this slice in the DWI sequence,
    //  erase their volumes indexes.
    //
    if(CorruptedVolumes.size()>0)
    {
        for(unsigned int j=0;j<CorruptedVolumes.size();j++)
        {
            unsigned int VolumeIndex = CorruptedVolumes[j];
            VolumeIndexes.erase(VolumeIndexes.begin()+VolumeIndex);
            GradientTable.erase(GradientTable.begin()+VolumeIndex);
            BValues.erase(BValues.begin()+VolumeIndex);

            // btkCoutVariable(VolumeIndex);
        }

    }

    ////////////////////////////////////////////////////////////////////////////
    //
    //  Extract volumes whose indexes are present in the VolumesIndexes vector
    //  Then a new sequence of one slice will be created without outliers.
    //
    JoinSequencePointer joiner = JoinSequenceFilter::New();
    joiner -> SetOrigin( 0.0 );
    joiner -> SetSpacing( m_InputSequence -> GetSpacing()[3]);

    for(unsigned int i=0; i< VolumeIndexes.size(); i++)
    {
        SequenceRegionType  SliceRegion = InputRegion;
        SliceRegion.SetSize(2,1);
        SliceRegion.SetSize(3,0);
        SliceRegion.SetIndex(SliceIndex);
        SliceRegion.SetIndex(3,VolumeIndexes[i]);

        ImageExtractorPointer extractor = ImageExtractor::New();
        extractor -> SetInput(InputSequence);
        extractor -> SetExtractionRegion(SliceRegion);
        extractor -> SetDirectionCollapseToSubmatrix( );
        extractor -> Update();

        joiner -> SetInput(i,extractor -> GetOutput());
    }
    joiner -> Update();

    SequencePointer SliceSequence = joiner -> GetOutput();
    SliceSequence -> SetGradientTable(GradientTable);
    SliceSequence -> SetBValues(BValues);

    ////////////////////////////////////////////////////////////////////////
    //
    // Change the region slice index to be equal to 0
    //
    SequenceRegionType NewRegion = SliceSequence->GetLargestPossibleRegion();
    NewRegion.SetIndex(2,0);
    SliceSequence -> SetRegions(NewRegion);

    return SliceSequence;

}

//----------------------------------------------------------------------------------------
template < typename TImage >
void
OutlierCorrectionFilter<TImage>::GetOutliers()  throw (ExceptionObject)
{
    btkCoutMacro(std::endl<<"Outliers detection: ");

    ////////////////////////////////////////////////////////////////////////////
    //
    // Initializations
    //
    m_OutliersIndexes.clear();


    ////////////////////////////////////////////////////////////////////////////
    //
    // Calcul similarities for the whole sequence
    // The similarity calculated is the normalized mutual information
    // between the slice of the B0 and the slice for all DWI volume
    //
    S2SSimilarityPointer SimilarityFilter = S2SSimilarityType::New();
    SimilarityFilter -> SetInputSequence(m_InputSequence);
    //SimilarityFilter -> SetVerboseMod(m_VerboseMod);
    try
    {
        SimilarityFilter -> Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        btkCerrMacro(excp);
    }

    SimilarityS2SType Similarity        = SimilarityFilter -> GetSimilarity();
    S2SVectorType     Mean              = SimilarityFilter -> GetMean();
    S2SVectorType     StandardDeviation = SimilarityFilter -> GetStandardDeviation();

    if(m_SimilarityFileName.compare("") != 0)
    {
        SimilarityFilter -> WriteData(m_SimilarityFileName);
    }


    ////////////////////////////////////////////////////////////////////////////
    //
    // Outliers detections:
    // A slice of a DWI volume is considered as an outlier if the similarity
    // for this slice is:
    //  > mean+m_Factor*StandardDeviation
    //  < mean-m_Factor*StandardDeviation
    // mean              = the mean similarity for this slice in the input sequence
    // StandardDeviation = the standard deviation for this slice in the input sequence
    //
    unsigned int NumberOfSlices = Similarity.size();
    unsigned int NumberOfVolumes= Similarity[0].size();

    SequenceIndexType index;
    index[0]=0;
    index[1]=0;
    for (unsigned short j=0; j < NumberOfSlices; j++)
    {
        for (unsigned short i=0; i < NumberOfVolumes; i++)
        {
            if(Similarity[j][i]>Mean[j]+m_Factor*StandardDeviation[j]
                    ||Similarity[j][i]<Mean[j]-m_Factor*StandardDeviation[j] )
            {
                index[2]=j;
                index[3]=i+1;
                m_OutliersIndexes.push_back(index);
                m_BoolIndexes[j][i]= true;

                btkCoutMacro("      Slice "<<index[2]<<" - Volume "<<index[3]<<" - Similarity value = "<<Similarity[j][i]);
            }
        }
    }
    if(m_OutliersIndexes.empty())
    {
        btkCoutMacro("      No outliers detected");
    }
}

//----------------------------------------------------------------------------------------

template < typename TImage >
void
OutlierCorrectionFilter<TImage>::CorrectOutliers(SequenceConstPointer InputSequence)
{
    //TODO clean

    btkCoutMacro(std::endl<<"Outliers correction: ");

    ////////////////////////////////////////////////////////////////////////////
    //
    // Initializations
    //
    SequenceRegionType Region4D    = InputSequence -> GetLargestPossibleRegion();
    SequenceIndexType  Index4D     = Region4D.GetIndex();
    SequenceSizeType   Size4D      = Region4D.GetSize();

    std::vector< unsigned short>           BValues       = InputSequence->GetBValues();
    std::vector< btk::GradientDirection >  GradientTable = InputSequence->GetGradientTable();

    unsigned int NumberOfSlices = Size4D[2];

    if(m_Method.compare("SH") ==0)
    {
        unsigned int nb_loop=0;
        unsigned int i;

        #pragma omp parallel for private(i) schedule(dynamic)
        for(i=0; i< NumberOfSlices; i++)
        {
            SequenceIndexType SliceIndex = Index4D;
            SliceIndex[2]=i;

            ////////////////////////////////////////////////////////////////////////////
            //
            // check if there are some outliers for this slice and get volume indexes.
            //
            std::vector<unsigned int> CorruptedVolumes;
            for(unsigned int j=0; j<m_OutliersIndexes.size(); j++ )
            {
                if(i == m_OutliersIndexes[j][2])
                {
                    CorruptedVolumes.push_back(m_OutliersIndexes[j][3]);
                }
            }

            ////////////////////////////////////////////////////////////////////////////
            //
            // Create a sequence of one slice without corrupted volumes if exist.
            //
            SequencePointer SliceSequence = ExtractOutlier(InputSequence,CorruptedVolumes,i);

            ////////////////////////////////////////////////////////////////////////////
            //
            // Spherical harmonic decomposition
            //
            DecompositionPointer Decomposition = DecompositionType::New();
            Decomposition -> SetInput(SliceSequence);
            Decomposition -> SetSphericalHarmonicsOrder(4);
            Decomposition -> Update();

            ////////////////////////////////////////////////////////////////////////////
            //
            // Get signal models
            //
            ModelPointer  Model = ModelType::New();
            Model -> SetInputModelImage(Decomposition->GetOutput());
            Model -> SetBValue(BValues[1]);
            Model -> SetSphericalResolution(300);
            Model -> Update();

            ////////////////////////////////////////////////////////////////////////////
            //
            // Replace outliers values by values estimated
            //

            SequenceSizeType Slice4DSize;
            Slice4DSize[0]=Size4D[0];
            Slice4DSize[1]=Size4D[1];
            Slice4DSize[2]=1;
            Slice4DSize[3]=Size4D[3];

            SequenceRegionType SliceRegion;
            SliceRegion.SetSize(Slice4DSize);
            SliceRegion.SetIndex(SliceIndex);

            SequenceIterator It (m_CorrectedSequence, SliceRegion);

            for(It.GoToBegin();!It.IsAtEnd();++It)
            {
                IndexContinuous index3D;
                SequenceIndexType   index4D = It.GetIndex();

                index3D[0] = index4D[0];
                index3D[1] = index4D[1];
                index3D[2] = 0;

                SequenceIndexType indexB0 = index4D;
                indexB0[3]=0; //Because we have to multiplie the estimated signal by the B0's voxels values.

                if(index4D[3]>0)
                {
                    It.Set(Model->SignalAt(index3D,GradientTable[index4D[3]]) * InputSequence->GetPixel(indexB0 ));
                }
            }
            nb_loop++;
            std::cout<<"\r  -> "<<nb_loop <<" slices treated (total "<< NumberOfSlices <<") ... "<<std::flush;
        }
        btkCoutMacro(" Done.");
        std::cout<<std::endl;
    }

    else if(m_Method.compare("WSH") ==0)
    {

        btkCoutMacro("  Sigma = "<<m_Sigma);
        btkCoutMacro("  Ellipsoid size = "<<m_EllipsoidSize);

        ////////////////////////////////////////////////////////////
        //
        // Remove B0
        //
        JoinSequencePointer InJoiner = JoinSequenceFilter::New();
        InJoiner -> SetOrigin( 0.0 );
        InJoiner -> SetSpacing( InputSequence -> GetSpacing()[3]);

        for(unsigned int i=1; i< m_Size4D[3]; i++)
        {
            SequenceRegionType  ImageRegion = InputSequence -> GetLargestPossibleRegion();
            ImageRegion.SetSize(3,0);
            ImageRegion.SetIndex(3,i);

            ImageExtractorPointer extractor = ImageExtractor::New();
            extractor -> SetInput(InputSequence);
            extractor -> SetExtractionRegion(ImageRegion);
            extractor -> SetDirectionCollapseToSubmatrix( );
            extractor -> Update();

            InJoiner -> SetInput(i-1,extractor -> GetOutput());
        }
        InJoiner -> Update();

        SequencePointer GradientSequence = InJoiner -> GetOutput();
        std::vector< btk::GradientDirection > NewGradientTable = InputSequence -> GetGradientTable();
        NewGradientTable.erase(NewGradientTable.begin(),NewGradientTable.begin()+1);
        GradientSequence -> SetGradientTable(NewGradientTable);


        ////////////////////////////////////////////////////////////
        //
        // Weighted Decomposition
        //
        WeighteEstimationFilterPointer EstimationFilter = WeightedEstimationFilter::New();
        EstimationFilter -> SetInputSequence(GradientSequence.GetPointer());
        EstimationFilter -> SetEllipsoidSize(m_EllipsoidSize);
        EstimationFilter -> SetBoolIndexes(m_BoolIndexes);
        EstimationFilter -> SetVerboseMod(m_VerboseMod);
        EstimationFilter -> SetSigma(m_Sigma);
        EstimationFilter -> Update();

        SequencePointer CorrectedSequence = EstimationFilter -> GetOutput();


        ////////////////////////////////////////////////////////////
        //
        // add B0
        //
        JoinSequencePointer OutJoiner = JoinSequenceFilter::New();
        OutJoiner -> SetOrigin( 0.0 );
        OutJoiner -> SetSpacing( InputSequence -> GetSpacing()[3]);

        SequenceRegionType  B0Region = InputSequence -> GetLargestPossibleRegion();
        B0Region.SetSize(3,0);
        B0Region.SetIndex(3,0);

        ImageExtractorPointer B0Extractor = ImageExtractor::New();
        B0Extractor -> SetInput(InputSequence);
        B0Extractor -> SetExtractionRegion(B0Region);
        B0Extractor -> SetDirectionCollapseToSubmatrix( );
        B0Extractor -> Update();

        OutJoiner -> SetInput(0,B0Extractor -> GetOutput());

        for(unsigned int i=0; i< m_Size4D[3]-1; i++)
        {
            SequenceRegionType  ImageRegion = CorrectedSequence -> GetLargestPossibleRegion();
            ImageRegion.SetSize(3,0);
            ImageRegion.SetIndex(3,i);

            ImageExtractorPointer extractor = ImageExtractor::New();
            extractor -> SetInput(CorrectedSequence);
            extractor -> SetExtractionRegion(ImageRegion);
            extractor -> SetDirectionCollapseToSubmatrix( );
            extractor -> Update();

            OutJoiner -> SetInput(i+1,extractor -> GetOutput());
        }
        OutJoiner -> Update();

        m_CorrectedSequence = OutJoiner -> GetOutput();
        btkCoutMacro(std::endl);

    }
    else if(m_Method.compare("RBF") ==0)
    {

        btkCoutMacro("  Correct only outlier slices ! (old)");

        double theta;
        double phi;
        float rspa = 1.0;

        btkCoutMacro("  rspa = "<<rspa<<" ( fixed )");
        btkCoutMacro("  rgra = "<<m_Rgra);

        ////////////////////////////////////////////////////////////////////////////
        //
        // RBF interpolator initialization
        //
        RBFInterpolatorPointer Interpolator = RBFInterpolatorType::New();
        Interpolator -> SetInputImage( m_CorrectedSequence );
        Interpolator -> SetGradientTable(GradientTable);

        ////////////////////////////////////////////////////////////////////////////
        //
        // Replace outliers values by values estimated
        //
        SequenceSizeType Slice3DSize;
        Slice3DSize[0]=Size4D[0];
        Slice3DSize[1]=Size4D[1];
        Slice3DSize[2]=1;
        Slice3DSize[3]=1; //treat each outliers separatly

        unsigned int i;
        #pragma omp parallel for private(i) schedule(dynamic)
        for(i=0; i < m_OutliersIndexes.size(); i++)
        {
            SequenceRegionType SliceRegion;

            SliceRegion.SetSize(Slice3DSize);
            SliceRegion.SetIndex(m_OutliersIndexes[i]);

            Interpolator -> Initialize(SliceRegion);
            SequenceIterator It (m_CorrectedSequence, SliceRegion);

            for(It.GoToBegin();!It.IsAtEnd();++It)
            {
                SequencePointType point;
                SequenceIndexType  index4D = It.GetIndex();

                Interpolator -> GetGradientDirection(index4D[3],theta,phi);
                m_CorrectedSequence -> TransformIndexToPhysicalPoint(index4D,point);

                float value = Interpolator -> EvaluateAt(index4D, theta, phi, rspa, m_Rgra, 0);
                It.Set((short)value);

            }
        }

    }
}

//----------------------------------------------------------------------------------------

template < typename TImage >
void OutlierCorrectionFilter<TImage>::Update() throw (ExceptionObject)
{

    SequenceRegionType Region4D = m_InputSequence -> GetLargestPossibleRegion();
    m_Size4D      = Region4D.GetSize();

    ////////////////////////////////////////////////////////////////////////////
    //
    // Allocate the corrected sequence
    //
    m_CorrectedSequence = btk::ImageHelper<TSequence>::DeepCopy(m_InputSequence);
    m_CorrectedSequence -> SetBValues(m_InputSequence->GetBValues());
    m_CorrectedSequence -> SetGradientTable(m_InputSequence->GetGradientTable());

    ////////////////////////////////////////////////////////////////////////////
    //
    // Initialize the boolean vector of outliers.
    //
    m_BoolIndexes.resize(m_Size4D[2]);
    for(unsigned int i=0; i<m_BoolIndexes.size(); i++)
    {
        m_BoolIndexes[i].resize(m_Size4D[3]-1,false);
    }

    ////////////////////////////////////////////////////////
    //
    // Outliers exlucsion and correction
    //
    this->GetOutliers();
    this->CorrectOutliers(m_InputSequence);

    ////////////////////////////////////////////////////////
}

} //namespace btk

#endif

