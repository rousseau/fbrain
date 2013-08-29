/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/08/2013
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


#include "btkDiffusionDataset.h"

namespace btk
{

DiffusionDataset::DiffusionDataset()
{
    // ----
}

//----------------------------------------------------------------------------------------

DiffusionDataset::~DiffusionDataset()
{

    // ----
}


//----------------------------------------------------------------------------------------

void DiffusionDataset::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    //Superclass::PrintSelf(os,indent);
    os<<indent<< "Origin: "<< this->m_Origin << std::endl;
    os<<indent<< "Spacing: "<< this->m_Spacing << std::endl;
    os<<indent<< "Direction: "<< this->m_Direction << std::endl;
    os<<indent<< "Region: "<< this->m_SequenceRegion<< std::endl;

}


//----------------------------------------------------------------------------------------

void DiffusionDataset::Initialize()
{
    btkCoutMacro(std::endl<<"Dataset initialization and Outliers detection: ");

    this->m_GradientTable       = m_InputSequence -> GetGradientTable();

    this->m_SequenceRegion    = m_InputSequence->GetLargestPossibleRegion();
    this->m_Spacing           = m_InputSequence->GetSpacing();
    this->m_Origin            = m_InputSequence->GetOrigin();
    this->m_Direction         = m_InputSequence->GetDirection();

    std::vector<GradientDirection> GradientTable = m_InputSequence->GetGradientTable();
    std::vector< unsigned short >  BValues       = m_InputSequence->GetBValues();

    DiffusionSequence::SizeType  size4D = this->m_SequenceRegion.GetSize();
    DiffusionSequence::IndexType index3D; //[0,0,0,0] at initialization?
    DiffusionSequence::SizeType  size3D = size4D;
    size3D[3]=0;

    DiffusionSequence::SizeType size = size4D;
    size[3]-=1; // To not consider de B0
    this->m_SequenceRegion.SetSize(size);

    ////////////////////////////////////////////////////////////////////////////
    // B0 extraction
    ////////////////////////////////////////////////////////////////////////////
    Region4DType B0Region;
    B0Region.SetIndex( index3D );
    B0Region.SetSize( size3D );

    ExtractorPointer B0Extractor = ExtractorType::New();
    B0Extractor -> SetInput( m_InputSequence );
    B0Extractor -> SetExtractionRegion( B0Region );
    B0Extractor -> SetDirectionCollapseToSubmatrix( );
    B0Extractor -> Update();

    m_B0 = B0Extractor->GetOutput();
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Similarity calculation to get outliers -> S2SSIMILARITY HAVE TO BE CHANGED ...
    ////////////////////////////////////////////////////////////////////////////
    //
    // Calcul similarities for the whole sequence
    // The similarity calculated is the normalized mutual information
    // between the slice of the B0 and the slice for all DWI volume
    //
    SimilarityPointer SimilarityFilter = SimilarityType::New();
    SimilarityFilter -> SetInputSequence(m_InputSequence.GetPointer());
    //SimilarityFilter -> SetVerboseMod(true);
    SimilarityFilter -> Update();

    std::vector<std::vector<float > > similarity        = SimilarityFilter -> GetSimilarity();
    std::vector<float>                mean              = SimilarityFilter -> GetMean();
    std::vector<float>                standardDeviation = SimilarityFilter -> GetStandardDeviation();
    ////////////////////////////////////////////////////////////////////////////

    btkCoutMacro("  -> Outliers detected: ");
    ////////////////////////////////////////////////////////////////////////////
    // Dwi slices extraction and outliers detection
    ////////////////////////////////////////////////////////////////////////////
    //
    // Outliers detections:
    // Reference: Dubois et Al, "Correction strategy for infants'
    //            diffusion-weighted images corrupted with motion"
    // A slice of a DWI volume is considered as an outlier if the similarity
    // of a slice is:
    //  > mean+3*StandardDeviation
    //  < mean-3*StandardDeviation
    // mean              = the mean similarity for this slice in the input sequence
    // StandardDeviation = the standard deviation for this slice in the input sequence
    //
    size3D[2]=1; // slice size
    bool hasOutliers = false;

    for(unsigned int j=0; j<size4D[3]-1; j++)
    {
        for(unsigned int i=0; i<size4D[2]; i++)
        {

            index3D[2]=i;
            index3D[3]=j+1;

            Region4DType sliceRegion;
            sliceRegion.SetIndex( index3D );
            sliceRegion.SetSize( size3D );

            ExtractorPointer Extractor = ExtractorType::New();
            Extractor -> SetInput( m_InputSequence );
            Extractor -> SetExtractionRegion( sliceRegion );
            Extractor -> SetDirectionCollapseToSubmatrix( );
            Extractor -> Update();

            DiffusionSlice::Pointer slice = DiffusionSlice::New();
            slice -> SetImage(Extractor->GetOutput());
            slice -> SetGradientDirection(GradientTable[j+1]);
            slice -> SetBValue(BValues[j+1]);

            if(!m_InitialTransforms.empty())
            {
                slice -> Transform(m_InitialTransforms[j]);
            }

            if(similarity[i][j]>mean[i]+3*standardDeviation[i] ||
                    similarity[i][j]<mean[i]-3*standardDeviation[i] )
            {
                btkCoutMacro("      Volume "<<index3D[3]<<" - Slice "<<index3D[2]<<" - Similarity value = "<<similarity[i][j]);
                slice -> SetOutlierStatus(true);
                bool hasOutliers = true;
            }
            m_Data.push_back(slice);
        }
        ////////////////////////////////////////////////////////////////////////////


    }
    if(hasOutliers==false)
    {
        btkCoutMacro("      No outliers detected ... ");
    }
    btkCoutMacro("");
}

//----------------------------------------------------------------------------------------

void DiffusionDataset::RemoveOutliers()
{
    ////////////////////////////////////////////////////////////////////////////
    // Delete slices which have an OutlierStatus = true
    ////////////////////////////////////////////////////////////////////////////
    for (DataVectorType::iterator it = m_Data.begin() ; it != m_Data.end(); ++it)
    {
        if((*it)->IsOutlier())
        {
            m_Data.erase(it);
        }
    }
}

} // end namespace btk
