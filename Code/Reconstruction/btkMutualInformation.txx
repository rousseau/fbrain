/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 03/06/2013
  Author(s): Frederic Champ (champ@unistra.fr)

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

#ifndef __btkMutualInformation_txx
#define __btkMutualInformation_txx

#include "btkMutualInformation.h"
#include "btkImageHelper.h"

namespace btk
{

//----------------------------------------------------------------------------------------
template <typename TImage >
MutualInformation< TImage >
::MutualInformation()
{
    m_VerboseMod=false;
    m_Method="Normalized";
    m_NumberOfBins=100;
    m_PercentageOfBins=0;
    m_binMax= 0;
    m_binMin = 0;
    m_histoSize = 0;

}

//----------------------------------------------------------------------------------------

template <typename TImage >
MutualInformation< TImage >
::~MutualInformation()
{
    //---
}

//----------------------------------------------------------------------------------------

template <typename TImage >
void
MutualInformation< TImage >
::Initialize()
{

    ////////////////////////////////////////////////////////////////////////////
    //
    // Initialize parameters for join histogram.
    //
    ConstIteratorType ItFixed(m_FixedImage, m_FixedImage -> GetLargestPossibleRegion());
    ConstIteratorType ItMoving(m_MovingImage, m_MovingImage -> GetLargestPossibleRegion());
    PixelType Fixed0  = ItFixed.Get();
    PixelType Moving0 = ItMoving.Get();

    if(Fixed0<Moving0)
    {
        m_binMax = Fixed0;
        m_binMin = Fixed0;
    }
    else
    {
        m_binMax = Moving0;
        m_binMin = Moving0;

    }
    for(ItFixed.GoToBegin(), ItMoving.GoToBegin(); !ItFixed.IsAtEnd() && !ItMoving.IsAtEnd(); ++ItFixed, ++ItMoving)
    {
        PixelType FixedVal = ItFixed.Get();
        PixelType MovingVal = ItMoving.Get();

        PixelType MaxVal = FixedVal;
        if(MovingVal > FixedVal)
        {
            MaxVal = MovingVal;
        }

        PixelType MinVal = FixedVal;
        if(MovingVal < FixedVal)
        {
            MinVal = MovingVal;
        }

        if(MaxVal>m_binMax)
        {
            m_binMax=MaxVal;

        }
        if(MinVal<m_binMin)
        {
            m_binMin=MinVal;
        }

    }
    if(m_PercentageOfBins !=0 && m_PercentageOfBins <= 100)
    {
        m_histoSize = m_binMax-m_binMin;
        m_histoSize*=m_PercentageOfBins/100.0;
    }
    else
    {
        m_histoSize=m_NumberOfBins;
    }

}

//----------------------------------------------------------------------------------------


template <typename TImage >
float
MutualInformation< TImage >
::GetValue() throw (ExceptionObject)
{

    ////////////////////////////////////////////////////////////////////////////
    //
    // Join image filter
    //
    JoinFilterPointer joinFilter = JoinFilterType::New();
    joinFilter->SetInput1( m_FixedImage);
    joinFilter->SetInput2( m_MovingImage);

    if(m_VerboseMod)
        std::cout<<"  -> Join images ... ";

    try
    {
        joinFilter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        btkCerrMacro(excp);

    }
    if(m_VerboseMod)
        btkCoutMacro(" Done.");

    ////////////////////////////////////////////////////////////////////////////
    //
    // Compute joint histogram
    //
    HistogramFilterPointer histogramFilter = HistogramFilterType::New();
    histogramFilter->SetInput(  joinFilter->GetOutput() );

    HistogramSizeType size(2);
    size[0] = m_histoSize;  // number of bins for the first  channel
    size[1] = m_histoSize;  // number of bins for the second channel
    histogramFilter->SetHistogramSize( size );

    HistogramMeasurementVectorType binMinimum(3);
    HistogramMeasurementVectorType binMaximum(3);
    binMinimum[0] = m_binMin;
    binMinimum[1] = m_binMin;
    binMinimum[2] = m_binMin;

    binMaximum[0] = m_binMax;
    binMaximum[1] = m_binMax;
    binMaximum[2] = m_binMax;

    histogramFilter->SetHistogramBinMinimum( binMinimum );
    histogramFilter->SetHistogramBinMaximum( binMaximum );

    if(m_VerboseMod)
        std::cout<<"  -> Compute joint histogram ... ";
    try
    {
        histogramFilter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        btkCerrMacro(excp);

    }
    if(m_VerboseMod)
        btkCoutMacro(" Done.");

    ////////////////////////////////////////////////////////////////////////////
    //
    // Calculate joint entropy
    //
    HistogramConstPointer Histogram = histogramFilter->GetOutput();
    HistogramIterator itr = Histogram->Begin();
    HistogramIterator end = Histogram->End();

    const double Sum = Histogram->GetTotalFrequency();
    double JointEntropy = 0.0;

    while( itr != end )
    {
        const double count = itr.GetFrequency();
        if( count > 0.0 )
        {
            const double probability = count / Sum;
            JointEntropy += - probability * vcl_log( probability ) / vcl_log( 2.0 );
        }
        ++itr;
    }
    if(m_VerboseMod)
        std::cout << "  -> Joint entropy = " << JointEntropy << " bits " << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    //
    // Calculate entropy of fixed image
    //
    size[0] = m_histoSize;  // number of bins for the first  channel
    size[1] =   1;  // number of bins for the second channel

    histogramFilter->SetHistogramSize( size );
    histogramFilter->Update();

    itr = Histogram->Begin();
    end = Histogram->End();

    double Entropy1 = 0.0;

    while( itr != end )
    {
        const double count = itr.GetFrequency();
        if( count > 0.0 )
        {
            const double probability = count / Sum;
            Entropy1 += - probability * vcl_log( probability ) / vcl_log( 2.0 );
        }
        ++itr;
    }
    if(m_VerboseMod)
        btkCoutMacro("  -> Fixed image Entropy = " << Entropy1 << " bits " );


    ////////////////////////////////////////////////////////////////////////////
    //
    // Calculate entropy of moving image
    //
    size[0] =   1;  // number of bins for the first channel
    size[1] = m_histoSize;  // number of bins for the second channel

    histogramFilter->SetHistogramSize( size );
    histogramFilter->Update();
    itr = Histogram->Begin();
    end = Histogram->End();

    double Entropy2 = 0.0;

    while( itr != end )
    {
        const double count = itr.GetFrequency();
        if( count > 0.0 )
        {
            const double probability = count / Sum;
            Entropy2 += - probability * vcl_log( probability ) / vcl_log( 2.0 );
        }
        ++itr;
    }
    if(m_VerboseMod)
        btkCoutMacro("  -> Moving image entropy = " << Entropy2 << " bits ");

    ////////////////////////////////////////////////////////////////////////////
    //
    // Calculate mutual information
    //
    float MutualInformation = Entropy1 + Entropy2 - JointEntropy;

    if(m_Method.compare("NonNormalized") == 0)
    {

        ////////////////////////////////////////////////////////////////////////////
        //
        // MutualInformation = Entropy1 + Entropy2 - JointEntropy;
        //
        if(m_VerboseMod)
            btkCoutMacro("  Mutual Information = " << MutualInformation << " bits ");
    }

    else if(m_Method.compare("Normalized") == 0)
    {

        ////////////////////////////////////////////////////////////////////////////
        //
        // Normalized Mutual Information, where the value of Mutual Information gets
        // divided by the mean entropy of the input images.
        //
        MutualInformation = 2.0 * MutualInformation / ( Entropy1 + Entropy2 );
        if(m_VerboseMod)
            btkCoutMacro("  Normalized Mutual Information 1 = " << MutualInformation);
    }
    else
    {
        btkCerrMacro(" -> "<<m_Method<<" : Unknown option for Method ( NonNormalized or Normalized)");
        btkCoutMacro("  Return non noramlized mutual information (default) ... ")

    }

    if(std::isnan(MutualInformation))
        MutualInformation=0.0;

    return -MutualInformation;
}


} // end namespace btk


#endif
