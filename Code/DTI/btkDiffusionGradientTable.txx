/*=========================================================================

Purpose: store gradient directions and make all operations required after
transforming the image.

=========================================================================*/
#ifndef __btkDiffusionGradientTable_txx
#define __btkDiffusionGradientTable_txx

#include "btkDiffusionGradientTable.h"
#include "vnl/vnl_inverse.h"

namespace btk
{

//----------------------------------------------------------------------
// Construct without computing moments
template<class TImage>
DiffusionGradientTable<TImage>::DiffusionGradientTable(void)
{
  m_NumberOfGradients = 0;
  m_Transform = 0;
}

//----------------------------------------------------------------------
// Destructor
template<class TImage>
DiffusionGradientTable<TImage>::
~DiffusionGradientTable()
{
}

template<class TInputImage>
void
DiffusionGradientTable<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "Number of gradients: " << m_NumberOfGradients << std::endl;
}

template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::LoadFromFile( const char* input )
{

  m_GradientTable.set_size( m_NumberOfGradients,3 );

  // Fill gradient table in cartesian coordinates

  FILE* fr;
  fr = fopen( input, "r" );

  float value;
  unsigned ncol = 1;
  unsigned nrow = 1;

  for (unsigned int j=1; j <= 3*m_NumberOfGradients; j++)
  {
    if ( ( j % m_NumberOfGradients ) == 0 )
    {
      fscanf( fr, "%f\n", &value);
      m_GradientTable(nrow-1,ncol-1)=value;
      nrow = 0;
      ncol++;
    } else
    {
      fscanf( fr, "%f ", &value);
      m_GradientTable(nrow-1,ncol-1)=value;
    }
    nrow++;
  }

  fclose (fr);

  std::cout << m_GradientTable << std::endl;

}

template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::RotateGradients()
{

  if (!m_Transform)
  {
    itkExceptionMacro(<<"Transform is not present");
  }

  vnl_matrix< double > Direction;
  vnl_matrix< double > Spacing(3,3);
  vnl_matrix< double > DirectionInv;
  vnl_matrix< double > SpacingInv;

  Direction = m_Image -> GetDirection().GetVnlMatrix().extract(3,3,0,0);
  DirectionInv = vnl_inverse(Direction);

  typename ImageType::SpacingType spacing = m_Image -> GetSpacing();

  Spacing.fill(0.0);
  Spacing(0,0) = spacing[0]; Spacing(1,1) = spacing[1]; Spacing(2,2) = spacing[2];
  SpacingInv = vnl_inverse(Spacing);

/*  std::cout << "Direccion = " << std::endl <<  Direction << std::endl;
  std::cout << "Direccion inverse= " << std::endl <<  DirectionInv << std::endl;
  std::cout << "Spacing = " << std::endl <<  Spacing << std::endl;
  std::cout << "Spacing inverse = " << std::endl <<  SpacingInv << std::endl; */


  vnl_matrix<double> R = m_Transform -> GetMatrix().GetVnlMatrix();
  vnl_matrix<double> Rinv(3,3);
  Rinv = vnl_inverse(R);

  for (unsigned int i=0; i < m_GradientTable.rows(); i++)
  {

    vnl_matrix<double> R = m_Transform -> GetMatrix().GetVnlMatrix();
    vnl_vector<double> grad = m_GradientTable.get_row(i);
    m_GradientTable.set_row(i, SpacingInv*DirectionInv*Rinv*Direction*Spacing*grad);

  }

//  std::cout << m_GradientTable << std::endl;

}

template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::TransformGradientsToImageCoordinates()
{

  vnl_matrix< double > Direction;
  vnl_matrix< double > Spacing(3,3);
  vnl_matrix< double > DirectionInv;
  vnl_matrix< double > SpacingInv;

  Direction = m_Image -> GetDirection().GetVnlMatrix().extract(3,3,0,0);
  DirectionInv = vnl_inverse(Direction);

  typename ImageType::SpacingType spacing = m_Image -> GetSpacing();

  Spacing.fill(0.0);
  Spacing(0,0) = spacing[0]; Spacing(1,1) = spacing[1]; Spacing(2,2) = spacing[2];
  SpacingInv = vnl_inverse(Spacing);

  std::cout << "Direccion = " << std::endl <<  Direction << std::endl;
  std::cout << "Direccion inverse= " << std::endl <<  DirectionInv << std::endl;
  std::cout << "Spacing = " << std::endl <<  Spacing << std::endl;
  std::cout << "Spacing inverse = " << std::endl <<  SpacingInv << std::endl;


  for (unsigned int i=0; i < m_GradientTable.rows(); i++)
  {
    vnl_vector<double> grad = m_GradientTable.get_row(i);
    m_GradientTable.set_row(i, SpacingInv*DirectionInv*grad);
  }

}


template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::SaveToFile( const char* output )
{

  // Fill gradient table in cartesian coordinates

  FILE* fw;
  fw = fopen( output, "w" );

  for (unsigned int i=0; i<=2; ++i)
  {
    for (unsigned int j=0; j <= m_NumberOfGradients; j++)
    {

      if ( j < m_NumberOfGradients )
      {
        fprintf( fw, "%f ", m_GradientTable(j,i) );
      }
      else
      {
        fprintf( fw, "%f\n", m_GradientTable(j,i) );
      }

    }
  }

  fclose (fw);

}


} // end namespace btk

#endif
