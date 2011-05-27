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
  /*

  In nifti format the gradients are in image coordinates,
  i.e. g = S*[i j k]', where S is the diagonal matrix of
  spacings and [i j k]' are the voxel indexes.

  WARNING: sometimes dcm2nii generates a table of gradients
  in RAS world coordinates

  */

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

}

template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::RotateGradientsInWorldCoordinates()
{


  if (!m_Transform)
  {
    itkExceptionMacro(<<"Transform is not present");
  }

  vnl_matrix< double > Direction;
  vnl_matrix< double > DirectionInv;

  Direction = m_Image -> GetDirection().GetVnlMatrix().extract(3,3,0,0);
  DirectionInv = vnl_inverse(Direction);

  // FIXME Here the transformation provided should be used to rotate the gradients
  // It's the calling program who should take care of providing the correct transform
  // After correcting this, all the programs using thsi class should be modified accordingly


  vnl_matrix<double> R = m_Transform -> GetMatrix().GetVnlMatrix();
  vnl_matrix<double> Rinv(3,3);
  Rinv = vnl_inverse(R);

  for (unsigned int i=0; i < m_GradientTable.rows(); i++)
  {

    vnl_matrix<double> R = m_Transform -> GetMatrix().GetVnlMatrix();
    vnl_vector<double> grad = m_GradientTable.get_row(i);
    m_GradientTable.set_row(i, DirectionInv*Rinv*Direction*grad);

  }

}

template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::RotateGradients()
{

  for (unsigned int i=0; i < m_GradientTable.rows(); i++)
  {
    vnl_vector<double> grad = m_GradientTable.get_row(i);
    m_GradientTable.set_row(i, m_RotationMatrix*grad);
  }

}

template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::TransformGradientsToImageCoordinates()
{

  /*

  Please be aware that this method assumes that gradients are
  in RAS world coordinates, since it's the way in which sometimes
  dcm2nii stores gradients. If gradients are expressed in LPS the
  result will be incorrect.

  */

  std::cout << "Table in RAS WC = " << std::endl;
  std::cout <<  m_GradientTable << std::endl << std::endl;

  vnl_matrix< double > Direction;
  vnl_matrix< double > DirectionInv;

  Direction = m_Image -> GetDirection().GetVnlMatrix().extract(3,3,0,0);

  // The direction matrix provided by ITK transforms points into
  // LPS, so a correction must be performed by multiplying 1st and
  // 2nd rows by -1 for obtaining points in RAS

  Direction.scale_row(0,-1);
  Direction.scale_row(1,-1);

  DirectionInv = vnl_inverse(Direction);

  for (unsigned int i=0; i < m_GradientTable.rows(); i++)
  {
    vnl_vector<double> grad = m_GradientTable.get_row(i);
    grad = DirectionInv*grad;
    grad.normalize();
    m_GradientTable.set_row(i, grad);
  }

  std::cout << "Table in image coordinates = " << std::endl;
  std::cout <<  m_GradientTable << std::endl << std::endl;

}


template < typename TInputImage >
void
DiffusionGradientTable<TInputImage>
::TransformGradientsToWorldCoordinates()
{

  /*

  Please be aware that this method assumes that gradients are
  in RAS world coordinates, since it's the way in which sometimes
  dcm2nii stores gradients.

  */

  std::cout << "Table in image coordinates = " << std::endl;
  std::cout <<  m_GradientTable << std::endl << std::endl;

  vnl_matrix< double > Direction;
  Direction = m_Image -> GetDirection().GetVnlMatrix().extract(3,3,0,0);

  // The direction matrix provided by ITK transforms points into
  // LPS, so a correction must be performed by multiplying 1st and
  // 2nd rows by -1 for obtaining points in RAS

  Direction.scale_row(0,-1);
  Direction.scale_row(1,-1);

  for (unsigned int i=0; i < m_GradientTable.rows(); i++)
  {
    vnl_vector<double> grad = m_GradientTable.get_row(i);
    grad = Direction*grad;
    grad.normalize();
    m_GradientTable.set_row(i, grad);
  }

  std::cout << "Table in RAS WC = " << std::endl;
  std::cout <<  m_GradientTable << std::endl << std::endl;

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
    for (unsigned int j=0; j < m_NumberOfGradients; j++)
    {

      if ( j < (m_NumberOfGradients-1) )
      {
        fprintf( fw, "%.14lf ", m_GradientTable(j,i) );
      }
      else
      {
        fprintf( fw, "%.14lf\n", m_GradientTable(j,i) );
      }

    }

  }

  fclose (fw);

}


} // end namespace btk

#endif
