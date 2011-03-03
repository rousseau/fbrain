/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 19/01/2011
  Author(s): Estanislao Oubel (oubel@unistra.fr)

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

#ifndef __btkDiffusionGradientTable_h
#define __btkDiffusionGradientTable_h

#include "itkEuler3DTransform.h"
#include "itkMacro.h"
#include "itkImage.h"
#include "itkSpatialObject.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"


namespace btk
{

  using namespace itk;

/** \class DiffusionGradientTable
 * \brief Manage gradient table
 *
 * Add long description here
 * Add long description here
 * Add long description here
 *
 * \ingroup Operators
 *
 */
template < class TImage >
class ITK_EXPORT DiffusionGradientTable : public Object
{
public:
  /** Standard class typedefs. */
  typedef DiffusionGradientTable<TImage>   Self;
  typedef Object                           Superclass;
  typedef SmartPointer<Self>               Pointer;
  typedef SmartPointer<const Self>         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(DiffusionGradientTable, Object);

  /** Extract the dimension of the image. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TImage::ImageDimension);

  /** Standard scalar type within this class. */
  typedef double                       ScalarType;

  /** Standard vector type within this class. */
  typedef Vector<ScalarType,itkGetStaticConstMacro(ImageDimension)> VectorType;

  /** Spatial Object type within this class. */
  typedef SpatialObject< itkGetStaticConstMacro(ImageDimension) > SpatialObjectType;

  /** Spatial Object member types used within this class. */
  typedef typename SpatialObjectType::Pointer      SpatialObjectPointer;
  typedef typename SpatialObjectType::ConstPointer SpatialObjectConstPointer;

  /** Standard matrix type within this class. */
  typedef Matrix<ScalarType,
                 itkGetStaticConstMacro(ImageDimension),
                 itkGetStaticConstMacro(ImageDimension)>   MatrixType;

  /** Standard image type within this class. */
  typedef TImage ImageType;

  /** Standard image type pointer within this class. */
  typedef typename ImageType::Pointer      ImagePointer;
  typedef typename ImageType::ConstPointer ImageConstPointer;

  /** Affine transform for mapping to and from principal axis */
  typedef Euler3DTransform<double> TransformType;
  typedef typename TransformType::Pointer      TransformPointer;

  // Load gradient table from file
  void LoadFromFile( const char* input );

  void RotateGradients();

  // Transform gradients in word coordinates to image coordinates
  void TransformGradientsToImageCoordinates();

  // Save gradient table to file
  void SaveToFile( const char* output );

  // Set/Get transform
  itkSetObjectMacro( Transform, TransformType );
  itkGetObjectMacro( Transform, TransformType );

  // Set/Get transform
  itkSetObjectMacro( Image, ImageType );
  itkGetObjectMacro( Image, ImageType );

  // Set/Get number of gradients
  itkSetMacro( NumberOfGradients, unsigned int );
  itkGetMacro( NumberOfGradients, unsigned int );

protected:
  DiffusionGradientTable();
  virtual ~DiffusionGradientTable();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  DiffusionGradientTable(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  vnl_matrix< double >      m_GradientTable;
  unsigned int              m_NumberOfGradients;
  TransformPointer          m_Transform;
  ImagePointer              m_Image;


};  // class DiffusionGradientTable

} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkDiffusionGradientTable.txx"
#endif

#endif /* __btkDiffusionGradientTable_h */
