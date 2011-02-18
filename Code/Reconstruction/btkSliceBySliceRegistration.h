/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/03/2010
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

#ifndef __btkSliceBySliceRegistration_h
#define __btkSliceBySliceRegistration_h

#include "itkImageRegistrationMethod.h"
#include "btkRigidRegistration.h"
#include "btkAffineRegistration.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

#include "btkUserMacro.h"

#include "itkCommand.h"

namespace btk
{

using namespace itk;

/** \class SliceBySliceRegistration
 * \brief Describe the class briefly here.
 *
 * Full class description
 * Full class description
 * Full class description

 * \sa ImageRegistrationMethod
 * \ingroup RegistrationFilters
 */

template <typename TImage>
class SliceBySliceRegistration :
public ImageRegistrationMethod <TImage,TImage>
{
public:
  /** Standard class typedefs. */
  typedef SliceBySliceRegistration  Self;
  typedef ImageRegistrationMethod<TImage,TImage>       Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SliceBySliceRegistration, ImageRegistrationMethod);

  /**  Type of the Fixed image. */
  typedef          TImage                               ImageType;
  typedef typename ImageType::Pointer                   ImagePointer;

  typedef typename ImageType::PointType                 PointType;
  typedef typename ImageType::RegionType                RegionType;
  typedef typename ImageType::SizeType               		SizeType;
  typedef typename ImageType::IndexType               	IndexType;

  /**  Type of the Transform . */
  typedef AffineTransform<
                          double,
                          ImageType::ImageDimension >   TransformType;
  typedef typename TransformType::Pointer               TransformPointer;
  typedef typename TransformType::ParametersType        ParametersType;

  typedef itk::Euler3DTransform< double >           RigidTransformType;

  typedef  std::vector<TransformPointer>            TransformPointerArray;


  typedef itk::RigidRegistration< ImageType >       RigidRegistrationType;
  typedef itk::AffineRegistration< ImageType >      AffineRegistrationType;

  typedef vnl_matrix<double> VnlMatrixType;
  typedef itk::Matrix<double,3,3> MatrixType;
  typedef vnl_vector<double> VnlVectorType;

  void StartRegistration();


  TransformPointerArray GetTransformArray()
  {
    return m_TransformArray;
  };

  void SetTransformArray( TransformPointerArray TransformArray)
  {
    m_TransformArray = TransformArray;
    m_TransformArrayIsSet = true;
  };

protected:
  SliceBySliceRegistration();
  virtual ~SliceBySliceRegistration() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Initialize by setting the interconnects between the components.
   */
  void Initialize() throw (ExceptionObject);

private:
  SliceBySliceRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename RigidRegistrationType::Pointer   rigidRegistration;
  typename AffineRegistrationType::Pointer  affineRegistration;

  TransformPointerArray       m_TransformArray;
  bool m_TransformArrayIsSet;

};


} // end namespace btk


#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSliceBySliceRegistration.txx"
#endif

#endif
