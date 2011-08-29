/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 16/08/2011
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

#ifndef __btkSliceBySliceTransform_h
#define __btkSliceBySliceTransform_h

#include "itkTransform.h"
#include "itkEuler3DTransform.h"
#include "itkImage.h"
#include "itkContinuousIndex.h"
#include "list"

namespace btk
{
using namespace itk;

template <class TScalarType,unsigned int NDimensions=3>
class SliceBySliceTransform  : public Transform<TScalarType,NDimensions,NDimensions>
{
public:
  /** Standard class typedefs. */
  typedef SliceBySliceTransform  Self;
  typedef Transform<TScalarType,NDimensions,NDimensions> Superclass;
  typedef Euler3DTransform< TScalarType > TransformType;

  typedef Image< short,NDimensions > ImageType;
  typedef typename ImageType::Pointer ImagePointerType;

  typedef ContinuousIndex<double, NDimensions > ContinuousIndexType;

  typedef SmartPointer< Self >   Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  typedef SmartPointer<TransformType> TransformPointer;
  typedef std::vector<TransformPointer> TransformPointerList;
  typedef typename TransformPointerList::const_iterator TransformPointerListConstIterator;
  typedef std::vector<TransformBase::Pointer> TransformBasePointerList;

  typedef typename Superclass::InputPointType InputPointType;
  typedef typename Superclass::OutputPointType OutputPointType;

  typedef typename Superclass::InputVectorType InputVectorType;
  typedef typename Superclass::OutputVectorType OutputVectorType;

  typedef typename Superclass::InputVnlVectorType InputVnlVectorType;
  typedef typename Superclass::OutputVnlVectorType OutputVnlVectorType;

  typedef typename Superclass::InputCovariantVectorType InputCovariantVectorType;
  typedef typename Superclass::OutputCovariantVectorType OutputCovariantVectorType;

  typedef typename Superclass::JacobianType JacobianType;

  typedef typename Superclass::ParametersType ParametersType;

  /** New method for creating an object using a factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( SliceBySliceTransform, Transform );

  /** Method to transform a point. */
  virtual OutputPointType TransformPoint(const InputPointType& p ) const;

  /**  Method to transform a vector. */
  virtual OutputVectorType TransformVector(const InputVectorType& p) const;

  /**  Method to transform a vnl_vector. */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &p ) const;

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType TransformCovariantVector( const InputCovariantVectorType &p ) const;

  /** Compute Jacobian (product of component Jacobians). Not thread safe!! */
  virtual const JacobianType & GetJacobian(const InputPointType &p ) const;

  /** NOT IMPLEMENTED: */
  virtual void GetJacobianWithRespectToParameters( const InputPointType  &p,
                                 JacobianType & jacobian) const
  {
    itkExceptionMacro("GetJacobianWithRespectToParameters "
                        "not yet implemented.");
  }

  /** NOT IMPLEMENTED: */
  virtual void GetJacobianWithRespectToPosition(
                                       const InputPointType & p,
                                       JacobianType &j ) const
  {
    itkExceptionMacro("GetJacobianWithRespectToPosition "
                        "not yet implemented.");
  }

  /** Set the image where the transformation is defined. */
  void SetImage( ImageType * image);

  /** Get the Transformation Parameters. */
  const ParametersType& GetParameters(void) const;

  void SetParameters( const ParametersType & parameters );

  /** Print self */
  void PrintSelf(std::ostream &os, Indent indent) const
  {
    Superclass::PrintSelf(os,indent);
  }

  /** Get the transformation associated to a slice. */
  TransformType * GetSliceTransform( unsigned int i ) const
  {
    return m_TransformList[i];
  }

  /** Set the transformation parameters for a slice. */
  void SetSliceParameters( unsigned int i, const ParametersType & parameters )
  {
    m_TransformList[i] -> SetParameters( parameters );
    this -> Modified();
  }

  /** Initialize with the identity. */
  void Initialize();

  /** Initialize with a transformation. */
  void Initialize( TransformType * t );

  /** Set the Fixed Parameters. */
  void SetFixedParameters( const ParametersType & fp );

  /** Get the Fixed Parameters. */
  const ParametersType & GetFixedParameters(void) const;

  /** Get the number of slices (transforms). */
  itkGetMacro( NumberOfSlices, unsigned int);


protected:
    /** Default constructor. Otherwise we get a run time warning from itkTransform. */
  SliceBySliceTransform() : Superclass( NDimensions, 0 ) {}

private:
  /** List of transforms. */
  TransformPointerList m_TransformList;
  ImagePointerType 		 m_Image;
  unsigned int 				 m_NumberOfSlices;
  unsigned int 				 m_ParametersPerSlice;

  /** Temporary storage for transformation Jacobian. */
  mutable JacobianType m_Jacobian;
};

}

# include "btkSliceBySliceTransform.txx"

#endif /* __btkSliceBySliceTransform_h */
