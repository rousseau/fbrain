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
  virtual OutputPointType TransformPoint(const InputPointType& p ) const
  {
    OutputPointType Tp( p );
    typename ImageType::IndexType index;
    m_Image -> TransformPhysicalPointToIndex( Tp , index);

    Tp = m_TransformList[ index[2] ] -> TransformPoint( Tp );

    return Tp;
  }

  /**  Method to transform a vector. */
  virtual OutputVectorType TransformVector(const InputVectorType& p) const
  {
    OutputVectorType Tp( p );
    for ( TransformPointerListConstIterator it = this->m_TransformList.begin(); it != this->m_TransformList.end(); ++it )
      {
      Tp = (*it)->TransformVector( Tp );
      }
    return Tp;
  }

  /**  Method to transform a vnl_vector. */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &p ) const
  {
    OutputVnlVectorType Tp( p );
    for ( TransformPointerListConstIterator it = this->m_TransformList.begin(); it != this->m_TransformList.end(); ++it )
      {
      Tp = (*it)->TransformVector( Tp );
      }
    return Tp;
  }

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType TransformCovariantVector( const InputCovariantVectorType &p ) const
  {
    OutputCovariantVectorType Tp( p );
    for ( TransformPointerListConstIterator it = this->m_TransformList.begin(); it != this->m_TransformList.end(); ++it )
      {
      Tp = (*it)->TransformCovariantVector( Tp );
      }
    return Tp;
  }

  /** Compute Jacobian (product of component Jacobians). Not thread safe!! */
  virtual const JacobianType & GetJacobian(const InputPointType &p ) const
  {

    this->m_Jacobian.SetSize( NDimensions, this->GetNumberOfParameters() );
    this->m_Jacobian.Fill(0.0);

    typename ImageType::IndexType index;
    m_Image -> TransformPhysicalPointToIndex( p , index);
    JacobianType jacobian = m_TransformList[ index[2] ] -> GetJacobian(p);

    unsigned int offset = index[2]*m_TransformList[0]->GetNumberOfParameters();

    for(unsigned int i = 0; i < NDimensions; i++)
      for(unsigned int j = 0; j < m_TransformList[0]->GetNumberOfParameters(); j++)
      {
        this->m_Jacobian[i][j] = jacobian[i][j+offset];
      }

    return this->m_Jacobian;

  }

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

  void SetImage( ImageType * image)
  {
    m_Image = image;
    typename ImageType::SizeType size = m_Image -> GetLargestPossibleRegion().GetSize();

    m_NumberOfSlices = size[2];
    m_TransformList.resize(m_NumberOfSlices);

    TransformPointer tmpTransform = TransformType::New();
    m_ParametersPerSlice = tmpTransform -> GetNumberOfParameters();
    this -> m_Parameters.SetSize( m_NumberOfSlices * m_ParametersPerSlice );

    InputPointType centerPoint;
    ContinuousIndexType centerIndex;

    centerIndex[0] = (size[0]-1)/2.0;
    centerIndex[1] = (size[1]-1)/2.0;

    for(unsigned int i=0; i<m_NumberOfSlices; i++)
    {
      centerIndex[2] =  i;
      m_Image -> TransformContinuousIndexToPhysicalPoint(centerIndex,centerPoint);

      m_TransformList[i] = TransformType::New();
      m_TransformList[i] -> SetIdentity();
      m_TransformList[i] -> SetCenter(centerPoint);
    }

  }

  /** Get the Transformation Parameters. */
  const ParametersType& GetParameters(void) const
  {
  for(unsigned int i=0; i<m_NumberOfSlices; i++)
  {
    ParametersType param = m_TransformList[i] -> GetParameters();

    for (unsigned int p=0; p<m_ParametersPerSlice; p++)
    {
      this->m_Parameters[p + i*m_ParametersPerSlice] = param[p];
    }

  }
  return this->m_Parameters;

  }

  void SetParameters( const ParametersType & parameters )
  {

    this -> m_Parameters = parameters;

    for(unsigned int i=0; i<m_NumberOfSlices; i++)
    {
      ParametersType param(m_ParametersPerSlice);

      for (unsigned int p=0; p<m_ParametersPerSlice; p++)
        param[p] = this->m_Parameters[p + i*m_ParametersPerSlice];

      m_TransformList[i] -> SetParameters( param );

    }
  }


  /** Print self */
  void PrintSelf(std::ostream &os, Indent indent) const
  {
    Superclass::PrintSelf(os,indent);
  }

  TransformType * GetSliceTransform( unsigned int i ) const
  {
    return m_TransformList[i];
  }


  void SetSliceParameters( unsigned int i, const ParametersType & parameters )
  {
    m_TransformList[i] -> SetParameters( parameters );
  }


  void Initialize( TransformType * t )
  {
    for(unsigned int i=0; i<m_NumberOfSlices; i++)
    {

      typename TransformType::Pointer inverseRigidTransform = TransformType::New();
      inverseRigidTransform -> SetIdentity();
      inverseRigidTransform -> SetCenter( t -> GetCenter() );
      inverseRigidTransform -> SetParameters( t -> GetParameters() );
      inverseRigidTransform -> GetInverse( inverseRigidTransform );
      inverseRigidTransform -> SetCenter( m_TransformList[i] -> GetCenter() );

      m_TransformList[i] -> SetParameters( inverseRigidTransform -> GetParameters() );
      m_TransformList[i] -> SetCenter( m_TransformList[i] -> GetCenter() );

    }
  }

  void SetFixedParameters( const ParametersType & fp )
  {
    this -> m_FixedParameters = fp;

    m_NumberOfSlices = this -> m_FixedParameters[0];

    if (m_TransformList.size() == 0 )
    {

      m_TransformList.resize(m_NumberOfSlices);

      TransformPointer tmpTransform = TransformType::New();
      m_ParametersPerSlice = tmpTransform -> GetNumberOfParameters();
      this -> m_Parameters.SetSize( m_NumberOfSlices * m_ParametersPerSlice );

      for(unsigned int i=0; i<m_NumberOfSlices; i++)
      {
        m_TransformList[i] = TransformType::New();
        m_TransformList[i] -> SetIdentity();
      }

    }

    InputPointType c;
    typedef typename ParametersType::ValueType ParameterValueType;

    for ( unsigned int j = 0; j<m_NumberOfSlices; j++)
    {
      for ( unsigned int i = 0; i < NDimensions; i++ )
        c[i] = this -> m_FixedParameters[j*NDimensions + i + 1];

      m_TransformList[j] -> SetCenter ( c );
    }

  }


  /** Get the Fixed Parameters. */
  const ParametersType & GetFixedParameters(void) const
  {
    this->m_FixedParameters.SetSize ( NDimensions * m_NumberOfSlices + 1 );

    this->m_FixedParameters[0] = m_NumberOfSlices;

    for ( unsigned int j = 0; j < m_NumberOfSlices; j++ )
    {
      for ( unsigned int i = 0; i < NDimensions; i++ )
        {
        this->m_FixedParameters[j*NDimensions + i + 1] = m_TransformList[j]-> GetCenter()[i];
        }
    }
    return this->m_FixedParameters;
  }


  itkGetMacro( NumberOfSlices, unsigned int);


protected:
    /** Default constructor. Otherwise we get a run time warning from itkTransform. */
  SliceBySliceTransform() : Superclass( NDimensions, 0 ) {}

private:
  /** List of transformations. */
  TransformPointerList m_TransformList;
  ImagePointerType m_Image;
  unsigned int m_NumberOfSlices;
  unsigned int m_ParametersPerSlice;


  /** Temporary storage for transformation Jacobian. */
  mutable JacobianType m_Jacobian;
};

}

//# include "btkSliceBySliceTransform.txx"

#endif /* __btkSliceBySliceTransform_h */
