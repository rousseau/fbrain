/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/08/2012
  Author(s): Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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

#ifndef BTK_SLICESINTERSECTIONITKCOSTFUNCTION_TXX
#define BTK_SLICESINTERSECTIONITKCOSTFUNCTION_TXX

#include "btkSlicesIntersectionITKCostFunction.hxx"

namespace btk
{
template< typename TImage >
SlicesIntersectionITKCostFunction<TImage>::SlicesIntersectionITKCostFunction()
    :m_VerboseMode(false),m_NumberOfParameters(6),m_MovingImageNum(0), m_MovingSliceNum(0)
{
    m_VNLCostFunction = NULL;
    //TODO : Add a constructor with parameters (m_NumberOfParameters)
    m_VNLCostFunction = new SlicesIntersectionVNLCostFunction<TImage>(m_NumberOfParameters);
}
//-------------------------------------------------------------------------------------------------
template< typename TImage>
SlicesIntersectionITKCostFunction<TImage>::~SlicesIntersectionITKCostFunction()
{
    if(m_VNLCostFunction !=NULL)
    {
        delete m_VNLCostFunction;
        m_VNLCostFunction = NULL;
    }
}
//-------------------------------------------------------------------------------------------------
template< typename TImage >
void SlicesIntersectionITKCostFunction<TImage>::Initialize()
{

        m_VNLCostFunction->SetImages(m_Images);
        m_VNLCostFunction->SetMasks(m_Masks);
        m_VNLCostFunction->SetTransforms(m_Transforms);
        m_VNLCostFunction->SetInverseTransforms(m_InverseTransforms);
        m_VNLCostFunction->SetVerboseMode(m_VerboseMode);

        m_VNLCostFunction->SetMovingImageNum(m_MovingImageNum);
        m_VNLCostFunction->SetMovingSliceNum(m_MovingSliceNum);
        //m_VNLCostFunction->SetCenterOfTransform(m_CenterOfTransform);

        m_VNLCostFunction->SetSlicesGroup(m_SlicesGroup);
        m_VNLCostFunction->SetGroupNum(m_GroupNum);

        m_VNLCostFunction->Initialize();

}

//-------------------------------------------------------------------------------------------------
template< typename TImage >
typename SlicesIntersectionITKCostFunction< TImage >::MeasureType
SlicesIntersectionITKCostFunction<TImage>::GetValue(const ParametersType &parameters) const
{
    MeasureType cost = m_VNLCostFunction->f(parameters);

    return ( cost );
}
//-------------------------------------------------------------------------------------------------
template< typename TImage >
void SlicesIntersectionITKCostFunction< TImage >::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
   derivative = m_VNLCostFunction->GetGradient(parameters);
}
//-------------------------------------------------------------------------------------------------
template< typename TImage >
unsigned int SlicesIntersectionITKCostFunction< TImage >::GetNumberOfParameters() const
{
    return m_NumberOfParameters;
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

}


#endif // BTKSLICESINTERSECTIONITKCOSTFUNCTION_TXX
