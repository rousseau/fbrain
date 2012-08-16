/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 13/08/2012
  Author(s): Julien Pontabry (pontabry@unistra.fr)
  
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

#include "btkTensorModel.h"

namespace btk
{

TensorModel::TensorModel() : m_BValue(1500)
{
    // ----
}

//----------------------------------------------------------------------------------------

TensorModel::~TensorModel()
{
    // ----
}

//----------------------------------------------------------------------------------------

void TensorModel::Update()
{
    // Create interpolate function on model image
    m_ModelImageFunction = InterpolateModelFunction::New();
    m_ModelImageFunction->SetInputImage(m_InputModelImage);
}

//----------------------------------------------------------------------------------------

float TensorModel::ModelAt(ContinuousIndex cindex, GradientDirection direction)
{
    ModelImage::PixelType tensor = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    ModelImage::PixelType::EigenValuesArrayType   eigenValues;
    ModelImage::PixelType::EigenVectorsMatrixType eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);

    float coeffx = (eigenValues[2] > 0.f) ? std::sqrt( 2.0 * eigenValues[2] ) : 0.f; // FIXME ?
    float coeffy = (eigenValues[1] > 0.f) ? std::sqrt( 2.0 * eigenValues[1] ) : 0.f; // FIXME ?
    float coeffz = (eigenValues[0] > 0.f) ? std::sqrt( 2.0 * eigenValues[0] ) : 0.f; // FIXME ?

    // tmp = uTEe
    float x = coeffx * direction[0];
    float y = coeffy * direction[1];
    float z = coeffz * direction[2];
    btk::GradientDirection tmp(
                    x*eigenVectors(2,0) + y*eigenVectors(1,0) + z*eigenVectors(0,0),
                    x*eigenVectors(2,1) + y*eigenVectors(1,1) + z*eigenVectors(0,1),
                    x*eigenVectors(2,2) + y*eigenVectors(1,2) + z*eigenVectors(0,2)
                );

    // Return the rau parameter of corresponding spherical harmonics direction
    return std::sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
}

//----------------------------------------------------------------------------------------

float TensorModel::ModelAt(PhysicalPoint point, GradientDirection direction)
{
    ContinuousIndex cindex;
    m_InputModelImage->TransformPhysicalPointToContinuousIndex(point, cindex);

    return ModelAt(cindex, direction);
}

//----------------------------------------------------------------------------------------

float TensorModel::SignalAt(ContinuousIndex cindex, GradientDirection direction)
{
    ModelImage::PixelType tensor = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    // s = exp[ -b uTDu ]
    return std::exp( -static_cast< float >(m_BValue) * (
                            (direction[0]*tensor(0,0) + direction[1]*tensor(1,0) + direction[2]*tensor(2,0))*direction[0] +
                            (direction[0]*tensor(0,1) + direction[1]*tensor(1,1) + direction[2]*tensor(2,1))*direction[1] +
                            (direction[0]*tensor(0,2) + direction[1]*tensor(1,2) + direction[2]*tensor(2,2))*direction[2]
                         ) );
}

//----------------------------------------------------------------------------------------

float TensorModel::SignalAt(PhysicalPoint point, GradientDirection direction)
{
    ContinuousIndex cindex;
    m_InputModelImage->TransformPhysicalPointToContinuousIndex(point, cindex);

    return SignalAt(cindex, direction);
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > TensorModel::MeanDirectionsAt(ContinuousIndex cindex)
{
    ModelImage::PixelType tensor = m_ModelImageFunction->EvaluateAtContinuousIndex(cindex);

    ModelImage::PixelType::EigenValuesArrayType   eigenValues;
    ModelImage::PixelType::EigenVectorsMatrixType eigenVectors;
    tensor.ComputeEigenAnalysis(eigenValues, eigenVectors);

    std::vector< btk::GradientDirection > meanDirections;
    meanDirections.push_back(btk::GradientDirection(eigenVectors(2,0), eigenVectors(2,1), eigenVectors(2,2)));

    return meanDirections;
}

//----------------------------------------------------------------------------------------

std::vector< btk::GradientDirection > TensorModel::MeanDirectionsAt(PhysicalPoint point)
{
    ContinuousIndex cindex;
    m_InputModelImage->TransformPhysicalPointToContinuousIndex(point, cindex);

    return MeanDirectionsAt(cindex);
}

//----------------------------------------------------------------------------------------

void TensorModel::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    Superclass::PrintSelf(os, indent);
}

} // namespace btk
