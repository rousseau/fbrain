/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 31/01/2014
  Author(s): François Rousseau
  
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

#include "btkPandoraBoxTransform.h"

namespace btk
{
void PandoraBoxTransform::ConvertParametersToRigidMatrix(itkTransformType::Pointer outputTransform, std::vector<float> & inputParameters, std::vector<float> & center)
{
  //Quite ugly (and slow) way to do that ! (but we keep the flexibility of using MatrixOffsetTransformBase instead of specific transform such as Euler or Affine)
  //Create a ITK-based Euler transform, compute matrix, and then get parameters
  //inputParameters : 3 rotation angles (degrees) and 3 translations (mm)
  //Angles must be in degrees here

  // constant for converting degrees into radians
  const double dtr = ( vcl_atan(1.0) * 4.0 ) / 180.0;

  itk::Euler3DTransform<double>::Pointer eulerTransform = itk::Euler3DTransform<double>::New();
  eulerTransform->SetRotation(dtr*inputParameters[0], dtr*inputParameters[1], dtr*inputParameters[2]);
  outputTransform->SetMatrix( eulerTransform->GetMatrix() );

  itkTransformType::OutputVectorType translation;
  translation[0] = inputParameters[3];
  translation[1] = inputParameters[4];
  translation[2] = inputParameters[5];

  outputTransform->SetTranslation( translation );
  
  //Check center, offset, translation, affine matrix...
  //std::cout<<"Matrix : \n"<<outputTransform->GetMatrix()<<std::endl;
  //std::cout<<"Translation : "<<outputTransform->GetTranslation()<<std::endl;
  //std::cout<<"Center : "<<outputTransform->GetCenter()<<std::endl;
  //std::cout<<"Offset : "<<outputTransform->GetOffset()<<std::endl;

}

void PandoraBoxTransform::ConvertRigidMatrixToParameters(std::vector<float> & outputParameters, itkTransformType::Pointer inputTransform)
{
  //Same as above (ConvertParametersToRigidMatrix)
  //Use a ITK-based Euler transform to compute the parameters

  // constant for converting degrees into radians
  const double rtd = 180.0 / ( vcl_atan(1.0) * 4.0 ) ;

  itk::Euler3DTransform<double>::Pointer eulerTransform = itk::Euler3DTransform<double>::New();
  eulerTransform->SetMatrix( inputTransform->GetMatrix() );

  outputParameters[0] = eulerTransform->GetAngleX() * rtd;
  outputParameters[1] = eulerTransform->GetAngleY() * rtd;
  outputParameters[2] = eulerTransform->GetAngleZ() * rtd;

  itkTransformType::OutputVectorType translation;
  translation = inputTransform->GetTranslation();
  outputParameters[3] = translation[0];
  outputParameters[4] = translation[1];
  outputParameters[5] = translation[2];
}
  
void PandoraBoxTransform::ConvertParametersToMatrix(itkTransformType::Pointer outputTransform, vnl_vector< double > & inputParameters, itkVector & center)
{
  unsigned int n = inputParameters.size();
  
  //Compute rotation matrix using Euler decomposition
  vnl_matrix< double > rotation(3,3,0);
  const double dtr = ( vcl_atan(1.0) * 4.0 ) / 180.0;
  
  itk::Euler3DTransform<double>::Pointer eulerTransform = itk::Euler3DTransform<double>::New();
  eulerTransform->SetRotation(dtr*inputParameters[0], dtr*inputParameters[1], dtr*inputParameters[2]);
  rotation = eulerTransform->GetMatrix().GetVnlMatrix();
  
  //Compute scale matrix
  vnl_matrix< double > scale(3,3,0);
  if(inputParameters.size() == 6)
  {
    scale.set_identity();
  }
  else
  {
    if(n >= 7)
    {
      scale[0][0] = inputParameters[6];
      if(n >= 8)
        scale[1][1] = inputParameters[7];
      else
        scale[1][1] = inputParameters[6];
      if(n >= 9)
        scale[2][2] = inputParameters[8];
      else
        scale[2][2] = inputParameters[6];
    }
  }
  
  //Compute skew matrix
  vnl_matrix< double > skew(3,3,0);
  skew.set_identity();
  if (n >= 10)
  {
    skew[0][1] = inputParameters[9];
    if (n >= 11)
      skew[0][2] = inputParameters[10];
    if (n >= 12)
      skew[1][2]= inputParameters[11];
  }
  
  //Update transform matrix
  outputTransform->SetMatrix(rotation*skew*scale);
  
  //std::cout<<outputTransform->GetMatrix()<<std::endl;
  
  //Update translation
  itkVector translation;
  translation[0] = inputParameters[3];
  translation[1] = inputParameters[4];
  translation[2] = inputParameters[5];
  outputTransform->SetTranslation(translation);
  
  //Update center
  outputTransform->SetCenter(center);
  
  //Check center, offset, translation, affine matrix...
  //std::cout<<"ConvertParametersToMatrix"<<std::endl;
  //std::cout<<"Matrix : \n"<<outputTransform->GetMatrix()<<std::endl;
  //std::cout<<"Translation : "<<outputTransform->GetTranslation()<<std::endl;
  //std::cout<<"Center : "<<outputTransform->GetCenter()<<std::endl;
  //std::cout<<"Offset : "<<outputTransform->GetOffset()<<std::endl;

}
  
  

} // namespace btk
