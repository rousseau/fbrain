/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 20/03/2013
  Author(s):Marc Schweitzer (marc.schweitzer(at)unistra.fr)
  
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

#ifndef BTK_RANDOMSLICEBYSLICETRANSFORMGENARATOR_H
#define BTK_RANDOMSLICEBYSLICETRANSFORMGENARATOR_H

#include "itkObject.h"
#include "itkImage.h"
#include "itkEuler3DTransform.h"

#include "btkEulerSliceBySliceTransform.h"
#include "btkMacro.h"


namespace btk
{

class RandomSliceBySliceTransformGenerator: public itk::Object
{
    public:
        typedef RandomSliceBySliceTransformGenerator Self;
        typedef itk::Object                          Superclass;
        typedef itk::SmartPointer< Self >            Pointer;
        typedef itk::SmartPointer< const Self >      ConstPointer;
        typedef btk::EulerSliceBySliceTransform< double, 3 > TransformType;
        typedef itk::Euler3DTransform< double > Rigid3DTransformType;
        typedef Rigid3DTransformType::ParametersType   Parameters;
        typedef itk::Image< float, 3> ImageType; // type of image is not very important here, we  set it to float to avoid useless template

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(RandomSliceBySliceTransformGenarator, itk::Object);

        void Update();

        btkSetMacro(MaxRotation, float);
        btkSetMacro(MaxTranslation,float);
        btkSetMacro(Level, unsigned int);

        btkGetMacro(Transform,TransformType::Pointer);


        btkSetMacro(ActiveParameters,Parameters);

        btkSetMacro(Image, ImageType::Pointer);

        btkSetMacro(VerboseMode,bool);
        btkGetMacro(VerboseMode,bool);


    protected:
        RandomSliceBySliceTransformGenerator();
        virtual ~RandomSliceBySliceTransformGenerator(){};

    private:

        TransformType::Pointer  m_Transform;
        unsigned int            m_NumberOfSlices;
        float                   m_MaxRotation;
        float                   m_MaxTranslation;
        float                   m_MinRotation;
        float                   m_MinTranslation;
        unsigned int            m_Level;
        Parameters              m_ActiveParameters;
        unsigned int            m_NumberOfParameters;
        ImageType::Pointer      m_Image; /** Image to be transformed, slice by slice transform use image inside */

        bool                    m_VerboseMode;


};

} // namespace btk

#endif // BTK_RANDOMSLICEBYSLICETRANSFORMGENARATOR_H
