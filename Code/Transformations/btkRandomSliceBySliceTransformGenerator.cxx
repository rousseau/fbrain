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

#include "btkRandomSliceBySliceTransformGenerator.h"
#include "btkMathFunctions.h"

namespace btk
{

RandomSliceBySliceTransformGenerator::RandomSliceBySliceTransformGenerator():m_NumberOfParameters(6),m_MaxRotation(5),m_MaxTranslation(5),
    m_VerboseMode(false)
{
    m_ActiveParameters.set_size(m_NumberOfParameters);
    m_ActiveParameters.Fill(1);
}

void RandomSliceBySliceTransformGenerator::Update()
{
    if(m_Image.IsNull())
    {
        return;
    }

    std::cout<<"Compute a random slice by slice transform..."<<std::endl;


    m_NumberOfSlices = m_Image->GetLargestPossibleRegion().GetSize()[2];

    std::cout<<"Number of Slices : "<<m_NumberOfSlices<<std::endl;
    std::cout<<"Processing..."<<std::flush;

    m_MaxRotation = btk::MathFunctions::DegreesToRadians(m_MaxRotation);

    double m_MinRotation =  -m_MaxRotation;
    double m_MinTranslation =  -m_MaxTranslation;

   int threshold = 0;
   switch(m_Level)
   {
   case 1:
       threshold = 8;
       break;

   case 2:
       threshold = 5;
       break;
   case 3:
       threshold = 3;
       break;

   case 4:
       threshold = 0;//all slice
       break;

   default:
       threshold = 8;
       break;

   }


    m_Transform = TransformType::New();
    m_Transform->SetImage(m_Image); // To compute the right number of slice and the right center
    m_Transform->Initialize();


    //Slice N°i random transformation :
    Rigid3DTransformType::Pointer RandomTransform = Rigid3DTransformType::New();
    srand(time(NULL));


    Rigid3DTransformType::ParametersType randomParams(RandomTransform->GetNumberOfParameters());


    for(int i = 0; i<m_NumberOfSlices; i++)
    {
        //randomParams = m_Transform->GetSliceParameters(i);
        randomParams.Fill(0);
        m_Transform->SetSliceParameters(i, randomParams);

        int rndm = rand()%10;
        RandomTransform->SetIdentity();
        Rigid3DTransformType::ParametersType randomParams(RandomTransform->GetNumberOfParameters());

        if(rndm >= threshold)
        {
            int rot_trans = rand()%3;

            //Rotation
            if(rot_trans == 0)
            {
                int paramR = rand()%3 ; // x y and z
                if(m_ActiveParameters[paramR] == 1)
                {
                    randomParams[paramR] = rand()/((double)RAND_MAX/std::abs(m_MaxRotation - m_MinRotation )) - std::abs(m_MinRotation);
                }

            }
            //Translation
            else if(rot_trans == 1)
            {

                int paramT = (rand()%3) + 3;
                if(m_ActiveParameters[paramT] == 1)
                {
                    randomParams[paramT] = rand()/((double)RAND_MAX/ std::abs(m_MaxTranslation - m_MinTranslation)) - std::abs(m_MinTranslation);
                }

            }
            //Rotation and translation
            else
            {
                int paramR = rand()%3;
                int paramT = (rand()%3) + 3;
                if(m_ActiveParameters[paramR] == 1)
                {
                    randomParams[paramR] = rand()/((double)RAND_MAX/std::abs(m_MaxRotation - m_MinRotation )) - std::abs(m_MinRotation);
                }
                if(m_ActiveParameters[paramT] == 1)
                {
                    randomParams[paramT] = rand()/((double)RAND_MAX/ std::abs(m_MaxTranslation - m_MinTranslation)) - std::abs(m_MinTranslation);
                }

            }


            RandomTransform->SetParameters(randomParams);
        }

        m_Transform->SetSliceParameters(i, RandomTransform->GetParameters());

        if(m_VerboseMode)
        {
            for(int x=0;x<m_NumberOfParameters;x++)
            {
                std::cout<<"Parameters "<<x<<" = "<<randomParams[x]<<std::endl;
            }
        }


        randomParams.Fill(0);

    }

    std::cout<<" Done !"<<std::endl;

}

} // namespace btk
