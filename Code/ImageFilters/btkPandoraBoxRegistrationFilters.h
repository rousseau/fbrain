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

#ifndef BTK_PANDORA_BOX_REGISTRATION_FILTERS_H
#define BTK_PANDORA_BOX_REGISTRATION_FILTERS_H

// STL includes
#include "vector"
#include <string>

#define _USE_MATH_DEFINES
#include <math.h>
//#include <random> //unfortunately, there are some compilation issues with ITK and C++11
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkMatrixOffsetTransformBase.h"
//#include "itkTransform.h"
#include "itkContinuousIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolationWeightFunction.h"
//#include "itkEuler3DTransform.h"
#include "itkPointSet.h"

#include "vnl/vnl_sparse_matrix.h"

#include "btkPandoraBoxCostFunction.h"
#include "btkSimplex.h"


namespace btk
{

class PandoraBoxRegistrationFilters
{
    public:
    typedef itk::Image< float, 3>                              itkFloatImage;
    typedef itkFloatImage::Pointer                             itkFloatImagePointer;
    typedef itk::ImageRegionIterator< itkFloatImage >          itkFloatIterator;
    typedef itk::ImageRegionIteratorWithIndex< itkFloatImage > itkFloatIteratorWithIndex;
    typedef itk::ContinuousIndex<double,3>                     itkContinuousIndex;

    typedef itk::Vector<float, 1>                              itkVector;
    typedef itk::Image<itkVector, 3>                           itkVectorImage;
    typedef itk::PointSet<itkVectorImage::PixelType, 3>        itkPointSet;
    typedef itk::ImageRegionIterator< itkVectorImage >         itkVectorIterator;


    //typedef itk::MatrixOffsetTransformBase<double,3,3>         itkTransformType;

    //Common functions for slice motion estimation
    //By default, the center of the transforms is 0,0,0 in the physical space.
    //It has to be checked whether this choice is good enough in our case, or whether we should move it to the center of the slice

    //Slice motion estimation using a 3D reference image
    //Register3DImages is the basic function, used in any affine registration approach
    static void Register3DImages(itkFloatImagePointer & movingImage, itkFloatImagePointer & movingMaskImage, itkFloatImagePointer & referenceImage, itkFloatImagePointer & referenceMaskImage, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, vnl_vector< double > toleranceVector, double tolerance);
    static void Register3DImages(PandoraBoxCostFunction & costFunction, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, vnl_vector< double > toleranceVector, double tolerance);

    //MultiStart3DRegistration implements a multi-start strategy by perturbing the input parameters
    static void GenerateRandomParameters(vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange);
    static void GenerateUniformlyDistributedParameters(vnl_vector< double > & inputParameters, std::vector< vnl_vector< double > > & outputParameters, vnl_vector< double > & parameterRange, int samplingRate);

    static void MultiStart3DRegistration(PandoraBoxCostFunction & costFunction, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, unsigned int numberOfPerturbations, vnl_vector< double > toleranceVector, double tolerance);

    //Coarse3DRegistration explores the domain (i.e. possible rotation values) of registration parameters
    static void Coarse3DRegistration(PandoraBoxCostFunction & costFunction, vnl_vector< double > & inputParameters, vnl_vector< double > & outputParameters, vnl_vector< double > & parameterRange, unsigned int numberOfIncrements, vnl_vector< double > toleranceVector, double tolerance);

    //Slice motion estimation using a observation model

    //Slice motion estimation using a set of slices

    //Outliers detection in each slice

    //Slice-based relative bias correction


    //protected:

    //private:

};

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkPandoraBoxRegistrationFilters.txx"
#endif

#endif // BTK_PANDORA_BOX_REGISTRATION_FILTERS_H
