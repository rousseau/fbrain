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

#ifndef BTK_PANDORA_BOX_RECONSTRUCTION_FILTERS_H
#define BTK_PANDORA_BOX_RECONSTRUCTION_FILTERS_H

// STL includes
#include "vector"

#define _USE_MATH_DEFINES
#include <math.h>

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransform.h"
#include "itkContinuousIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolationWeightFunction.h"
#include "itkSubtractImageFilter.h"

#include "vnl/vnl_sparse_matrix.h"

namespace btk
{

class PandoraBoxReconstructionFilters
{
    public:
    typedef itk::Image< float, 3>                              itkFloatImage;
    typedef itkFloatImage::Pointer                             itkFloatImagePointer;
    typedef itk::ImageRegionIterator< itkFloatImage >          itkFloatIterator;
    typedef itk::ImageRegionIteratorWithIndex< itkFloatImage > itkFloatIteratorWithIndex;
    typedef itk::ContinuousIndex<double,3>                     itkContinuousIndex;

    typedef itk::MatrixOffsetTransformBase<double,3,3>         itkTransformType;

    //Convert 3D image to a stack of 3D images (slices)
    static void Convert3DImageToSliceStack(std::vector<itkFloatImagePointer> & outputStack, itkFloatImagePointer & inputImage);

    //Convert a stack of 3D images to a 3D image
    static void ConvertSliceStackTo3DImage(itkFloatImagePointer & outputImage, std::vector<itkFloatImagePointer> & inputStack);

    //Project a 3D image into a stack of slices
    static void Project3DImageToSliceStack(std::vector<itkFloatImagePointer> & outputStack, itkFloatImagePointer & inputImage, std::vector<itkFloatImagePointer> & inputStack, std::vector<itkTransformType::Pointer> & affineSBSTransforms);

    //Compute PSF
    static void ComputePSFImage(itkFloatImagePointer & PSFImage, itkFloatImage::SpacingType HRSpacing, itkFloatImage::SpacingType LRSpacing);

    //Compute parameters of the observation model : H, Y, X (Y = HX)
    static void ComputerObservationModelParameters(vnl_sparse_matrix<float> & H, vnl_vector<float> & Y, vnl_vector<float> & X, itkFloatImagePointer & HRImage, std::vector< std::vector<itkFloatImagePointer> > & maskStacks, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkTransformType::Pointer> > & inverseAffineSBSTransforms, itkFloatImagePointer & PSFImage);

    //Injection
    static void ImageFusionByInjection(itkFloatImagePointer & outputImage, itkFloatImagePointer & maskImage, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkTransformType::Pointer> > & affineSBSTransforms);

    //Simulate observations using the observation model Y=HX
    static void SimulateObservations(vnl_sparse_matrix<float> & H, vnl_vector<float> & X, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkFloatImagePointer> > & outputStacks);

    //Convert a 3D image to a vnl vector (HRimage -> X)
    static void Convert3DImageToVNLVector(vnl_vector<float> & X, itkFloatImagePointer & inputImage);

    //Convert a vnl vector to a 3D image (X -> HRimage)
    static void ConvertVNLVectorTo3DImage(itkFloatImagePointer & outputImage, vnl_vector<float> & X);

    //Convert a vnl vector into a image stack (Y -> stacks)
    static void ConvertVNLVectorToSliceStack(std::vector< std::vector<itkFloatImagePointer> > & outputStacks, vnl_vector<float> & Y);

    //Iterative back projection
    static void IterativeBackProjection(itkFloatImagePointer & outputImage, itkFloatImage::SpacingType & outputSpacing, std::vector< std::vector<itkFloatImagePointer> > & inputStacks, std::vector< std::vector<itkFloatImagePointer> > & maskStacks, std::vector< std::vector<itkTransformType::Pointer> > & affineSBSTransforms, std::vector< std::vector<itkTransformType::Pointer> > & inverseAffineSBSTransforms, unsigned int maxIterations);

    //SR (L1,L2,robust) + Reg(local, patch, tv)

    //Slice motion estimation using a 3D reference image

    //Slice motion estimation using a observation model

    //Slice motion estimation using a set of slices

    //protected:

    //private:

};

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkPandoraBoxReconstructionFilters.txx"
#endif

#endif // BTK_PANDORA_BOX_RECONSTRUCTION_FILTERS_H
