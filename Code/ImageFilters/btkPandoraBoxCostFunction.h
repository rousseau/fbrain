/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 12/12/2014
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

#ifndef BTK_PANDORA_BOX_COST_FUNCTION_H
#define BTK_PANDORA_BOX_COST_FUNCTION_H

// STL includes
#include "vector"

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkBSplineInterpolationWeightFunction.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkEuler3DTransform.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "btkPandoraBoxTransform.h"
#include "btkJointHistogram.h"


namespace btk
{
  
  class PandoraBoxCostFunction
  {
  public:
    typedef itk::Image< short, 3>                              itkShortImage;
    typedef itkShortImage::Pointer                             itkShortImagePointer;
    typedef itk::ImageRegionIterator< itkShortImage >          itkShortIterator;
    
    typedef itk::Image< float, 3>                              itkFloatImage;
    typedef itkFloatImage::Pointer                             itkFloatImagePointer;
    typedef itk::ImageRegionIterator< itkFloatImage >          itkFloatIterator;
    typedef itk::ImageRegionIteratorWithIndex< itkFloatImage > itkFloatIteratorWithIndex;

    
    //typedef itk::IdentityTransform<double, 3>                                    itkIdentityTransform;
    //typedef itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>  itkBSplineInterpolator;
    
    typedef itk::MatrixOffsetTransformBase<double,3,3>               itkTransformType;
    //itk::TransformFactory<itkTransformType>::RegisterTransform();
    
    typedef itk::ContinuousIndex<double,3>     itkContinuousIndex;
    typedef itk::MinimumMaximumImageCalculator <itkFloatImage>       itkImageCalculatorFilter;

    itkFloatImagePointer movingImage;
    itkFloatImagePointer referenceImage;
    itkFloatImagePointer movingMask;
    itkFloatImagePointer referenceMask;
    itk::Vector<double, 3> center;
    
    //Use two interpolators can slow down the registration but it is more rigourous to take into account both (image and mask)
    itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>::Pointer bsInterpolatorMovingImage;
    itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>::Pointer bsInterpolatorMovingMask;

    itkTransformType::Pointer transform;

    btk::JointHistogram jointHistogram;

    PandoraBoxCostFunction()
    {
      this->center[0] = 0;
      this->center[1] = 0;
      this->center[2] = 0;
    }
    
    void SetReferenceImage(itkFloatImagePointer & inputImage){
      referenceImage = inputImage;
      
      transform = itkTransformType::New();
    }
    void SetReferenceMask(itkFloatImagePointer & inputImage){
      referenceMask = inputImage;
    }
    
    void SetMovingImage(itkFloatImagePointer & inputImage){
      movingImage = inputImage;
      
      //Use currently a linear interpolation in the moving image
      bsInterpolatorMovingImage = itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>::New();
      bsInterpolatorMovingImage->SetSplineOrder(1);
      bsInterpolatorMovingImage->SetInputImage(movingImage);
    }
    
    void SetMovingMask(itkFloatImagePointer & inputImage){
      movingMask = inputImage;
      
      //Use nearest neighbour in the moving mask
      bsInterpolatorMovingMask = itk::BSplineInterpolateImageFunction<itkFloatImage, double, double>::New();
      bsInterpolatorMovingMask->SetSplineOrder(0);
      bsInterpolatorMovingMask->SetInputImage(movingMask);
    }
    
    void SetCenter(itk::Vector<double, 3> & c)
    {
      center[0] = c[0];
      center[1] = c[1];
      center[2] = c[2];
    }
    
    void InitializeJointHistogram(unsigned int nx, unsigned int ny){

        jointHistogram.SetNumberOfBins(nx,ny);

        //std::cout<<"Initialization of joint histogram (ax,bx,ay,by)"<<std::endl;
        itkImageCalculatorFilter::Pointer refCalculator = itkImageCalculatorFilter::New();
        refCalculator->SetImage(referenceImage);
        refCalculator->Compute();
        itkImageCalculatorFilter::Pointer movCalculator = itkImageCalculatorFilter::New();
        movCalculator->SetImage(referenceImage);
        movCalculator->Compute();

        jointHistogram.SetAx( (jointHistogram.GetNumberOfBinsX() - 1) * 1.0 / (refCalculator->GetMaximum() - refCalculator->GetMinimum()) );
        jointHistogram.SetAy( (jointHistogram.GetNumberOfBinsY() - 1) * 1.0 / (movCalculator->GetMaximum() - movCalculator->GetMinimum()) );
        jointHistogram.SetBx( - jointHistogram.GetAx() * refCalculator->GetMinimum() );
        jointHistogram.SetBy( - jointHistogram.GetAy() * movCalculator->GetMinimum() );

    }

    void FillJointHistogram(vnl_vector<double>  params){

        jointHistogram.ClearJointHistogram();

        btk::PandoraBoxTransform::ConvertParametersToMatrix(transform, params, this->center);

        itkFloatIteratorWithIndex itReference(referenceImage,referenceImage->GetLargestPossibleRegion());
        itkFloatIteratorWithIndex itMask(referenceMask,referenceMask->GetLargestPossibleRegion());
        itkFloatImage::IndexType refIndex; //index of the current voxel in the reference image
        itkFloatImage::PointType refPoint; //physical point location of reference image

        itkFloatImage::PointType transformedPoint; //Physical point location after applying transform
        itkContinuousIndex       inputContIndex;   //continuous index in the 3D image

        double res = 0;
        double weightedSum=0;
        double weight = 0;

        //loop over slice voxels
        for(itReference.GoToBegin(), itMask.GoToBegin(); !itReference.IsAtEnd(); ++itReference, ++itMask)
        {
          if(itMask.Get() > 0)
          {
            //Coordinate (index) of the current pixel in the current image
            refIndex = itReference.GetIndex();

            //Coordinate in the physical world (mm)
            referenceImage->TransformIndexToPhysicalPoint(refIndex,refPoint);

            //Apply affine transform
            transformedPoint = transform->TransformPoint(refPoint);

            //Coordinate in the 3D image (continuous index)
            movingImage->TransformPhysicalPointToContinuousIndex(transformedPoint,inputContIndex);

            //Should we symmetrized the criterion?
            //Take max(w1,w2) instead of w1, or w1*w2

            //Simple version, taking only into account for reference mask -------------------------------------------
            if(bsInterpolatorMovingImage->IsInsideBuffer(inputContIndex))
            {
              double interpolatedValueImage = bsInterpolatorMovingImage->EvaluateAtContinuousIndex(inputContIndex);
              weight = itMask.Get();
              //weightedSum += weight;
              //res += weight * pow(interpolatedValueImage - itReference.Get(), 2.0);
              jointHistogram.AddSample( itReference.Get(), interpolatedValueImage, weight );
            }
          }
        }

    }

    void FillJointHistogramOpenMP(vnl_vector<double>  params){

        jointHistogram.ClearJointHistogram();
        vnl_matrix< double >  data = jointHistogram.GetData();
        btk::PandoraBoxTransform::ConvertParametersToMatrix(transform, params, this->center);

        //solution : créer un vecteur de joint histogram
        //remplir chaque joint histogram independamment
        //faire la somme et l'affecter au membre joint histogram
        std::vector< btk::JointHistogram > jhVector;
        itkFloatImage::SizeType    size    = referenceImage->GetLargestPossibleRegion().GetSize();
        int x,y,z;
        double weightedSum=0;

        #pragma omp parallel
        {
            const int nthreads = omp_get_num_threads();
            const int ithread = omp_get_thread_num();

            #pragma omp single
            jhVector.resize(nthreads);

            #pragma omp single
            for(int i=0; i<nthreads; i++){
                jhVector[i].SetNumberOfBins(jointHistogram.GetNumberOfBinsX(),jointHistogram.GetNumberOfBinsY());
                jhVector[i].SetAx(jointHistogram.GetAx());
                jhVector[i].SetAy(jointHistogram.GetAy());
                jhVector[i].SetBx(jointHistogram.GetBx());
                jhVector[i].SetBy(jointHistogram.GetBy());
            }

            #pragma omp for private(x,y,z)
            for(z=0; z < (int)size[2]; z++)
            for(y=0; y < (int)size[1]; y++)
            for(x=0; x < (int)size[0]; x++)
            {
              itkFloatImage::IndexType refIndex; //index of the current voxel in the reference image
              refIndex[0] = x;
              refIndex[1] = y;
              refIndex[2] = z;

              if(referenceMask->GetPixel( refIndex ) > 0)
              {
                  itkFloatImage::PointType refPoint; //physical point location of reference image
                  referenceImage->TransformIndexToPhysicalPoint(refIndex,refPoint);

                  itkFloatImage::PointType transformedPoint; //Physical point location after applying transform
                  transformedPoint = transform->TransformPoint(refPoint);

                  itkContinuousIndex       inputContIndex;   //continuous index in the 3D image
                  movingImage->TransformPhysicalPointToContinuousIndex(transformedPoint,inputContIndex);

                  if(bsInterpolatorMovingImage->IsInsideBuffer(inputContIndex))
                  {
                    double interpolatedValueImage = bsInterpolatorMovingImage->EvaluateAtContinuousIndex(inputContIndex);
                    double weight = referenceMask->GetPixel( refIndex );

                    double currentReferenceValue = referenceImage->GetPixel( refIndex );

                    jhVector[ithread].AddSample(currentReferenceValue, interpolatedValueImage, weight);
                  }
              }
            }

            //Put everything together to get the final joint histogram
            #pragma omp for reduction(+:weightedSum)
            for(int i=0; i<nthreads; i++)
                weightedSum += jhVector[i].GetNumberOfSamples();

            //Adding all joint histograms can be done with or without openmp
            //#pragma omp single
            //for(int k=0; k<nthreads; k++)
            //    jointHistogram.SetData( jointHistogram.GetData() + jhVector[k].GetData() );

            #pragma omp for
            for(int i=0; i<jointHistogram.GetNumberOfBinsX(); i++)
            {
                for(int j=0; j<jointHistogram.GetNumberOfBinsY(); j++)
                    for(int k=0; k<nthreads; k++)
                        jointHistogram.m_Data(i,j) = jointHistogram.m_Data(i,j) + jhVector[k].m_Data(i,j);

            }

        }

        jointHistogram.SetNumberOfSamples(weightedSum);

    }

    virtual double operator () (vnl_vector<double> params)=0;
    //protected:
    
    //private:
    
  };
  
  class PandoraBoxCostFunctionMI : public PandoraBoxCostFunction
  {
  public:

    virtual double operator () (vnl_vector<double>  params)
    {
      this->FillJointHistogram(params);
      return - this->jointHistogram.MutualInformation();
    }
  };

  class PandoraBoxCostFunctionMIOpenMP : public PandoraBoxCostFunction
  {
  public:

    virtual double operator () (vnl_vector<double>  params)
    {
      this->FillJointHistogramOpenMP(params);
      return - this->jointHistogram.MutualInformation();
    }
  };

  class PandoraBoxCostFunctionNMI : public PandoraBoxCostFunction
  {
  public:

    virtual double operator () (vnl_vector<double>  params)
    {
      this->FillJointHistogram(params);
      return - this->jointHistogram.NormalizedMutualInformation();
    }
  };

  class PandoraBoxCostFunctionNMIOpenMP : public PandoraBoxCostFunction
  {
  public:

    virtual double operator () (vnl_vector<double>  params)
    {
      this->FillJointHistogramOpenMP(params);
      return - this->jointHistogram.NormalizedMutualInformation();
    }
  };

  class PandoraBoxCostFunctionMSE : public PandoraBoxCostFunction
  {
  public:
    
    
    virtual double operator () (vnl_vector<double>  params)
    {
      btk::PandoraBoxTransform::ConvertParametersToMatrix(transform, params, this->center);
            
      itkFloatIteratorWithIndex itReference(referenceImage,referenceImage->GetLargestPossibleRegion());
      itkFloatIteratorWithIndex itMask(referenceMask,referenceMask->GetLargestPossibleRegion());
      itkFloatImage::IndexType refIndex; //index of the current voxel in the reference image
      itkFloatImage::PointType refPoint; //physical point location of reference image
      
      itkFloatImage::PointType transformedPoint; //Physical point location after applying transform
      itkContinuousIndex       inputContIndex;   //continuous index in the 3D image

      double res = 0;
      double weightedSum=0;
      double weight = 0;

      //loop over slice voxels
      for(itReference.GoToBegin(), itMask.GoToBegin(); !itReference.IsAtEnd(); ++itReference, ++itMask)
      {
        if(itMask.Get() > 0)
        {
          //Coordinate (index) of the current pixel in the current image
          refIndex = itReference.GetIndex();
        
          //Coordinate in the physical world (mm)
          referenceImage->TransformIndexToPhysicalPoint(refIndex,refPoint);
        
          //Apply affine transform
          transformedPoint = transform->TransformPoint(refPoint);
        
          //Coordinate in the 3D image (continuous index)
          movingImage->TransformPhysicalPointToContinuousIndex(transformedPoint,inputContIndex);
        
          //Should we symmetrized the criterion?
          //Take max(w1,w2) instead of w1, or w1*w2

          //Simple version, taking only into account for reference mask -------------------------------------------
          if(bsInterpolatorMovingImage->IsInsideBuffer(inputContIndex))
          {
            double interpolatedValueImage = bsInterpolatorMovingImage->EvaluateAtContinuousIndex(inputContIndex);
            weight = itMask.Get();
            weightedSum += weight;
            res += weight * pow(interpolatedValueImage - itReference.Get(), 2.0);
          }

          //max(wref,wmov)-----------------------------------------------------------------------------------------
          /*
          double interpolatedValueMask = bsInterpolatorMovingMask->EvaluateAtContinuousIndex(inputContIndex);
          double interpolatedValueImage = bsInterpolatorMovingImage->EvaluateAtContinuousIndex(inputContIndex);
          weight = itMask.Get();
          if(weight < interpolatedValueMask)
              weight = interpolatedValueMask;
          weightedSum += weight;
          res += weight * pow(interpolatedValueImage - itReference.Get(), 2.0);
            */

          //Compute only the cost function if masks intersect
          /*
          double interpolatedValueMask = bsInterpolatorMovingMask->EvaluateAtContinuousIndex(inputContIndex);
          if(interpolatedValueMask > 0)
          {
            double interpolatedValueImage = bsInterpolatorMovingImage->EvaluateAtContinuousIndex(inputContIndex);
            //Weighted MSE using image mask value

            //weight = itMask.Get() * interpolatedValueMask; // very bad option ... this is very unstable
            weight = itMask.Get();
            weightedSum += weight;
            res += weight * pow(interpolatedValueImage - itReference.Get(), 2.0);
          }
          */
        }
      }
      return res / weightedSum;
    }
  };
    
  class PandoraBoxCostFunctionMSEOpenMP : public PandoraBoxCostFunction
  {
  public:


    virtual double operator () (vnl_vector<double>  params)
    {
      btk::PandoraBoxTransform::ConvertParametersToMatrix(transform, params, this->center);

//      itkFloatIteratorWithIndex itReference(referenceImage,referenceImage->GetLargestPossibleRegion());
//      itkFloatIteratorWithIndex itMask(referenceMask,referenceMask->GetLargestPossibleRegion());
//      itkFloatImage::IndexType refIndex; //index of the current voxel in the reference image
//      itkFloatImage::PointType refPoint; //physical point location of reference image

//      itkFloatImage::PointType transformedPoint; //Physical point location after applying transform
//      itkContinuousIndex       inputContIndex;   //continuous index in the 3D image

      double res = 0;
      double weightedSum=0;

      int x,y,z;
      itkFloatImage::SizeType    size    = referenceImage->GetLargestPossibleRegion().GetSize();

      #pragma omp parallel for private(x,y,z) reduction(+: res,weightedSum)
      for(z=0; z < (int)size[2]; z++)
      for(y=0; y < (int)size[1]; y++)
      for(x=0; x < (int)size[0]; x++)
      {
        itkFloatImage::IndexType refIndex; //index of the current voxel in the reference image
        refIndex[0] = x;
        refIndex[1] = y;
        refIndex[2] = z;
        if(referenceMask->GetPixel( refIndex ) > 0)
        {
            itkFloatImage::PointType refPoint; //physical point location of reference image
            referenceImage->TransformIndexToPhysicalPoint(refIndex,refPoint);

            itkFloatImage::PointType transformedPoint; //Physical point location after applying transform
            transformedPoint = transform->TransformPoint(refPoint);

            itkContinuousIndex       inputContIndex;   //continuous index in the 3D image
            movingImage->TransformPhysicalPointToContinuousIndex(transformedPoint,inputContIndex);

            if(bsInterpolatorMovingImage->IsInsideBuffer(inputContIndex))
            {
                double interpolatedValueImage = bsInterpolatorMovingImage->EvaluateAtContinuousIndex(inputContIndex);
                double weight = referenceMask->GetPixel( refIndex );
                double mse = weight * pow(interpolatedValueImage - referenceImage->GetPixel( refIndex ), 2.0);

                weightedSum += weight;
                res += mse;
            }
        }
      }
      return res / weightedSum;
    }
  };

} // namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkPandoraBoxCostFunction.txx"
#endif

#endif // BTK_PANDORA_BOX_COST_FUNCTION_H
