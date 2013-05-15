/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 19/11/2012
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


#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkResampleImageFilter.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkTransform.h"
#include "itkEuler3DTransform.h"
#include "itkAmoebaOptimizer.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegistrationMethod.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "btkImageHelper.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkCenteredEulerSliceBySliceTransform.h"
#include "btkCommandIterationUpdate.h"
#include "btkSlicesIntersectionITKCostFunction.hxx"

/* OTHERS */
#include "iostream"
#include "sstream"
#include "fstream"
#include <tclap/CmdLine.h>

typedef float PixelType;
const unsigned int Dimension = 3;
typedef itk::Image<PixelType, Dimension> itkImage;
typedef itk::Image<unsigned char, Dimension> MaskType;
typedef itk::ImageRegionIteratorWithIndex< itkImage >  Iterator;
typedef itk::Euler3DTransform<double> TransformType;
typedef btk::EulerSliceBySliceTransform< double, 3, PixelType> Transform;
//typedef btk::CenteredEulerSliceBySliceTransform<double, 3, PixelType> Transform;
typedef itk::ResampleImageFilter<itkImage, itkImage> Resampler;
typedef itk::MatrixOffsetTransformBase<double> MatrixTransformType;
typedef itk::Transform<double> itkTransform;
typedef btk::SlicesIntersectionITKCostFunction<itkImage> CostFunction;
typedef Transform::ParametersType Parameters;


class CostFunctionTester
{
    public:
        CostFunctionTester(){}
        ~CostFunctionTester(){}

        double GetValue(Parameters _x, int _image, int _slice)
        {
            std::vector< unsigned int > SlicesGroup;
             itkImage::SizeType sizeMov = m_Images[_image]->GetLargestPossibleRegion().GetSize();
             unsigned int nbGroup = sizeMov[2];
             SlicesGroup.resize(sizeMov[2]);

             for(unsigned int smov = 0; smov < sizeMov[2]; smov++ )
            {
                SlicesGroup[smov] = smov%nbGroup;
            }

            CostFunction::Pointer f = CostFunction::New();
            f->SetVerboseMode(false);
            f->SetNumberOfParameters(6);
            f->SetImages(this->m_Images);
            f->SetMasks(this->m_Masks);
            f->SetTransforms(this->m_Transforms);
            f->SetInverseTransforms(this->m_InverseTransform);
            f->SetMovingImageNum(_image);// firstImage
            //f->SetMovingSliceNum(_slice);
            f->SetGroupNum(_slice);
            f->SetSlicesGroup(SlicesGroup);
            f->Initialize();

            return f->GetValue(_x);
        }

        Parameters GetDerivative(Parameters _x , int _image, int _slice)
        {
            std::vector< unsigned int > SlicesGroup;
             itkImage::SizeType sizeMov = m_Images[_image]->GetLargestPossibleRegion().GetSize();
             unsigned int nbGroup = sizeMov[2];
             SlicesGroup.resize(sizeMov[2]);
            for(unsigned int smov = 0; smov < sizeMov[2]; smov++ )
            {
                SlicesGroup[smov] = smov%nbGroup;
            }

            CostFunction::Pointer f = CostFunction::New();
            f->SetVerboseMode(false);
            f->SetNumberOfParameters(6);
            f->SetImages(this->m_Images);
            f->SetMasks(this->m_Masks);
            f->SetTransforms(this->m_Transforms);
            f->SetInverseTransforms(this->m_InverseTransform);
            f->SetMovingImageNum(_image);// firstImage
            //f->SetMovingSliceNum(_slice);
            f->SetGroupNum(_slice);
            f->SetSlicesGroup(SlicesGroup);

            f->Initialize();
            Parameters _gx(6);
            f->GetDerivative(_x, _gx);
            return _gx;

        }

        void SetImages(std::vector<itkImage::Pointer> _im)
        {
            m_Images = _im;
        }

        void SetMasks(std::vector<MaskType::Pointer> _masks)
        {
            m_Masks = _masks;
        }

        void SetTransforms(std::vector<Transform::Pointer> _t)
        {
            m_Transforms = _t;
        }

        void SetInverseTransforms(std::vector<Transform::Pointer> _it)
        {
            m_InverseTransform =  _it;
        }

    private:
        std::vector<itkImage::Pointer> m_Images;
        std::vector<Transform::Pointer> m_Transforms;
        std::vector<MaskType::Pointer> m_Masks;
        std::vector<Transform::Pointer> m_InverseTransform;


};


int main(int argc, char * argv[])
{

    // TCLAP :
    TCLAP::CmdLine cmd("Test cost function performed by slice intersection, save text file with cost and gradient", ' ', "Unversioned");
    TCLAP::MultiArg<std::string> inputArg("i","input","Low-resolution image file",true,"string",cmd);
    TCLAP::MultiArg<std::string> maskArg("m","mask","Low-resolution mask file",true,"string",cmd);
    TCLAP::ValueArg<unsigned int> sliceArg("s","slice","slice number",true,1,"uint",cmd);
    TCLAP::ValueArg<unsigned int> imageArg("n","image","image number",true,0,"uint",cmd);

    std::vector<std::string> input,masks;

    cmd.parse(argc,argv);
    input = inputArg.getValue();
    masks = maskArg.getValue();
    unsigned int slice = sliceArg.getValue();
    unsigned int imageNum = imageArg.getValue();
    std::vector< itkImage::Pointer > inputImages;
    std::vector< MaskType::Pointer >maskImages;

    inputImages = btk::ImageHelper<itkImage>::ReadImage(input);

    maskImages = btk::ImageHelper<MaskType>::ReadImage(masks);

    std::vector<Transform::Pointer> transforms(input.size());
    std::vector<Transform::Pointer> inverseTransforms(input.size());



    for(unsigned int i = 0; i< transforms.size(); i++)
    {
        inverseTransforms[i] = Transform::New();
        inverseTransforms[i]->SetImage(inputImages[i]);
        inverseTransforms[i]->Initialize();
        inverseTransforms[i]->SetIdentity();

        transforms[i] = Transform::New();
        transforms[i]->SetImage(inputImages[i]);
        transforms[i]->Initialize();
        transforms[i]->SetIdentity();


        typedef itk::DiscreteGaussianImageFilter<itkImage, itkImage> gaussianFilter;
        gaussianFilter::ArrayType sigma;
        sigma[0] = 0.;
        sigma[1] = 0.;
        sigma[2] = 0.0001;



    }


    int NumberOfSlices = inputImages[imageNum]->GetLargestPossibleRegion().GetSize()[2];

    CostFunctionTester* cft = new CostFunctionTester();
    cft->SetImages(inputImages);
    cft->SetMasks(maskImages);
    cft->SetTransforms(transforms);
    cft->SetInverseTransforms(inverseTransforms);

    std::vector<double> Min(6), Max(6), Step(6);

    for(unsigned int i = 0; i< 6; i++)
    {

        if(i <3)
        {

            Min[i] = -20.0;
            Max[i] = 20.0;
            Step[i] = 1.0;

        }
        else
        {
            Min[i] = -20.0;
            Max[i] = 20.0;
            Step[i] = 2.0;
        }
    }

    Parameters X(6),gX(6);
    Parameters MinX(6);
    MinX.Fill(0.0);
    X.Fill(0.0);
    double CostValue = 0.0;
    double MinValue = 0.0;

    std::ofstream Rx, Ry, Rz, Tx, Ty, Tz;
    std::ofstream gRx, gRy, gRz, gTx, gTy, gTz;
    std::ofstream RxRy, RxRz, RxTx, RxTy, RxTz;
    std::ofstream RyRy, RyTx, RyTy, RyTz;
    std::ofstream RzRy, RyRz, RzTx, RzTy, RzTz;
    std::ofstream TxTy, TxTz, TzTy;

//    X[1] = 0;
//    CostValue = cft->GetValue(X, 0, slice);

//    return 0;



/**************************************************/
//    TzTy.open("TzTy.txt");
//    for(double i = Min[5]; i< Max[5]; i+= Step[5])
//    {
//        for(double j = Min[4]; j< Max[4]; j+= Step[4])
//        {
//            X.Fill(0.0);
//            X[3] = transforms[0]->GetSliceParameters(slice)[3];
//            X[4] = transforms[0]->GetSliceParameters(slice)[4];
//            X[5] = transforms[0]->GetSliceParameters(slice)[5];
//            X[8] = i;
//            X[7] = j;
//            CostValue = cft->GetValue(X, 0, slice);
//            if(TzTy.is_open())
//            {
//                TzTy<<i<<" "<<j<<" "<<CostValue<<std::endl;
//            }
//            else
//            {
//                std::cout<<"Unable to open TzTy.txt"<<std::endl;
//            }

//        }
//    }

////        for(unsigned int slice = 0; slice<NumberOfSlices; slice++)
////        {


//   return 0;
    Rx.open("Rx.txt");
    gRx.open("gRx.txt");
    for(double i = Min[0]; i <Max[0] +Step[0] ; i+=Step[0])
    {
        X.Fill(0.0);
//        X[3] = transforms[0]->GetSliceParameters(slice)[3];
//        X[4] = transforms[0]->GetSliceParameters(slice)[4];
//        X[5] = transforms[0]->GetSliceParameters(slice)[5];
        X[0]= i;
        CostValue = cft->GetValue(X, imageNum, slice);

        gX = cft->GetDerivative(X,imageNum,slice);


        //std::cout<<"Parameters : "<<X<<" Value : "<<CostValue<<std::endl;

        if(Rx.is_open())
        {
            Rx<<i<<" "<<CostValue<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open Rx.txt"<<std::endl;
        }

        if(gRx.is_open())
        {
            gRx<<i<<" "<<gX[0]<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open gRx.txt"<<std::endl;
        }

    }
    Rx.close();
    gRx.close();

    //rotY

    Ry.open("Ry.txt");
    gRy.open("gRy.txt");
    for(double i = Min[1]; i <Max[1] + Step[1] ; i+=Step[1])
    {
        X.Fill(0.0);
//        X[3] = transforms[0]->GetSliceParameters(slice)[3];
//        X[4] = transforms[0]->GetSliceParameters(slice)[4];
//        X[5] = transforms[0]->GetSliceParameters(slice)[5];
        X[1]= i;
        CostValue = cft->GetValue(X, 0, slice);
        gX = cft->GetDerivative(X,0,slice);
        //std::cout<<"Parameters : "<<X<<" Value : "<<CostValue<<std::endl;

        if(Ry.is_open())
        {
            Ry<<i<<" "<<CostValue<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open Ry.txt"<<std::endl;
        }
        if(gRy.is_open())
        {
            gRy<<i<<" "<<gX[1]<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open gRy.txt"<<std::endl;
        }
    }
    Ry.close();
    gRy.close();

    //rotZ
    Rz.open("Rz.txt");
    gRz.open("gRz.txt");
    for(double i = Min[2]; i <Max[2] + Step[2]; i+=Step[2])
    {
        X.Fill(0.0);
//        X[3] = transforms[0]->GetSliceParameters(slice)[3];
//        X[4] = transforms[0]->GetSliceParameters(slice)[4];
//        X[5] = transforms[0]->GetSliceParameters(slice)[5];
        X[2]= i;
        CostValue = cft->GetValue(X, imageNum, slice);
        gX = cft->GetDerivative(X,imageNum,slice);
        //std::cout<<"Parameters : "<<X<<" Value : "<<CostValue<<std::endl;

        if(Rz.is_open())
        {
            Rz<<i<<" "<<CostValue<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open Rz.txt"<<std::endl;
        }
        if(gRz.is_open())
        {
            gRz<<i<<" "<<gX[2]<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open gRz.txt"<<std::endl;
        }
    }
    Rz.close();
    gRz.close();

    //traX
    Tx.open("Tx.txt");
    gTx.open("gTx.txt");
    for(double i = Min[3]; i <Max[3] + Step[3]; i+=Step[3])
    {
        X.Fill(0.0);
//        X[3] = transforms[0]->GetSliceParameters(slice)[3];
//        X[4] = transforms[0]->GetSliceParameters(slice)[4];
//        X[5] = transforms[0]->GetSliceParameters(slice)[5];
        X[3]= i;
        CostValue = cft->GetValue(X, imageNum, slice);
        gX = cft->GetDerivative(X,imageNum,slice);
        //std::cout<<"Parameters : "<<X<<" Value : "<<CostValue<<std::endl;
        if(Tx.is_open())
        {
            Tx<<i<<" "<<CostValue<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open Tx.txt"<<std::endl;
        }
        if(gTx.is_open())
        {
            gTx<<i<<" "<<gX[3]<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open gTx.txt"<<std::endl;
        }
    }
    Tx.close();
    gTx.close();

    //traY
    Ty.open("Ty.txt");
    gTy.open("gTy.txt");
    for(double i = Min[4]; i <Max[4]+ Step[4]; i+=Step[4])
    {
        X.Fill(0.0);
//        X[3] = transforms[0]->GetSliceParameters(slice)[3];
//        X[4] = transforms[0]->GetSliceParameters(slice)[4];
//        X[5] = transforms[0]->GetSliceParameters(slice)[5];
        X[4]= i;
        CostValue = cft->GetValue(X, imageNum, slice);
        gX = cft->GetDerivative(X,imageNum,slice);
        //std::cout<<"Parameters : "<<X<<" Value : "<<CostValue<<std::endl;
        if(Ty.is_open())
        {
            Ty<<i<<" "<<CostValue<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open Ty.txt"<<std::endl;
        }
        if(gTy.is_open())
        {
            gTy<<i<<" "<<gX[4]<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open gTy.txt"<<std::endl;
        }
    }
    Ty.close();
    gTy.close();

    //traZ
    Tz.open("Tz.txt");
    gTz.open("gTz.txt");
    for(double i = Min[5]; i <Max[5]+ Step[5]; i+=Step[5])
    {
        X.Fill(0.0);
//        X[3] = transforms[0]->GetSliceParameters(slice)[3];
//        X[4] = transforms[0]->GetSliceParameters(slice)[4];
//        X[5] = transforms[0]->GetSliceParameters(slice)[5];
        X[5]= i;
        CostValue = cft->GetValue(X, imageNum, slice);
        gX = cft->GetDerivative(X,imageNum,slice);
        //std::cout<<"Parameters : "<<X<<" Value : "<<CostValue<<std::endl;

        if(Tz.is_open())
        {
            Tz<<i<<" "<<CostValue<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open Tz.txt"<<std::endl;
        }
        if(gTz.is_open())
        {
            gTz<<i<<" "<<gX[5]<<std::endl;
        }
        else
        {
            std::cout<<"Unable to open gTz.txt"<<std::endl;
        }
    }
    Tz.close();
    gTz.close();

    //}

    delete cft;
}
