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

#ifndef BTK_BTKSIMULATELOWRESOLUTIONIMAGESFILTER_HXX
#define BTK_BTKSIMULATELOWRESOLUTIONIMAGESFILTER_HXX


#include "itkFixedArray.h"
#include "itkTransform.h"
#include "itkImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageToImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkSize.h"
#include "itkImageDuplicator.h"
#include "itkContinuousIndex.h"

#include "vnl/vnl_sparse_matrix.h"

#include "btkLinearInterpolateImageFunctionWithWeights.h"
#include "btkOrientedSpatialFunction.h"
#include "btkMacro.h"
#include "btkEulerSliceBySliceTransform.h"


namespace btk
{

/**
 * OLD Code for compute H and compute simulated image with existing High resolution images
 * Use for Simulating low resolution image
 */

template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType=double>
class SimulateLowResolutionImagesFilter : public itk::ImageToImageFilter<TInputImage, TOutputImage>
{
    public:
        /** Standard class typedefs. */
        typedef SimulateLowResolutionImagesFilter                  Self;
        typedef itk::ImageToImageFilter<TInputImage,TOutputImage>  Superclass;
        typedef SmartPointer<Self>                                 Pointer;
        typedef SmartPointer<const Self>                           ConstPointer;

        typedef TInputImage                             InputImageType;
        typedef TOutputImage                            OutputImageType;
        typedef typename InputImageType::Pointer        InputImagePointer;
        typedef typename InputImageType::ConstPointer   InputImageConstPointer;
        typedef typename OutputImageType::Pointer       OutputImagePointer;

        typedef typename InputImageType::RegionType     InputImageRegionType;
        typedef std::vector<InputImageRegionType>       InputImageRegionVectorType;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(SimulateLowResolutionImagesFilter, ImageToImageFilter);

        /** Number of dimensions. */
        itkStaticConstMacro(ImageDimension, unsigned int,
                            TOutputImage::ImageDimension);
        itkStaticConstMacro(InputImageDimension, unsigned int,
                            TInputImage::ImageDimension);


        /** Transform typedef. */
        typedef btk::EulerSliceBySliceTransform<double,3> TransformType;

        /** Image size typedef. */
        typedef Size<itkGetStaticConstMacro(ImageDimension)> SizeType;

        /** Image index typedef. */
        typedef typename TOutputImage::IndexType IndexType;

        /** Image continuous index typedef. */
        typedef ContinuousIndex<double, ImageDimension> ContinuousIndexType;

        /** Image point typedef. */
        typedef typename TOutputImage::PointType    PointType;

        /** Image pixel value typedef. */
        typedef typename TOutputImage::PixelType   PixelType;
        typedef typename TInputImage::PixelType    InputPixelType;

        /** Typedef to describe the output image region type. */
        typedef typename TOutputImage::RegionType OutputImageRegionType;

        /** Image spacing, origin and direction typedef */
        typedef typename TOutputImage::SpacingType   SpacingType;
        typedef typename TOutputImage::PointType     OriginPointType;
        typedef typename TOutputImage::DirectionType DirectionType;

        /** base type for images of the current ImageDimension */
        typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

        /**Const iterator typedef. */
        typedef ImageRegionConstIteratorWithIndex< InputImageType >  ConstIteratorType;

        /**Const iterator typedef. */
        typedef ImageRegionConstIteratorWithIndex< OutputImageType >  OutputIteratorType;

        /**Neighborhood iterator typedef. */
        typedef NeighborhoodIterator< InputImageType >   NeighborhoodIteratorType;

        /**Blurring function typedef. */
        typedef OrientedSpatialFunction<double, 3, PointType> FunctionType;

        /**Interpolator typedef. */
        typedef LinearInterpolateImageFunctionWithWeights<TOutputImage, double> InterpolatorType;

        /** Duplicator typedef. */
        typedef ImageDuplicator< InputImageType > DuplicatorType;

        /**Face calculator typedef. */
        typedef NeighborhoodAlgorithm
        ::ImageBoundaryFacesCalculator< OutputImageType > FaceCalculatorType;
        typedef typename FaceCalculatorType::FaceListType FaceListType;

        typedef vnl_matrix<double> VnlMatrixType;
        typedef itk::Matrix<double,3,3> MatrixType;
        typedef vnl_vector<double> VnlVectorType;

        typedef vnl_sparse_matrix<double> VnlSparseMatrixType;
    public:
        /** Constructor */
        SimulateLowResolutionImagesFilter();
        /** Update method (compute images) */
        void Update();
        /** Add Input, you should call it to set inputs */
        void AddInput(InputImageType* _arg)
        {
          m_Images.push_back(_arg);
          m_InputImageRegion.push_back(_arg->GetLargestPossibleRegion());
          this -> SetInput(_arg);

          // Add transforms for this image
          m_Transforms.resize( m_Transforms.size() + 1 );
          m_Transforms[m_Transforms.size()-1] = TransformType::New();
          m_Transforms[m_Transforms.size()-1]->SetImage(_arg);

          typename DuplicatorType::Pointer duplicator = DuplicatorType::New();
          duplicator -> SetInputImage (_arg);
          duplicator -> Update();
          m_SimulatedImages.push_back( duplicator -> GetOutput() );
          m_SimulatedImages[m_SimulatedImages.size()-1]->FillBuffer(0.0);

        }

        /** Set transformation of the ith image */
        void SetTransform( int i, TransformType* transform )
        {
          m_Transforms[i] = transform;
          m_Transforms[i]->SetImage(m_Images[i]);
        }

        /** Get the ith simulated image */
        InputImageType* GetSimulatedImage( int i )
        {
          return m_SimulatedImages[i];
        }
        /** Set the reference image (the reconstructed one) */
        void SetReferenceImage(InputImagePointer _ref)
        {
            this->m_ReferenceImage = _ref;
        }
        btkGetMacro(ReferenceImage,InputImagePointer);

    protected:
        /** Compute the simulated images */
        void UpdateSimulatedImages();
        /** Create H */
        void ComputeH();

    private :
        SimulateLowResolutionImagesFilter( const Self& ); //purposely not implemented
        void operator=( const Self& ); //purposely not implemented


        SizeType                    m_Size;              // Size of the output image
        std::vector< typename TransformType::Pointer> m_Transforms;
        InputImageRegionVectorType  m_InputImageRegion;

        OutputImageRegionType       m_OutputImageRegion;

        std::vector<InputImagePointer>     m_Images;
        std::vector<InputImagePointer>     m_SimulatedImages;
        InputImagePointer                  m_ReferenceImage;
        bool    m_SimulatedImagesUpdated;

        VnlSparseMatrixType m_H;
        VnlSparseMatrixType m_Ht;
        VnlSparseMatrixType m_Hbp;
        VnlVectorType       m_y;
        VnlVectorType       m_ysim;
        VnlVectorType       m_x;

        float m_Lambda;

        // Precomputed values for optimization

        VnlSparseMatrixType HtH;
        VnlVectorType HtY;
        double YtY;
        unsigned int m_PSF;

};

} // namespace btk



#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSimulateLowResolutionImagesFilter.txx"
#endif

#endif // BTK_BTKSIMULATELOWRESOLUTIONIMAGESFILTER_HXX
