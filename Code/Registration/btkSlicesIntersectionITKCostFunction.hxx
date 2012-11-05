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

#ifndef BTKSLICESINTERSECTIONITKCOSTFUNCTION_HXX
#define BTKSLICESINTERSECTIONITKCOSTFUNCTION_HXX

#include "itkSingleValuedCostFunction.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkEuler3DTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"

#include "btkMacro.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkSlicesIntersectionVNLCostFunction.hxx"

#include  "cfloat"


namespace btk
{
template<class TImage>
class SlicesIntersectionITKCostFunction: public itk::SingleValuedCostFunction
{
    public:
        /** Standard class typedefs. */
        typedef SlicesIntersectionITKCostFunction   Self;
        typedef itk::SingleValuedCostFunction               Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef TImage ImageType;
        typedef typename ImageType::PixelType VoxelType;
        typedef itk::Image<unsigned char, 3> MaskType;
        typedef typename ImageType::RegionType ImageRegion;
        typedef itk::LinearInterpolateImageFunction< ImageType, double> Interpolator;

        typedef itk::ImageRegionIteratorWithIndex< ImageType >  Iterator;
        typedef itk::ImageRegionConstIteratorWithIndex< ImageType >  ConstIterator;
        typedef itk::ImageRegionIteratorWithIndex< MaskType >  MaskIterator;
        //typedef itk::MatrixOffsetTransformBase<double, 3> TransformType;
        typedef itk::Euler3DTransform<double> TransformType;
        typedef btk::EulerSliceBySliceTransform<double,3,VoxelType> SliceBySliceTransformType;

        /** Run-time type information (and related methods). */
        itkTypeMacro(SlicesIntersectionITKCostFunction, CostFunction);

        /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
        typedef double MeasureType;

        /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
        typedef Superclass::ParametersType      ParametersType;
        typedef Superclass::ParametersValueType ParametersValueType;

        /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
        typedef itk::Array< ParametersValueType > DerivativeType;

        /** This method returns the value of the cost function corresponding
      * to the specified parameters.    */
        virtual MeasureType GetValue(const ParametersType & parameters) const;

        /** This method returns the derivative of the cost function corresponding
      * to the specified parameters.   */
        virtual void GetDerivative(const ParametersType & parameters,
                                   DerivativeType & derivative) const;

        /** Return the number of parameters required to compute
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
        virtual unsigned int GetNumberOfParameters(void) const ;



        /** Get/Set Methods for VerboseMode, default false */
        btkGetMacro(VerboseMode, bool);
        btkSetMacro(VerboseMode, bool);

        /** Set Method for input images */
        btkSetMacro(Images, std::vector<typename ImageType::Pointer>);

        /** Set method for input masks */
        btkSetMacro(Masks, std::vector<MaskType::Pointer> );

        /** Set method for input transforms */
        btkSetMacro(Transforms, std::vector<typename SliceBySliceTransformType::Pointer>);

        /** Set method for input inverse transforms */
        btkSetMacro(InverseTransforms, std::vector<typename SliceBySliceTransformType::Pointer>);

        /** Set method for moving slice num (third slice index) */
        btkSetMacro(MovingSliceNum, int);
        /** Set method for moving image num in the images vector*/
        btkSetMacro(MovingImageNum, int);

        /** Set/Get method for the Number of parameters (default 6, euler transform 3D)
            Changing this number will introduce some error */
        btkSetMacro(NumberOfParameters,unsigned int);


        /** Intialization */
        void Initialize();

        /** New method for creating an object using a factory. */
        itkNewMacro(Self);



    protected:

        SlicesIntersectionITKCostFunction();
        virtual ~SlicesIntersectionITKCostFunction() {}

    private :

        int m_NumberOfIntersectedVoxels;
        VoxelType m_SumOfIntersectedVoxels;
        bool m_VerboseMode;

        typename ImageType::Pointer m_ReferenceImage;
        typename ImageType::Pointer m_MovingImage;

        typename Interpolator::Pointer m_ReferenceInterpolator;
        typename Interpolator::Pointer m_MovingInterpolator;

        MaskType::Pointer m_ReferenceMask;
        MaskType::Pointer m_MovingMask;

        TransformType::Pointer m_Transform;
        TransformType::ParametersType m_Parameters;

        //typename ImageType::PointType
        std::vector<double> m_CenterOfMovingSlice;
        int m_MovingSliceNum;
        int m_MovingImageNum;
        int m_NumberOfImages;

        std::vector<typename ImageType::Pointer> m_Images;
        std::vector<typename SliceBySliceTransformType::Pointer> m_Transforms; // Inverse SliceBySliceTransforms
        std::vector<typename SliceBySliceTransformType::Pointer> m_InverseTransforms; // Inverse SliceBySliceTransforms
        std::vector<MaskType::Pointer> m_Masks;
        std::vector<typename Interpolator::Pointer> m_Interpolators;
        float m_lambda;

        SlicesIntersectionVNLCostFunction<ImageType> m_VNLCostFunction;
        unsigned int m_NumberOfParameters;



};
}
#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSlicesIntersectionITKCostFunction.txx"
#endif




#endif // BTKSLICESINTERSECTIONITKCOSTFUNCTION_HXX
