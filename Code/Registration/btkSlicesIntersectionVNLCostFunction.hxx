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

#ifndef BTKSlicesIntersectionVNLCostFunction_HXX
#define BTKSlicesIntersectionVNLCostFunction_HXX

#include "vnl/vnl_cost_function.h"
#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkEuler3DTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkContinuousIndex.h"
#include "itkVector.h"
#include "itkLineConstIterator.h"
#include "itkNeighborhoodIterator.h"
#include "itkStatisticsImageFilter.h"


#include "btkMacro.h"
#include "btkEulerSliceBySliceTransform.h"
#include "btkCenteredEulerSliceBySliceTransform.h"
#include "btkImageHelper.h"
#include "btkMathFunctions.h"


#include  "cfloat"
#include "sstream"

namespace btk
{
template <class TImage>
class SlicesIntersectionVNLCostFunction: public vnl_cost_function
{
public:

    /** Max of the cost function */
    static const double MAX_COSTFUNCTION_VALUE = DBL_MAX;

    /** typedefs  */
    typedef TImage ImageType;
    typedef typename ImageType::PixelType VoxelType;
    typedef itk::Image<unsigned char, 3> MaskType;
    typedef typename ImageType::RegionType ImageRegion;
    typedef itk::LinearInterpolateImageFunction< ImageType, double> Interpolator;
    typedef itk::NeighborhoodIterator<ImageType> Neighborhood;
    typedef itk::StatisticsImageFilter<ImageType> Statistics;


    typedef itk::ImageRegionIteratorWithIndex< ImageType >  Iterator;
    typedef itk::ImageRegionConstIteratorWithIndex< ImageType >  ConstIterator;
    typedef itk::ImageRegionIteratorWithIndex< MaskType >  MaskIterator;

    //typedef itk::Euler3DTransform<double> TransformType;
    typedef itk::CenteredEuler3DTransform<double> TransformType;
    //typedef btk::EulerSliceBySliceTransform<double,3,VoxelType> SliceBySliceTransformType;
    typedef btk::CenteredEulerSliceBySliceTransform<double, 3, VoxelType> SliceBySliceTransformType;

    typedef itk::ContinuousIndex<double, 3 > ContinuousIndexType;

    /** Get/Set Method for the verboseMode */
    btkGetMacro(VerboseMode, bool);
    btkSetMacro(VerboseMode, bool);

    /** Set method for input images */
    btkSetMacro(Images, std::vector<typename ImageType::Pointer>);

    /** Set method for input masks */
    btkSetMacro(Masks, std::vector<MaskType::Pointer> );

    /** Set method for Transforms */
    btkSetMacro(Transforms, std::vector<typename SliceBySliceTransformType::Pointer>);
    /** Set method for Inverse transforms */
    btkSetMacro(InverseTransforms, std::vector<typename SliceBySliceTransformType::Pointer>);

    /** Set the num of the slice (third component of size) */
    btkSetMacro(MovingSliceNum, int);

    /** Set the moving Image Num (images vector position) */
    btkSetMacro(MovingImageNum, int);

    /** Set/Get Number of points on the intersection line */
    btkSetMacro(NumberOfPointsOfLine, int);
    btkGetMacro(NumberOfPointsOfLine, int);

    /** Constructor */
    SlicesIntersectionVNLCostFunction(unsigned int dim);
    SlicesIntersectionVNLCostFunction(){}

    /** Destructor */
    ~SlicesIntersectionVNLCostFunction();

    /** Cost Function */
    double f(const vnl_vector<double> &x) const;

    virtual vnl_vector<double>GetGradient(vnl_vector<double> const& x,double stepsize = 10e-8) const;
     /** Intialization */
     void Initialize();


//     /** Set the reference slice (for the reference coordinate system) */
//     btkSetMacro(ReferenceCoordSlice, int);

//     /** Set the reference stack (for the reference coordinate system) */
//     btkSetMacro(ReferenceCoordImage, int);


protected:

     /** ComputeDifference between two voxels */
      inline double SquareDifference(double ReferenceVoxel, double MovingVoxel) const
      {
          return (double)( (ReferenceVoxel - MovingVoxel ) * (ReferenceVoxel - MovingVoxel ) );
      }



private :

    bool m_VerboseMode;

    typename ImageType::Pointer m_ReferenceImage;
    typename ImageType::Pointer m_MovingImage;

    typename Interpolator::Pointer m_ReferenceInterpolator;
    typename Interpolator::Pointer m_MovingInterpolator;

    MaskType::Pointer m_ReferenceMask;
    MaskType::Pointer m_MovingMask;

    TransformType::Pointer m_Transform;
    //typename SliceBySliceTransformType::Pointer m_Transform;

    TransformType::ParametersType m_Parameters;

    //typename ImageType::PointType
    std::vector<double> m_CenterOfMovingSlice;
    int m_MovingSliceNum;
    int m_MovingImageNum;
    int m_NumberOfImages;
//    int m_ReferenceCoordSlice;
//    int m_ReferenceCoordImage;

    std::vector<typename ImageType::Pointer> m_Images;
    std::vector<typename SliceBySliceTransformType::Pointer> m_Transforms;
    std::vector<typename SliceBySliceTransformType::Pointer> m_InverseTransforms;
    std::vector<MaskType::Pointer> m_Masks;
    std::vector<typename Interpolator::Pointer> m_Interpolators;
    int m_NumberOfPointsOfLine;


};

}



#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSlicesIntersectionVNLCostFunction.txx"
#endif


#endif // BTKSlicesIntersectionVNLCostFunction_HXX
