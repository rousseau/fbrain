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
#include "btkImageHelper.h"
#include "btkMathFunctions.h"


#include "cfloat"
#include "sstream"
#include "numeric"


//#include "itkJoinImageFilter.h"
//#include "itkImageToHistogramFilter.h"
//#include "itkMinimumMaximumImageCalculator.h"

namespace btk
{
/**
 * Class for compute a cost function on intersection of slices.
 * This class is specificly designed for motion correction on fetal MRI.
 * This cost function can be used by a vnl_optimizer, or an ITK::Optimizer
 * (using the specific adaptor class btk::SliceIntersectionITKCostFunction).
 *
 * @author Marc Schweitzer
 * \ingroup Registration
 */
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

    typedef itk::Euler3DTransform<double> TransformType;
   typedef btk::EulerSliceBySliceTransform<double,3,VoxelType> SliceBySliceTransformType;

    typedef itk::ContinuousIndex<double, 3 > ContinuousIndexType;

//    typedef itk::JoinImageFilter< ImageType, ImageType >  JoinFilterType;

//    typedef typename JoinFilterType::OutputImageType               VectorImageType;

//    typedef typename itk::Statistics::ImageToHistogramFilter< VectorImageType >  HistogramFilterType;

//    typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxFilter;

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

    /** Set/Get the center of the transform (default, middle of the volume ) */
    btkSetMacro(CenterOfTransform, typename ImageType::PointType);
    btkGetMacro(CenterOfTransform, typename ImageType::PointType);

    /** Constructor */
    SlicesIntersectionVNLCostFunction(unsigned int dim);
    SlicesIntersectionVNLCostFunction(){}

    /** Destructor */
    ~SlicesIntersectionVNLCostFunction();

    /** Cost Function */
    double f(const vnl_vector<double> &x) const;

    virtual vnl_vector<double>GetGradient(vnl_vector<double> const& x,double stepsize = 0.05) const;
     /** Intialization */
     void Initialize();

     /** Get If the last evaluation of cost function has intersection or not */
     btkGetMacro(Intersection, bool);



protected:

     /** Compute squared Difference between two voxels */
      inline double SquaredDifference(VoxelType ReferenceVoxel, VoxelType MovingVoxel) const
      {
          return (double)( (ReferenceVoxel - MovingVoxel ) * (ReferenceVoxel - MovingVoxel ) );
      }
      /** Compute Absolute Difference between two voxels */
      inline double AbsoluteDifference(VoxelType ReferenceVoxel, VoxelType MovingVoxel) const
      {
          return (double)( std::abs(ReferenceVoxel - MovingVoxel ));
      }

	  //Cette fonction est complètement inutile !!!!!!!!!!! -----------------------------------------------
      /** Compute the square root difference between two voxels */
      inline double RootSquaredDifference(VoxelType ReferenceVoxel, VoxelType MovingVoxel) const
      {
          return (double)( std::sqrt((ReferenceVoxel - MovingVoxel) * (ReferenceVoxel - MovingVoxel )));
      }
      /** Found starting and ending point of an intersection line (if existing) */
      bool FoundIntersectionPoints(unsigned int _fixedImage, unsigned int _fixedSlice, typename ImageType::PointType &Point1, typename ImageType::PointType &Point2) const;




private :

    bool m_VerboseMode; /** Boolean to activate, desactivate the verbose mode */

    TransformType::Pointer m_X; /** Transformation with parameters found by the optimizer */
    TransformType::Pointer m_InverseX; /** Inverse of the transformation */

    TransformType::ParametersType m_Parameters; /** Parameters found by optimizer */

    int m_MovingSliceNum; /** Moving slice num */
    int m_MovingImageNum;/** Moving image num */
    int m_NumberOfImages;/** Number of images */


    std::vector<typename ImageType::Pointer> m_Images; /** Vector of pointer of images */
    std::vector<typename SliceBySliceTransformType::Pointer> m_Transforms; /** Vector of pointer of transformations. Warning: if you want to use transformation of moving image use m_X instead! */
    std::vector<typename SliceBySliceTransformType::Pointer> m_InverseTransforms; /** Vector of inverse transforms */
    std::vector<MaskType::Pointer> m_Masks; /** Vector of Masks */
    std::vector<typename Interpolator::Pointer> m_Interpolators; /** Vector of interpolators */
    int m_NumberOfPointsOfLine; /** Number Of points of the intersection line (currently unused) */

    typename ImageType::PointType m_CenterOfTransform;
    mutable bool m_Intersection; /** Return true if the moving slice has intersections or not. This variable is mutable while we modify it into a const function (f()) */


};

}



#ifndef ITK_MANUAL_INSTANTIATION
#include "btkSlicesIntersectionVNLCostFunction.txx"
#endif


#endif // BTKSlicesIntersectionVNLCostFunction_HXX
