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

#ifndef BTKMOTIONCORRECTIONBYINTERSECTION_H
#define BTKMOTIONCORRECTIONBYINTERSECTION_H

/* ITK */
#include "itkImage.h"
#include "itkEuler3DTransform.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkContinuousIndex.h"

/* BTK */
#include "btkEulerSliceBySliceTransform.h"
#include "btkSlicesIntersectionITKCostFunction.hxx"
#include "btkMathFunctions.h"
#include "btkOptimizer.h"
#include "btkRigidRegistration.h"
#include "btkSimulatedAnnealingOptimizer.h"
#include "btkSmartStepGradientDescentOptimizer.h"


/* OTHERS */
#include "cfloat"
#include "cmath"
#include "algorithm"
#include "ctime"


namespace btk
{
/**
 * \brief Perform motion correction using intersection of slices
 *
 * This class perform a registration based on intersection of othogonals slices.
 * This class is templated over the input images.
 *
 * \ingroup Reconstruction
 */

template< typename TImage>
class MotionCorrectionByIntersection
{
public:
    /** Typedefs  */
    static const double MAX_DOUBLE = DBL_MAX;
    typedef TImage ImageType;
    typedef typename ImageType::PixelType Pixel;
    typedef typename ImageType::RegionType ImageRegion;
    typedef itk::Image<unsigned char, 3> MaskType;
    typedef btk::EulerSliceBySliceTransform<double, 3, Pixel> TransformType;
    //typedef itk::CenteredEuler3DTransform<double> Rigid3DTransform;
    typedef itk::Euler3DTransform<double> Rigid3DTransform;
    typedef typename TransformType::ParametersType ParametersType;


    /** Set/Get Methods for input Images */
    btkSetMacro(Images, std::vector< typename ImageType::Pointer >);
    btkGetMacro(Images, std::vector< typename ImageType::Pointer >);


    /** Set/Get methods for input masks */
    btkSetMacro(Masks, std::vector< MaskType::Pointer >);
    btkGetMacro(Masks, std::vector< MaskType::Pointer >);

    /** Get Method for transforms */
    btkGetMacro(Transforms, std::vector< typename TransformType::Pointer >);

    /** Get Method for the inverse transforms */
    btkGetMacro(InverseTransforms, std::vector< typename TransformType::Pointer >);

    /** Set/Get method for VerboseMode, default off (false) */
    btkSetMacro(VerboseMode, bool);
    btkGetMacro(VerboseMode, bool);

    /** Set/Get methods for number Max of loop, default 1 */
    btkSetMacro(MaxLoop,int);
    btkGetMacro(MaxLoop,int);

    /** Set Use of slice exclusion (default OFF) */
    btkSetMacro(UseSliceExclusion, bool);


    /** Get The outliers 2D vector */
    std::vector< std::vector< bool > > GetOutliers()
    {
        return m_Outliers;
    }

    /** Initialize method */
    virtual void Initialize();

    /** Update method, here start the optimization */
    virtual void Update();

    /** Method who correct slices with an to great error */
    void SlicesExclusion();

    /** Constructor */
    MotionCorrectionByIntersection();

    /** Destructor */
    virtual ~MotionCorrectionByIntersection(){}

	/* Francois debugging */
	double ComputeCostFunctionValueForOneSlice(unsigned int i, unsigned int smov);
	double ComputeOverallCostFunctionValue();

protected:
    /** UpdateInfos method, called after each iteration for updating transformations parameters informations */
    virtual void UpdateInfos();



private:
    /** return motion parameters for testing the algorithm (Only use for developement and testing) */
    typename TransformType::ParametersType SimulateMotion(double _Rmin, double _Rmax, double _Tmin, double _Tmax);
    /** Apply a Gaussian filter for multi-resolution */
    void BlurImages(double level);

    /** Parameters of Rigid Transformation to compute */
    vnl_vector<double> m_X;

    std::vector< typename ImageType::Pointer > m_Images; /** Input images */
    std::vector< typename ImageType::Pointer > m_BlurredImages; /** Input images */
    std::vector< MaskType::Pointer > m_Masks; /** Input Masks */
    std::vector<typename TransformType::Pointer> m_Transforms; /** Transforms to be computed */
    std::vector<typename TransformType::Pointer> m_InverseTransforms; /** Inverse Transform */
    bool m_VerboseMode; /** Verbose Mode: true : on, false : off */
    int m_NumberOfImages; /** Number of input images */
    int m_ReferenceImage; /** Reference image num */

    itk::Array<double> m_ScaleX; /** Scale parameters for optimizer */

    int m_ReferenceStack; /** Reference stack (image) */
    int m_ReferenceSlice; /** Reference slice */

    int m_MaxLoop; /** Number max of iteration*/

    std::vector<std::vector<double> > m_BestError; /** Array of best error (for detecting outliers) */

    bool m_VerboseDbg; /** Debug mode (display of further informations) */
    double m_CurrentError; /** Error at the current optimizer position */

    bool m_UseSliceExclusion; /** Perform or not Slice Exclusion (detection of outliers) */

    std::vector< std::vector<bool> > m_Outliers; /** Vector 2D of bool, if the value is true it is a outlier slice */

    unsigned int m_NumberOfParameters;

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkMotionCorrectionByIntersection.txx"
#endif

#endif // BTKMOTIONCORRECTIONBYINTERSECTION_H
