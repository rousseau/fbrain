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
#include "itkPowellOptimizer.h"
#include "itkAmoebaOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkExhaustiveOptimizer.h"
#include "itkDiscreteGaussianImageFilter.h"

/* BTK */
#include "btkEulerSliceBySliceTransform.h"
#include "btkSlicesIntersectionITKCostFunction.hxx"
#include "btkMathFunctions.h"
#include "btkCenteredEulerSliceBySliceTransform.h"
#include "btkOptimizer.h"



/* VNL */

#include "vnl/algo/vnl_powell.h"

#include  "cfloat"
#include "cmath"
#include "algorithm"
#include "ctime"


namespace btk
{
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
    //typedef btk::EulerSliceBySliceTransform<double, 3, Pixel> Transform;
    typedef btk::CenteredEulerSliceBySliceTransform<double, 3, Pixel> Transform;


    /** Set/Get Methods for input Images */
    btkSetMacro(Images, std::vector< typename ImageType::Pointer >);
    btkGetMacro(Images, std::vector< typename ImageType::Pointer >);


    /** Set/Get methods for input masks */
    btkSetMacro(Masks, std::vector< MaskType::Pointer >);
    btkGetMacro(Masks, std::vector< MaskType::Pointer >);

    /** Get Method for transforms */
    btkGetMacro(Transforms, std::vector< typename Transform::Pointer >);

    /** Set/Get method for VerboseMode, default off (false) */
    btkSetMacro(VerboseMode, bool);
    btkGetMacro(VerboseMode, bool);

    /** Set/Get methods for number Max of loop, default 1 */
    btkSetMacro(MaxLoop,int);
    btkGetMacro(MaxLoop,int);

    /** Set Use of slice exclusion (default OFF) */
    btkSetMacro(UseSliceExclusion, bool);

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


protected:
    /** UpdateInfos method, called after each iteration for updating transformations parameters informations */
    virtual void UpdateInfos();



private:
    /** return motion parameters for testing the algorithm (Only use for developement) */
    typename Transform::ParametersType SimulateMotion(double _Rmin, double _Rmax, double _Tmin, double _Tmax);
    /** Apply a Gaussian filter for multi-resolution */
    void BlurImages(double level);
    /** Parameters of Rigid Transformation to compute */
    vnl_vector<double> m_X;

    std::vector< typename ImageType::Pointer > m_Images; /** Input images */
    std::vector< typename ImageType::Pointer > m_BlurredImages; /** Input images */
    std::vector< MaskType::Pointer > m_Masks; /** Input Masks */
    std::vector<typename Transform::Pointer> m_Transforms; /** Transforms to be computed */
    std::vector<typename Transform::Pointer> m_InverseTransforms; /** Inverse Transform */
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



};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkMotionCorrectionByIntersection.txx"
#endif

#endif // BTKMOTIONCORRECTIONBYINTERSECTION_H
