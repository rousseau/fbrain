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

/* BTK */
#include "btkEulerSliceBySliceTransform.h"
#include "btkSlicesIntersectionITKCostFunction.hxx"
#include "btkMathFunctions.h"
#include "btkCenteredEulerSliceBySliceTransform.h"



/* VNL */

#include "vnl/algo/vnl_powell.h"

#include  "cfloat"
#include "cmath"
#include "algorithm"


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

    /** Set/Get methods for number Max of loop, default 5 */
    btkSetMacro(MaxLoop,int);
    btkGetMacro(MaxLoop,int);

    /** Set Use of slice exclusion (default OFF) */
    btkSetMacro(UseSliceExclusion, bool);

    /** Initialize method */
    virtual void Initialize();

    /** Update method, here start the optimization */
    virtual void Update();

    /** Method who correct slices with an two great error */
    void SlicesExclusion();

    /** Constructor */
    MotionCorrectionByIntersection();

    /** Destructor */
    virtual ~MotionCorrectionByIntersection(){}


protected:
    /** UpdateInfos method, called after each iteration for updating transformations parameters informations */
    virtual void UpdateInfos();


private:
    /** Parameters of Rigid Transformation to compute */
    vnl_vector<double> m_X;

    std::vector< typename ImageType::Pointer > m_Images;
    std::vector< MaskType::Pointer > m_Masks;
    std::vector<typename Transform::Pointer> m_Transforms;
    std::vector<typename Transform::Pointer> m_InverseTransforms;
    bool m_VerboseMode;
    int m_NumberOfImages;
    int m_ReferenceImage;

    itk::Array<double> m_ScaleX;

    int m_ReferenceStack;
    int m_ReferenceSlice;

    int m_MaxLoop;

    std::vector<std::vector<double> > m_BestError;

    bool m_VerboseDbg;
    double m_CurrentError;

    bool m_UseSliceExclusion;

};


}

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkMotionCorrectionByIntersection.txx"
#endif

#endif // BTKMOTIONCORRECTIONBYINTERSECTION_H
