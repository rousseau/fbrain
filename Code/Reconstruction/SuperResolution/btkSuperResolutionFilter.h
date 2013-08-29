/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 16/03/2012
  Author(s): Schweitzer Marc (marc.schweitzer@unistra.fr)

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

#ifndef __BTK_SUPERRESOLUTIONFILTER_H__
#define __BTK_SUPERRESOLUTIONFILTER_H__

/* ITK */
#include "itkObject.h"
#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkTransform.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_conjugate_gradient.h"
#include "vnl/vnl_matops.h"

/* BTK */
#include "btkMacro.h"
#include "btkSRHMatrixComputation.hxx"
#include "btkPSF.h"
//#include "btkSuperResolutionCostFunctionITKWrapper.h"
//NOTE : Use ITK Wrapper when you want to use btkSuperResolutioCostFunction with an
// itk optimizer, and use VNL Wrapper when you want to use a VNL Optimizer
#include "btkSuperResolutionCostFunctionVNLWrapper.hxx"

/* others */

#include "vector"

namespace btk
{
/**
 * @class SuperResolutionFilter
 * @brief SuperResolutionFilter is a filter that handle the super-resolution process.
 * - First it call the H matrix filter (for compute H)
 * - Then it call the cost function that compute Min(f(y-Hx) + lambda g(x))
 * - The cost function is connected to the vnl optimizer
 * - Finally it generate the output (super-resolution image)
 * - Optionaly low resolution simulated image can be generated (for testing)
 *
 * The implemented method is similar to the one described in:
 *
 * F. Rousseau,  K. Kim,  C. Studholme,  M. Koob,  J.-L. Dietemann
 * On Super-Resolution for Fetal Brain MRI, Medical Image Computing and Computer
 * Assisted Intervention Pékin, Chine, pp. 355--362, Springer-Verlag, Lecture
 * Notes in Computer Science, Vol. 6362, doi:10.1007/978-3-642-15745-5, 2010
 * @ingroup SuperResolution
 * @author Marc Schweitzer
 */
class SuperResolutionFilter: public itk::Object
{
    public:
        /** Typedefs */
        typedef SuperResolutionFilter           Self;
        typedef itk::Object                     Superclass;
        typedef itk::SmartPointer< Self >       Pointer;
        typedef itk::SmartPointer< const Self > ConstPointer;

        typedef itk::Transform< double, 3 >     TransformType;

        typedef float                           PixelType;

        typedef float                           PrecisionType;

        typedef itk::Image< PixelType, 3 >          ImageType; // TODO : template ?

        typedef itk::Image< unsigned char, 3 >  ImageMaskType;

        typedef ImageType::RegionType           RegionType;

        typedef itk::ImageRegionIteratorWithIndex< ImageType >    Iterator;

        typedef btk::SRHMatrixComputation< ImageType >  H_Filter;

        typedef itk::ImageMaskSpatialObject< ImageType::ImageDimension > MaskType;

        typedef btk::SuperResolutionCostFunctionVNLWrapper< ImageType > VNLCostFunction;
        typedef itk::ImageRegionIteratorWithIndex< ImageType > itkIteratorWithIndex;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(btk::SuperResolutionFilter, btk::SuperResolutionFilter);


        void Update();

        /** Get/Set LR Images */
        btkSetMacro(Images, std::vector< ImageType::Pointer > );
        btkGetMacro(Images, std::vector< ImageType::Pointer >);

        /** Get/Set LR Masks */
        btkSetMacro(Masks, std::vector< ImageMaskType::Pointer > );
        btkGetMacro(Masks, std::vector< ImageMaskType::Pointer >);

        /** Get/Set Transforms */
        btkSetMacro(Transforms, std::vector< TransformType::Pointer >);
        btkGetMacro(Transforms, std::vector< TransformType::Pointer >);

        /** Get/Set Inverse Transforms */
        btkSetMacro(InverseTransforms, std::vector< TransformType::Pointer >);
        btkGetMacro(InverseTransforms, std::vector< TransformType::Pointer >);

        /** Get/Set Reference Image */
        btkSetMacro(ReferenceImage, ImageType::Pointer);
        btkGetMacro(ReferenceImage, ImageType::Pointer);

        /** Get Simulated Images */
        btkGetMacro(SimulatedImages, std::vector< ImageType::Pointer >);

        /** Get Output */
        btkGetMacro(Output,ImageType::Pointer);

        /** Set Lambda */
        btkSetMacro(Lambda,float);

        /** Use simulated images */

        void ComputeSimulatedImages(bool _b)
        {
            m_ComputeSimulations = _b;
        }

        /** Initialize must be called before Update() */
        virtual void Initialize();

        /** Add Transformation */
        void AddTransform(TransformType* _t)
        {
            if(m_Transforms.size() >= m_Images.size())
            {
                btkException("Size of transformations are greater than size of input images !");
            }

            m_Transforms.push_back(_t);
        }

        /** Add Inverse transformation */
        void AddInverseTransform(TransformType* _t)
        {
            if(m_InverseTransforms.size() >= m_Images.size())
            {
                btkException("Size of transformations are greater than size of input images !");
            }

            m_InverseTransforms.push_back(_t);
        }

        /** Set PSF to the H matrix filter  */
        void SetPSF(btk::PSF::Pointer _p)
        {
            m_H_Filter->SetPSF(_p);
        }


    protected:
        /** Simulate LR Images with the precalculated H */
        virtual void SimulateLRImages();
        /** Generate the output image */
        virtual void GenerateOutputData();
        /** Constructor */
        SuperResolutionFilter();
        /** Destructor */
        virtual ~SuperResolutionFilter();

    private:

        H_Filter::Pointer                       m_H_Filter;

        std::vector< TransformType::Pointer >   m_Transforms;

        std::vector< TransformType::Pointer >   m_InverseTransforms;

        std::vector< ImageType::Pointer >       m_Images;

        std::vector< RegionType >               m_Regions;

        std::vector< ImageType::Pointer >       m_SimulatedImages;

        std::vector< ImageMaskType::Pointer >   m_Masks;

        std::vector< MaskType::Pointer >        m_MasksObject;

        ImageType::Pointer                      m_ReferenceImage;

        vnl_vector< double >                    m_X;

        vnl_vector< PrecisionType >             m_Xfloat;

        ImageType::RegionType                   m_ReferenceRegion;

        btk::PSF::Pointer                       m_PSF;

        vnl_sparse_matrix< PrecisionType >*     m_H;

        vnl_vector< PrecisionType >*            m_Y;

        bool                                    m_ComputeSimulations;

        unsigned int                            m_NumberOfImages;

        ImageType::Pointer                      m_Output;

        float                                   m_Lambda;


};//end class

}//end namespace

#endif
