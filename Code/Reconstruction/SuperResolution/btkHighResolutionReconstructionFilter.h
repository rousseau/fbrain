/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 22/03/2012
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
#ifndef __BTK_HIGHRESOLUTIONRECONSTRUCTIONFILTER_H__
#define __BTK_HIGHRESOLUTIONRECONSTRUCTIONFILTER_H__

/* ITK */
#include "itkImage.h"
#include "itkImageMaskSpatialObject.h"
#include "itkIdentityTransform.h"
#include "itkTransformFactory.h"
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

/* BTK */
#include "btkMacro.h"
#include "btkSliceBySliceTransform.h"
//Typedefs:
#include "btkSuperResolutionType.h"


/* VNL */
#include "vnl/vnl_sparse_matrix.h"

/* OTHERS */
#include "iostream"

namespace btk
{
class HighResolutionReconstructionFilter
{
    //TODO: Maybe we should create sub-class of this one for doing several type of reconstruction !
    // ex : HighResolutionIterativeBackProjecton, HighResolutionSR...



public:


    HighResolutionReconstructionFilter();
    ~HighResolutionReconstructionFilter();

    virtual void Update() = 0;

    virtual itkImage::Pointer GetOutput()
    {
        return m_OutputHRImage;
    }

    // GETTER/SETTER :


    //TODO:Remove useless GET/SET or getter/setter for parameters which are compute and used  only in this class
    btkGetMacro(TransformsLR,std::vector< itkTransformBase::Pointer >);
    btkSetMacro(TransformsLR,std::vector< itkTransformBase::Pointer >);

    btkGetMacro(InverseTransformsLR,std::vector< itkTransformBase::Pointer >);
    btkSetMacro(InverseTransformsLR,std::vector< itkTransformBase::Pointer >);

    btkGetMacro(TransformsLRAffine,std::vector< itkAffineTransform::Pointer >);
    btkSetMacro(TransformsLRAffine,std::vector< itkAffineTransform::Pointer >);

    btkGetMacro(InverseTransformsLRAffine,std::vector< itkAffineTransform::Pointer >);
    btkSetMacro(InverseTransformsLRAffine,std::vector< itkAffineTransform::Pointer >);

    btkGetMacro(TransformsLRSbS,std::vector< btkSliceBySliceTransform::Pointer >);
    btkSetMacro(TransformsLRSbS,std::vector< btkSliceBySliceTransform::Pointer >);

    btkGetMacro(InverseTransformsLRSbS,std::vector< btkSliceBySliceTransform::Pointer >);
    btkSetMacro(InverseTransformsLRSbS,std::vector< btkSliceBySliceTransform::Pointer >);

    btkGetMacro(ImagesLR,std::vector< itkImage::Pointer > );
    btkSetMacro(ImagesLR,std::vector< itkImage::Pointer > );

    btkGetMacro(SimulatedImagesLR,std::vector< itkImage::Pointer > );
    btkSetMacro(SimulatedImagesLR,std::vector< itkImage::Pointer > );


    btkGetMacro(ImagesMaskLR,std::vector< itkImageMask::Pointer >  );
    btkSetMacro(ImagesMaskLR,std::vector< itkImageMask::Pointer >  );

    btkGetMacro(MasksLR,std::vector< itkMask::Pointer >);
    btkSetMacro(MasksLR,std::vector< itkMask::Pointer >);

    btkGetMacro(ImageHR, itkImage::Pointer);
    btkSetMacro(ImageHR, itkImage::Pointer);

    btkGetMacro(ImageMaskHR, itkImage::Pointer);
    btkSetMacro(ImageMaskHR, itkImage::Pointer);

    btkGetMacro(ReferenceImage, itkImage::Pointer);
    btkSetMacro(ReferenceImage, itkImage::Pointer);


    btkSetMacro(H,vnl_sparse_matrix< float >);
    btkGetMacro(H,vnl_sparse_matrix< float >);

    btkSetMacro(X,vnl_vector< float >);
    btkGetMacro(X,vnl_vector< float >);


    btkSetMacro(Y,vnl_vector< float >);
    btkGetMacro(Y,vnl_vector< float >);

    btkGetMacro(Offset,std::vector< unsigned int >);
    btkSetMacro(Offset,std::vector< unsigned int >);

    btkGetMacro(PSF, std::vector< itkImage::Pointer >);
    btkSetMacro(PSF, std::vector< itkImage::Pointer >);

    btkGetMacro(PaddingValue,float);
    btkSetMacro(PaddingValue,float);

    btkGetMacro(InterpolationOrderPSF,int);
    btkSetMacro(InterpolationOrderPSF,int);

    btkGetMacro(InterpolationOrderIBP,int);
    btkSetMacro(InterpolationOrderIBP,int);

    btkGetMacro(PsfType,int);
    btkSetMacro(PsfType,int);

    btkGetMacro(Nloops,int);
    btkSetMacro(Nloops,int);

    btkGetMacro(TransformType,TRANSFORMATION_TYPE);
    btkSetMacro(TransformType,TRANSFORMATION_TYPE);


protected:
    // accessible for subclasses
    std::vector< itkImage::Pointer >         m_ImagesLR;
    std::vector< itkImage::Pointer >         m_SimulatedImagesLR;
    std::vector< itkTransformBase::Pointer  >   m_TransformsLR;
    std::vector< itkTransformBase::Pointer >    m_InverseTransformsLR;

    std::vector< itkAffineTransform::Pointer  >   m_TransformsLRAffine;
    std::vector< itkAffineTransform::Pointer >    m_InverseTransformsLRAffine;

    std::vector< btkSliceBySliceTransform::Pointer  >   m_TransformsLRSbS;
    std::vector< btkSliceBySliceTransform::Pointer >    m_InverseTransformsLRSbS;

    std::vector< itkImageMask::Pointer >     m_ImagesMaskLR;
    std::vector< itkMask::Pointer >          m_MasksLR;
    itkImage::Pointer                        m_ImageHR;
    itkImage::Pointer                        m_CurrentImageHR;
    itkImage::Pointer                        m_OutputHRImage;
    itkImage::Pointer                        m_ImageMaskHR;
    itkMask::Pointer                         m_MaskHR;
    itkImage::Pointer                        m_ReferenceImage;

    TRANSFORMATION_TYPE                      m_TransformType;

    vnl_sparse_matrix< float >               m_H;
    vnl_vector< float >                      m_X;
    vnl_vector< float >                      m_Y;
    std::vector< unsigned int >              m_Offset;
    std::vector< itkImage::Pointer >         m_PSF;
    float                                    m_PaddingValue;
    int                                      m_InterpolationOrderPSF;
    int                                      m_InterpolationOrderIBP;
    int                                      m_PsfType;
    int                                      m_Nloops;

private:



};
}

#endif
