/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 09/01/2013
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


/* ITK */
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkImageFileReader.h"
#include "itkAffineTransform.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkTransformFactory.h"
#include "itkByteSwapper.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_det.h"
#include "vnl/vnl_inverse.h"
#include "vnl/algo/vnl_real_eigensystem.h"
#include "vnl/algo/vnl_qr.h"

/* BTK */

#include "btkImageHelper.h"
#include "btkIOTransformHelper.h"
#include "btkOrientedRASImage.h"
#include "btkMacro.h"

/* OTHERS */
#include "iostream"
#include "vector"
#include <tclap/CmdLine.h>
#include "fstream"

int main(int argc, char * argv[])
{
    typedef btk::OrientedRASImage<double, 3> ImageType;
    typedef vnl_matrix_fixed<double, 4, 4> MatrixType;
    typedef std::vector<MatrixType> MatrixStack;
    typedef itk::MatrixOffsetTransformBase<double, 3, 3> MatrixOffsetType;

    TCLAP::CmdLine cmd("Convert a FSL matrix to a ITK Affine Transform, specify --inverse do the opposite", ' ', "Unversioned");
    TCLAP::ValueArg<std::string> FSLMatArg("m","matrix","Affine FSL matrix file (.mat)",true,"","string",cmd);
    TCLAP::ValueArg<std::string> inputImArg("i","input","input image (image on which transformation may be applied)",true,"string","",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","Reference image",true,"","string",cmd);
    TCLAP::ValueArg<std::string> ITKtranArg  ("t","transformation","ITK Affine transformation file name (.txt)",true,"","string",cmd);
    TCLAP::SwitchArg invArg  ("","inverse","ITK to FSL",cmd);

    std::string   FSLFileName;
    std::string   ITKFileName;
    std::string   inputImageFile;
    std::string   refImageFile;

    // Parse the argv array.
    cmd.parse( argc, argv );

    FSLFileName=FSLMatArg.getValue();
    inputImageFile=inputImArg.getValue();
    ITKFileName=ITKtranArg.getValue();
    refImageFile=refArg.getValue();

    ImageType::Pointer mov_image = ImageType::New();
    mov_image = btk::ImageHelper<ImageType>::ReadImage(inputImageFile);

    ImageType::Pointer ref_image = ImageType::New();
    ref_image = btk::ImageHelper<ImageType>::ReadImage(refImageFile);

    // Set up the matrix stack
    std::vector<MatrixType> vmat;
    MatrixType fsl_mat, spc_ref_mat, spc_mov_mat, swp_ref_mat, swp_mov_mat, ref_RAS, mov_RAS, out_mat;


    // Set the ref/mov matrices
    ref_RAS = ref_image->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();
    mov_RAS = mov_image->GetVoxelSpaceToRASPhysicalSpaceMatrix().GetVnlMatrix();

    // Set the swap matrices
    swp_ref_mat.set_identity();
    if(vnl_det(ref_RAS) > 0)
    {
        swp_ref_mat(0,0) = -1.0;
        swp_ref_mat(0,3) = (ref_image->GetBufferedRegion().GetSize(0) - 1) * ref_image->GetSpacing()[0];
    }

    swp_mov_mat.set_identity();
    if(vnl_det(mov_RAS) > 0)
    {
        swp_mov_mat(0,0) = -1.0;
        swp_mov_mat(0,3) = (mov_image->GetBufferedRegion().GetSize(0) - 1) * mov_image->GetSpacing()[0];
    }

    // Set the spacing matrices
    spc_ref_mat.set_identity();
    spc_mov_mat.set_identity();
    for(size_t i = 0; i < 3; i++)
    {
        spc_ref_mat(i,i) = ref_image->GetSpacing()[i];
        spc_mov_mat(i,i) = mov_image->GetSpacing()[i];
    }

    // FSL to ITK
    if(!invArg.getValue())
    {
        MatrixType mat;

        std::ifstream fin(FSLFileName.c_str());
        for(size_t i = 0; i < 4; i++)
            for(size_t j = 0; j < 4; j++)
                if(fin.good())
                {
                    fin >> mat[i][j];
                }
                else
                {
                    throw "Unable to read matrix";
                }
        fin.close();

        vmat.push_back(mat);

        fsl_mat = vmat.back();

        //------------------

        out_mat =
                mov_RAS * vnl_inverse(spc_mov_mat) * swp_mov_mat *
                vnl_inverse(fsl_mat) *
                swp_ref_mat * spc_ref_mat * vnl_inverse(ref_RAS);


    }
    //ITK to FSL
    else
    {
        typedef itk::MatrixOffsetTransformBase<double, 3, 3> MatrixOffsetType;
        typedef itk::AffineTransform<double, 3> AffTran;
        itk::TransformFactory<MatrixOffsetType>::RegisterTransform();
        itk::TransformFactory<AffTran>::RegisterTransform();

        // TODO: Check if we can use btk::IOTransformHelper::Read
        itk::TransformFileReader::Pointer fltReader = itk::TransformFileReader::New();
        fltReader->SetFileName(ITKFileName);
        fltReader->Update();

        itk::TransformBase *base = fltReader->GetTransformList()->front();
        MatrixOffsetType *transform = dynamic_cast<MatrixOffsetType *>(base);

        MatrixType mat;
        mat.set_identity();
        if(transform)
        {
            for(size_t r = 0; r < 3; r++)
            {
                for(size_t c = 0; c < 3; c++)
                {
                    mat(r,c) = transform->GetMatrix()(r,c);
                }
                mat(r,3) = transform->GetOffset()[r];
            }
            mat(2,0) *= -1; mat(2,1) *= -1;
            mat(0,2) *= -1; mat(1,2) *= -1;
            mat(0,3) *= -1; mat(1,3) *= -1;
            vmat.push_back(mat);
        }

        //-------------
        fsl_mat = vmat.back();



        out_mat = vnl_inverse(vnl_inverse(swp_mov_mat) * spc_mov_mat* vnl_inverse(mov_RAS) *
                            fsl_mat *
                            ref_RAS*
                            vnl_inverse(spc_ref_mat)*
                            vnl_inverse(swp_ref_mat));


    }

    // Put it on the stack
    vmat.pop_back();
    vmat.push_back(out_mat);


    //-------
    //Write :
    if(invArg.getValue())
    {
        MatrixType mat = vmat.back();

        std::ofstream fout(FSLFileName.c_str());
        for(size_t i = 0; i < 4; i++)
            for(size_t j = 0; j < 4; j++)
                fout << mat[i][j] << (j < 3 ? " " : "\n");

        fout.close();
    }
    else
    {
        // Get the current matrix
        MatrixType mat = vmat.back();

        // Flip the entries that must be flipped
        mat(2,0) *= -1; mat(2,1) *= -1;
        mat(0,2) *= -1; mat(1,2) *= -1;
        mat(0,3) *= -1; mat(1,3) *= -1;

        // Create an ITK affine transform
        MatrixOffsetType::Pointer transform_to_write = MatrixOffsetType::New();

        // Populate its matrix
        MatrixOffsetType::MatrixType amat = transform_to_write->GetMatrix();
        MatrixOffsetType::OffsetType aoff = transform_to_write->GetOffset();

        for(size_t r = 0; r < 3; r++)
        {
            for(size_t c = 0; c < 3; c++)
            {
                if(mat(r,c) < 1.0e-06 && mat(r,c) > -1.0e-06)
                {
                    mat(r,c) = 0.0;
                }
                amat(r,c) = mat(r,c);
            }
            aoff[r] = mat(r,3);
        }

        transform_to_write->SetMatrix(amat);
        transform_to_write->SetOffset(aoff);

        // Write the transform

        btk::IOTransformHelper< MatrixOffsetType>::WriteTransform(transform_to_write, ITKFileName);

    }



}
