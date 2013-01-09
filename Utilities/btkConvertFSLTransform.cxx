/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 23/08/2012
  Author(s): Youssef Taleb, Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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

/* BTK */

#include "btkImageHelper.h"
#include "btkIOTransformHelper.h"

/* OTHERS */
#include "iostream"
#include "vector"
#include <tclap/CmdLine.h>
#include "fstream"

int main(int argc, char * argv[])
{
    /* Typedefs */
    typedef itk::Image<float, 3>                              ImageType;
    typedef itk::MatrixOffsetTransformBase< double,3 >        TransformType;
    typedef itk::Matrix<double,3,3>                           MatrixType;
    typedef itk::Matrix<double,4,4>                           TransformMatrixType;
    //TCLAP Commands for arguments

    TCLAP::CmdLine cmd("Convert a FSL matrix to a ITK Affine Transform, specify --inverse do the opposite", ' ', "Unversioned");
    TCLAP::ValueArg<std::string> FSLMatArg("m","matrix","Affine FSL matrix file (.mat)",true,"","string",cmd);
    TCLAP::ValueArg<std::string> inputImArg("i","input","input image (image on which transformation may be applied)",true,"string","",cmd);
    TCLAP::ValueArg<std::string> refArg("r","reference","Reference image",true,"","string",cmd);
    TCLAP::ValueArg<std::string> ITKtranArg  ("t","transformation","ITK Affine transformation file name (.txt)",true,"","string",cmd);
    TCLAP::SwitchArg invArg  ("","inverse","ITK to FSL",cmd);

    std::string   FSLFileName;
    std::fstream  FSLFile;
    std::string   ITKFileName;
    std::string   inputImageFile;
    std::string   refImageFile;

    // Parse the argv array.
    cmd.parse( argc, argv );

    FSLFileName=FSLMatArg.getValue();
    inputImageFile=inputImArg.getValue();
    ITKFileName=ITKtranArg.getValue();
    refImageFile=refArg.getValue();

    try
    {
        // Read reference image
        ImageType::Pointer refImage = ImageType::New();
        refImage=btk::ImageHelper<ImageType>::ReadImage(refImageFile);

        // Read input image
        ImageType::Pointer inputImage = ImageType::New();
        inputImage=btk::ImageHelper<ImageType>::ReadImage(inputImageFile);

        TransformMatrixType affine_fsl,m_itk;
        TransformType::Pointer itk_transform=TransformType::New();

        // FSL to ITK
        if(!invArg.getValue())
        {
            // Read FSL transform from file
            FSLFile.open(FSLFileName.c_str(), std::ios::in);

            for (int i = 0; i < 4; i++)
            {
                for (int j = 0; j < 4; j++)
                {
                    {
                        FSLFile >> affine_fsl(i,j)  ;
                    }
                }
            }
            FSLFile.close();
        }
        // ITK TO FSL
        else
        {
            // Read ITK transform from file and store it into 4*4 matrix
                itk_transform = btk::IOTransformHelper<TransformType>::ReadTransform(ITKFileName);
            m_itk.SetIdentity();
            m_itk.GetVnlMatrix().update((itk_transform->GetMatrix()).GetVnlMatrix());
            m_itk.GetVnlMatrix().set_column(3, (itk_transform->GetTranslation()).GetVnlVector());

            for (int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    if((i==2 & j<2) ||(i<2 &j==2))
                    {
                        m_itk(i,j)=-m_itk(i,j);
                    }

                }
            }
        }

        TransformMatrixType m_ref, m_mov;
        // Creation of reference image RAS Physical Space Matrix

        vnl_matrix<double> m_dir_ref, m_ras_matrix_ref;
        vnl_diag_matrix<double> m_scale_ref, m_lps_to_ras;
        vnl_vector<double> v_origin_ref, v_ras_offset_ref;

        // Compute the matrix
        m_dir_ref = refImage->GetDirection().GetVnlMatrix();
        m_scale_ref.set(refImage->GetSpacing().GetVnlVector());
        m_lps_to_ras.set(vnl_vector<double>(3, 1.0));
        m_lps_to_ras[0] = -1;
        m_lps_to_ras[1] = -1;
        m_ras_matrix_ref = m_lps_to_ras * m_dir_ref * m_scale_ref;

        // Compute the vector
        v_origin_ref = refImage->GetOrigin().GetVnlVector();
        v_ras_offset_ref = m_lps_to_ras * v_origin_ref;


        vnl_vector<double> refcol(4, 1.0);
        refcol.update(v_ras_offset_ref);
        m_ref.SetIdentity();
        m_ref.GetVnlMatrix().update(m_ras_matrix_ref);
        m_ref.GetVnlMatrix().set_column(3, refcol);


        // Creation of moving image RAS Physical Space Matrix
        vnl_matrix<double> m_dir, m_ras_matrix;
        vnl_diag_matrix<double> m_scale;
        vnl_vector<double> v_origin, v_ras_offset;

        // Compute the matrix
        m_dir = inputImage->GetDirection().GetVnlMatrix();
        m_scale.set(inputImage->GetSpacing().GetVnlVector());
        m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

        // Compute the vector
        v_origin = inputImage->GetOrigin().GetVnlVector();
        v_ras_offset = m_lps_to_ras * v_origin;

        vnl_vector<double> incol(4, 1.0);
        incol.update(v_ras_offset);
        m_mov.SetIdentity();
        m_mov.GetVnlMatrix().update(m_ras_matrix);
        m_mov.GetVnlMatrix().set_column(3, incol);

        TransformMatrixType  m_spcref, m_spcmov, m_swpref, m_swpmov, m_out;

        // Set the ref/mov matrices
        m_ref = m_ref.GetVnlMatrix();
        m_mov = m_mov.GetVnlMatrix();

        // Set the swap matrices
        m_swpref.SetIdentity();

        m_swpref(0,0) = -1.0;
        m_swpref(0,3) = (refImage->GetBufferedRegion().GetSize(0) - 1) * refImage->GetSpacing()[0];

        m_swpmov.SetIdentity();
        m_swpmov(0,0) = -1.0;
        m_swpmov(0,3) = (inputImage->GetBufferedRegion().GetSize(0) - 1) * inputImage->GetSpacing()[0];

        // Set the spacing matrices
        m_spcref.SetIdentity();
        m_spcmov.SetIdentity();
        for(size_t i = 0; i < 3; i++)
        {
            m_spcref(i,i) = refImage->GetSpacing()[i];
            m_spcmov(i,i) = inputImage->GetSpacing()[i];
        }

        // Compute the output matrix

        //FSL to ITK
        if(!invArg.getValue())
        {
            // M_OUT=MOVING_IMAGE_DIRECTION * (MOVING_IMAGE_SPACING)^-1 * MOVING_IMAGE_SWAP * FSL_MAT * REF_IMAGE_SWAP * REF_IMAGE_SPACING * (REF_IMAGE_DIRECTION)^-1
            m_out = m_mov*m_spcmov.GetInverse();
            m_out = m_out*m_swpmov;
            m_out = m_out*affine_fsl.GetInverse();
            m_out = m_out*m_swpref;
            m_out = m_out*m_spcref;
            m_out = m_out*m_ref.GetInverse();
        }
        //ITK TO FSL
        else
        {

            // M_OUT= (MOVING_IMAGE_SWAP)^-1 * MOVING_IMAGE_SPACING * (MOVING_IMAGE_DIRECTION)^-1 * ITK_MAT * REF_IMAGE_DIRECTION * (REF_IMAGE_SPACING)^-1 * (REF_IMAGE_SWAP)^-1
            m_out = m_swpmov.GetInverse();
            m_out = m_out*m_spcmov;
            m_out = m_out*m_mov.GetInverse();
            m_out = m_out*m_itk;
            m_out = m_out*m_ref;
            m_out = m_out*m_spcref.GetInverse();
            m_out = m_out*m_swpref.GetInverse();
            m_out = m_out.GetInverse();

        }

        //FSL To ITK
        if(!invArg.getValue())
        {

            // Compute transform
            TransformType::Pointer transform = TransformType::New();
            TransformType::ParametersType parameters(12);

            for (int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    if((i==2 & j<2) ||(i<2 &j==2))
                    {
                        parameters[3*i+j]=-m_out(i,j);
                    }

                    else
                    {
                        parameters[3*i+j]=m_out(i,j);
                    }
                }
            }

            for(int i=9;i<12;i++)
            {
                parameters[i]=m_out(i-9,3);
            }

            transform->SetParameters(parameters);
            ImageType::PointType center;
            center[0]=0;
            center[1]=0;
            center[2]=0;
            transform->SetCenter(center);
            btk::IOTransformHelper<TransformType>::WriteTransform(transform,ITKFileName);
        }
        //ITK to FSL
        else
        {
            FSLFile.open(FSLFileName.c_str(), std::ios::out);

            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 4; j++)
                    FSLFile << m_out[i][j] << (j < 3 ? " " : "\n");

            FSLFile.close();
        }

    }
    catch(itk::ExceptionObject &error)
    {
        std::cout << "ITK error: " << error << std::endl;
        return EXIT_FAILURE;
    }
    catch(std::string &message)
    {
        std::cout << "Error: " << message << std::endl;
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;

}
