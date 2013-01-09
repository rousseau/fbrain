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

#ifndef BTKORIENTEDRASIMAGE_H
#define BTKORIENTEDRASIMAGE_H

#include "itkImage.h"

namespace btk
{

/**
 * Oriented image with RAS physical coordinates (as opposed to LPS)
 */
template <class TPixel, unsigned int Dimension>
class ITK_EXPORT OrientedRASImage : public itk::Image<TPixel, Dimension>
{
    public:
        /** Standard class typedefs */
        typedef OrientedRASImage               Self;
        typedef itk::Image<TPixel, Dimension>  Superclass;
        typedef itk::SmartPointer<Self>  Pointer;
        typedef itk::SmartPointer<const Self>  ConstPointer;
        typedef itk::WeakPointer<const Self>  ConstWeakPointer;
        typedef itk::Matrix<double, Dimension+1, Dimension+1> TransformMatrixType;

        /** Method for creation through the object factory. */
        itkNewMacro(Self);

        /** Run-time type information (and related methods). */
        itkTypeMacro(OrientedRASImage, Image);

        /** Index typedef support. An index is used to access pixel values. */
        typedef typename Superclass::IndexType  IndexType;

        /** Direction typedef support. The direction cosines of the image. */
        typedef typename Superclass::DirectionType  DirectionType;

        /** Spacing typedef support.  Spacing holds the size of a pixel.  The
   * spacing is the geometric distance between image samples. */
        typedef typename Superclass::SpacingType SpacingType;

        typedef typename Superclass::AccessorType        AccessorType;
        typedef typename Superclass::AccessorFunctorType AccessorFunctorType;
        typedef typename Superclass::IOPixelType         IOPixelType;

        /** Tyepdef for the functor used to access a neighborhood of pixel pointers.*/
        typedef itk::NeighborhoodAccessorFunctor< Self >
        NeighborhoodAccessorFunctorType;

        /** Return the NeighborhoodAccessor functor. This method is called by the
   * neighborhood iterators. */
        NeighborhoodAccessorFunctorType GetNeighborhoodAccessor()
        { return NeighborhoodAccessorFunctorType(); }

        /** Return the NeighborhoodAccessor functor. This method is called by the
   * neighborhood iterators. */
        const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
        { return NeighborhoodAccessorFunctorType(); }


        /** \brief Get the continuous index from a physical point
   *
   * Returns true if the resulting index is within the image, false otherwise.
   * \sa Transform */
        template<class TCoordRep>
        bool TransformRASPhysicalPointToContinuousIndex(
                const itk::Point<TCoordRep, Dimension>& point,
                itk::ContinuousIndex<TCoordRep, Dimension>& index   ) const
        {
            itk::Point<TCoordRep, Dimension> p_lps = point;
            p_lps[0] = -point[0]; p_lps[1] = -point[1];
            return Superclass::TransformPhysicalPointToContinuousIndex(p_lps, index);
        }

        /** Get the index (discrete) from a physical point.
   * Floating point index results are truncated to integers.
   * Returns true if the resulting index is within the image, false otherwise
   * \sa Transform */
        template<class TCoordRep>
        bool TransformRASPhysicalPointToIndex(
                const itk::Point<TCoordRep, Dimension>& point,
                IndexType & index                                ) const
        {
            itk::Point<TCoordRep, Dimension> p_lps = point;
            p_lps[0] = -point[0]; p_lps[1] = -point[1];
            return Superclass::TransformPhysicalPointToIndex(p_lps, index);
        }

        /** Get a physical point (in the space which
   * the origin and spacing infomation comes from)
   * from a continuous index (in the index space)
   * \sa Transform */
        template<class TCoordRep>
        void TransformContinuousIndexToRASPhysicalPoint(
                const itk::ContinuousIndex<TCoordRep, Dimension>& index,
                itk::Point<TCoordRep, Dimension>& point        ) const
        {
            Superclass::TransformContinuousIndexToPhysicalPoint(index, point);
            point[0] = -point[0];
            point[1] = -point[1];
        }

        /** Get a physical point (in the space which
   * the origin and spacing infomation comes from)
   * from a discrete index (in the index space)
   *
   * \sa Transform */
        template<class TCoordRep>
        void TransformIndexToRASPhysicalPoint(
                const IndexType & index,
                itk::Point<TCoordRep, Dimension>& point ) const
        {
            Superclass::TransformIndexToPhysicalPoint(index, point);
            point[0] = -point[0];
            point[1] = -point[1];
        }

        /** Take a vector or covariant vector that has been computed in the
   * coordinate system parallel to the image grid and rotate it by the
   * direction cosines in order to get it in terms of the coordinate system of
   * the image acquisition device.  This implementation in the Image
   * multiply the array (vector or covariant vector) by the matrix of Direction
   * Cosines. The arguments of the method are of type FixedArray to make
   * possible to use this method with both Vector and CovariantVector.
   * The Method is implemented differently in the itk::Image.
   *
   * \sa Image
   */
        template<class TCoordRep>
        void TransformLocalVectorToRASPhysicalVector(
                const itk::FixedArray<TCoordRep, Dimension> & inputGradient,
                itk::FixedArray<TCoordRep, Dimension> & outputGradient ) const
        {
            Superclass::TransformLocalVectorToPhysicalVector(inputGradient, outputGradient);
            outputGradient[0] = -outputGradient[0];
            outputGradient[1] = -outputGradient[1];
        }

        /**
   * Get a matrix that maps points voxel coordinates to RAS coordinates
   */
        TransformMatrixType GetVoxelSpaceToRASPhysicalSpaceMatrix()
        {
            // Generate intermediate terms
            vnl_matrix<double> m_dir, m_ras_matrix;
            vnl_diag_matrix<double> m_scale, m_lps_to_ras;
            vnl_vector<double> v_origin, v_ras_offset;

            // Compute the matrix
            m_dir = this->GetDirection().GetVnlMatrix();
            m_scale.set(this->GetSpacing().GetVnlVector());
            m_lps_to_ras.set(vnl_vector<double>(Dimension, 1.0));
            m_lps_to_ras[0] = -1;
            m_lps_to_ras[1] = -1;
            m_ras_matrix = m_lps_to_ras * m_dir * m_scale;

            // Compute the vector
            v_origin = this->GetOrigin().GetVnlVector();
            v_ras_offset = m_lps_to_ras * v_origin;

            // Create the larger matrix
            TransformMatrixType mat;
            vnl_vector<double> vcol(Dimension+1, 1.0);
            vcol.update(v_ras_offset);
            mat.SetIdentity();
            mat.GetVnlMatrix().update(m_ras_matrix);
            mat.GetVnlMatrix().set_column(Dimension, vcol);

            return mat;
        };

        /**
   * Set a matrix that maps points voxel coordinates to RAS coordinates
   */
        void SetVoxelSpaceToRASPhysicalSpaceMatrix(vnl_matrix<double> mat)
        {
            // Generate intermediate terms
            vnl_matrix<double> m_dir, m_ras_matrix, m_dist;
            vnl_diag_matrix<double> m_ras_to_lps, m_scale;
            vnl_vector<double> v_origin ;
            vnl_vector<double> m_spacing(Dimension, 0.0);

            // Get the dim x dim submatrix from mat
            vnl_matrix<double> smat(Dimension,Dimension,0.0);
            for (size_t i=0; i< Dimension; i++)
                for (size_t j=0; j< Dimension; j++)
                    smat[i][j] = mat[i][j];
            //smat = mat.get_n_rows(0, Dimension).get_n_columns(0, Dimension);
            // Get the origin
            m_ras_to_lps.set(vnl_vector<double>(Dimension, 1.0));
            m_ras_to_lps[0] = -1;
            m_ras_to_lps[1] = -1;
            vnl_vector<double> v_ras_offset(Dimension,0.0);
            v_ras_offset.fill(0.0);
            for (size_t i=0; i< Dimension; i++)
                v_ras_offset[i] = mat[i][Dimension];
            v_origin = m_ras_to_lps * v_ras_offset;

            // Get the Spacing
            // First, create a matrix of the form [1 0 0; 0 1 0; 0 0 1; 0 0 0] to get distances between consecutive voxels
            // along each axis. When RAS mat will be applied to this matrix, we'll have 3 distance vectors
            vnl_diag_matrix<double> offsetmat(Dimension+1, Dimension);
            offsetmat.fill(0.0);
            for (size_t i=0; i < Dimension+1; i++)
                offsetmat[i]=1.0;
            m_dist = mat * offsetmat;
            // Then compute magnitude of the distance vectors, that's our spacing
            for (size_t i=0; i< Dimension; i++)
            {
                vnl_vector<double> distcol(m_dist.get_column(i));
                m_spacing[i] = distcol.magnitude();
            }
            m_scale.set(m_spacing);

            // Get the direction
            m_scale.invert_in_place();
            m_dir = m_ras_to_lps * smat * m_scale;

            // Set everything
            itk::Matrix<double, Dimension, Dimension> dir;
            dir.SetIdentity();
            for (size_t i=0; i< Dimension; i++)
                for (size_t j=0; j< Dimension; j++)
                    dir[i][j] = m_dir[i][j];
            this->SetDirection(dir);
            double origin[Dimension];
            for (size_t i=0; i< Dimension; i++)
                origin[i] = v_origin[i];
            this->SetOrigin(origin);
            double spacing[Dimension];
            for (size_t i=0; i< Dimension; i++)
                spacing[i] = m_spacing[i];
            this->SetSpacing(spacing);

        };

        /**
   * Get a matrix that maps points in the x * spacing + origin space to
   * the RAS space
   */
        TransformMatrixType GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix()
        {
            TransformMatrixType mat;
            mat.SetIdentity();

            for(size_t i = 0; i < Dimension; i++)
            {
                double ras_flip = (i < 2) ? -1 : 1;
                mat[i][Dimension] = ras_flip * this->GetOrigin()[i];
                for(size_t j = 0; j < Dimension; j++)
                {
                    mat[i][j] = ras_flip * this->GetDirection()(i,j) * this->GetSpacing()[i];
                    mat[i][Dimension] -= ras_flip * this->GetDirection()(i,j) * this->GetOrigin()[i];
                }
            }

            return mat;
        }



    protected:
        OrientedRASImage() {};
        virtual ~OrientedRASImage() {};

    private:
        OrientedRASImage(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented
};

} //namespace itk



#endif // BTKORIENTEDRASIMAGE_H
