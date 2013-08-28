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

#ifndef BTK_ORIENTEDRASIMAGE_H
#define BTK_ORIENTEDRASIMAGE_H

#include "itkImage.h"

namespace btk
{

/** @class OrientedRASImage
 * @brief
 * Oriented image with RAS physical coordinates, as opposite of itk::Image with LPS physical coordinates.
 * itk::Image is the base class.
 * @author Marc Schweitzer
 * @ingroup Images
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
        {
            return NeighborhoodAccessorFunctorType();
        }

        /** Return the NeighborhoodAccessor functor. This method is called by the
         * neighborhood iterators. */
        const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
        {
            return NeighborhoodAccessorFunctorType();
        }


        /** \brief Get the continuous index from a physical point
         * Returns true if the resulting index is within the image, false otherwise.
         * \sa Transform */
        template<class TCoordRep>
        bool TransformRASPhysicalPointToContinuousIndex(
                const itk::Point<TCoordRep, Dimension>& point,
                itk::ContinuousIndex<TCoordRep, Dimension>& index) const
        {
            itk::Point<TCoordRep, Dimension> point_LPS = point;
            point_LPS[0] = -point[0]; point_LPS[1] = -point[1];
            return Superclass::TransformPhysicalPointToContinuousIndex(point_LPS, index);
        }

        /** Get the index (discrete) from a physical point.
         * Floating point index results are truncated to integers.
         * Returns true if the resulting index is within the image, false otherwise
         * \sa Transform */
        template<class TCoordRep>
        bool TransformRASPhysicalPointToIndex(
                const itk::Point<TCoordRep, Dimension>& point,
                IndexType & index) const
        {
            itk::Point<TCoordRep, Dimension> point_LPS = point;
            point_LPS[0] = -point[0]; point_LPS[1] = -point[1];
            return Superclass::TransformPhysicalPointToIndex(point_LPS, index);
        }

        /** Get a physical point (in the space which
         * the origin and spacing infomation comes from)
         * from a continuous index (in the index space)
         * \sa Transform */
        template<class TCoordRep>
        void TransformContinuousIndexToRASPhysicalPoint(
                const itk::ContinuousIndex<TCoordRep, Dimension>& index,
                itk::Point<TCoordRep, Dimension>& point) const
        {
            Superclass::TransformContinuousIndexToPhysicalPoint(index, point);
            point[0] = -point[0];
            point[1] = -point[1];
        }

        /** Get a physical point (in the space which
         * the origin and spacing infomation comes from)
         * from a discrete index (in the index space)
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
         * \sa Image */
        template<class TCoordRep>
        void TransformLocalVectorToRASPhysicalVector(
                const itk::FixedArray<TCoordRep, Dimension> & inputGradient,
                itk::FixedArray<TCoordRep, Dimension> & outputGradient ) const
        {
            Superclass::TransformLocalVectorToPhysicalVector(inputGradient, outputGradient);
            outputGradient[0] = -outputGradient[0];
            outputGradient[1] = -outputGradient[1];
        }

        /** Get a matrix that maps points voxel coordinates to RAS coordinates */
        TransformMatrixType GetVoxelSpaceToRASPhysicalSpaceMatrix();

        /** Set a matrix that maps points voxel coordinates to RAS coordinates */
        void SetVoxelSpaceToRASPhysicalSpaceMatrix(vnl_matrix<double> mat);

        /** Get a matrix that maps points in the x * spacing + origin space to the RAS space */
        TransformMatrixType GetSpacingOriginPhysicalSpaceToRASPhysicalSpaceMatrix();



    protected:
        /** Constructor */
        OrientedRASImage();
        /** Destructor */
        virtual ~OrientedRASImage();
        /** Print Self */
        void PrintSelf(std::ostream& os, itk::Indent indent) const;

    private:
        OrientedRASImage(const Self&); //purposely not implemented
        void operator=(const Self&); //purposely not implemented
};

} //namespace btk

#ifndef ITK_MANUAL_INSTANTIATION
#include "btkOrientedRASImage.txx"
#endif


#endif // BTKORIENTEDRASIMAGE_H
