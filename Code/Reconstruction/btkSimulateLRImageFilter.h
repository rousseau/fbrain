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

#ifndef __BTK_SIMULATELRIMAGEFILTER_H__
#define __BTK_SIMULATELRIMAGEFILTER_H__

//ITK Includes
#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"
#include "itkStatisticsImageFilter.h"



//VNL includes
#include "vnl/vnl_sparse_matrix.h"



namespace btk
{
class SimulateLRImageFilter
{
public:
    typedef float PixelType;
    typedef itk::Image< PixelType, 3>         itkImage;
    typedef itk::Image< PixelType, 2>         SliceType;
    typedef itkImage::Pointer                 itkPointer;
    typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
    typedef itk::ImageRegionIterator< itkImage > itkIterator;
    typedef itk::ImageRegionIteratorWithIndex< itkImage > itkIteratorWithIndex;
    typedef itk::ContinuousIndex<double,3>     itkContinuousIndex;


    typedef itk::StatisticsImageFilter<itkImage>        itkStatisticsImageFilter;

    SimulateLRImageFilter();
    ~SimulateLRImageFilter();
    void Update();

    std::vector< itkPointer > GetOutput()
    {
        return m_simulatedOutputImages;
    }

    //GETTER/SETTER

    void SetLRImages(std::vector< itkPointer > images)
    {
        m_lRImages = images;
        m_simulatedOutputImages.resize(m_lRImages.size());
    }
    std::vector< itkPointer > GetLRImages()
    {
        return m_lRImages;
    }

    void setH(vnl_sparse_matrix<float> & H)
    {
        m_H = H;
    }
    vnl_sparse_matrix<float> GetH()
    {
        return m_H;
    }

    void setH(vnl_vector<float> & X)
    {
        m_X = X;
    }
    vnl_vector<float> GetX()
    {
        return m_X;
    }


    void SetOffset(std::vector<unsigned int> offset)
    {
        m_offset = offset;
    }
    std::vector< unsigned int > GetOffset()
    {
        return m_offset;
    }


private:

    std::vector< itkPointer > m_lRImages;

    std::vector< itkPointer > m_simulatedOutputImages;

    vnl_sparse_matrix<float>  m_H;

    vnl_vector<float> m_X;

    std::vector<unsigned int> m_offset;


};

}// end namespace

#endif
