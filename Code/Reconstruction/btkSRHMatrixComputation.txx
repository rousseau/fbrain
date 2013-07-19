/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 27/05/2013
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

#include "btkSRHMatrixComputation.hxx"

namespace btk
{
//-------------------------------------------------------------------------------------------------

template < class TImage >
SRHMatrixComputation< TImage >::SRHMatrixComputation():m_PSFType(0)
{
}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void SRHMatrixComputation< TImage >::Initialize()
{
    m_IsHComputed = false;

    m_Regions.resize(m_Images.size());
    m_PSFs.resize(m_Images.size());
    for(unsigned int i = 0; i< m_Images.size(); i++)
    {
        m_Regions[i] = m_Images[i]->GetLargestPossibleRegion();
        if(m_PSFType == 0)
        {
            m_PSFs[i] = btk::GaussianPSF::New();
        }
        else if(m_PSFType == 1)
        {
            m_PSFs[i] = btk::BoxCarPSF::New();
        }
        else if(m_PSFType == 2)
        {
            m_PSFs[i] = btk::GaussianPSF::New();
        }
        else
        {
            m_PSFs[i] = btk::GaussianPSF::New();
        }

    }
}
//-------------------------------------------------------------------------------------------------
template< class TImage >
void SRHMatrixComputation< TImage >::Update()
{


    this->Initialize();

    m_OutputImageRegion = m_ReferenceImage -> GetLargestPossibleRegion();
    IndexType start_hr  = m_OutputImageRegion.GetIndex();
    SizeType  size_hr   = m_OutputImageRegion.GetSize();

    //m_XSize : size of the SR image (used in other functions)
    m_XSize.width  = size_hr[0];
    m_XSize.height = size_hr[1];
    m_XSize.depth  = size_hr[2];

    IndexType end_hr;
    end_hr[0] = start_hr[0] + size_hr[0] - 1 ;
    end_hr[1] = start_hr[1] + size_hr[1] - 1 ;
    end_hr[2] = start_hr[2] + size_hr[2] - 1 ;

    // Differential continuous indexes to perform the neighborhood iteration
    SpacingType spacing_lr = m_Images[0] -> GetSpacing();
    SpacingType spacing_hr = m_ReferenceImage -> GetSpacing();

    //spacing_lr[2] is assumed to be the lowest resolution
    //compute the index of the PSF in the LR image resolution
    std::vector<ContinuousIndexType> deltaIndexes;
    int npoints =  spacing_lr[2] / (2.0 * spacing_hr[2]) ;
    ContinuousIndexType delta;
    delta[0] = 0.0; delta[1] = 0.0;


    for (int i = -npoints ; i <= npoints; i++ )
    {
      delta[2] = i * 0.5 / (double) npoints;
      deltaIndexes.push_back(delta);
    }

    // Set size of matrices
    unsigned int ncols = m_OutputImageRegion.GetNumberOfPixels();

    m_X.set_size( ncols );



    unsigned int nrows = 0;
    for(unsigned int im = 0; im < m_Images.size(); im++)
    {
        nrows += m_Regions[im].GetNumberOfPixels();
    }

    m_H.set_size(nrows, ncols);
    m_Y.set_size(nrows);
    m_Y.fill(0.0);

    std::cout<<"size X : "<<ncols<<std::endl;
    std::cout<<"size Y : "<<nrows<<std::endl;

    unsigned int im;
    #pragma omp parallel for private(im) schedule(dynamic)

    for(im = 0; im < m_Images.size(); im++)
    {

      // Interpolator for HR image
      InterpolatorPointer interpolator = InterpolatorType::New();
      interpolator -> SetInputImage( m_ReferenceImage );

      SpacingType inputSpacing = m_Images[im] -> GetSpacing();

      // PSF definition
//      typename FunctionType::Pointer function = FunctionType::New();
//      function -> SetPSF( 0 );
//      function -> SetDirection( m_Images[im] -> GetDirection() );
//      function -> SetSpacing( m_Images[im] -> GetSpacing() );



      m_PSFs[im]->SetDirection(m_Images[im] -> GetDirection());
      //m_PSFs[im]->SetSpacing(m_Images[im] -> GetSpacing());// spacing of psf should be at least equal to reference image ;
      m_PSFs[im]->SetSpacing(m_ReferenceImage->GetSpacing());

      typename ImageType::SizeType size;

      size[0] = (int)ceil(m_Images[im]->GetSpacing()[0] / m_ReferenceImage->GetSpacing()[0]) +3;
      size[1] = (int)ceil(m_Images[im]->GetSpacing()[1] / m_ReferenceImage->GetSpacing()[1]) +3;
      size[2] = (int)ceil(m_Images[im]->GetSpacing()[2] / m_ReferenceImage->GetSpacing()[2]) +3;

      m_PSFs[im]->ConstructImage(size);

      //Define the ROI of the current LR image
      IndexType inputIndex = m_Regions[im].GetIndex();
      SizeType  inputSize  = m_Regions[im].GetSize();

      //Define all indexes needed for iteration over the slices
      IndexType lrIndex;              //index of a point in the LR image im
      IndexType lrDiffIndex;          //index of this point in the current ROI of the LR image im
      unsigned int lrLinearIndex;     //index lineaire de ce point dans le vecteur

      IndexType hrIndex;
      IndexType hrDiffIndex;
      ContinuousIndexType hrContIndex;
      unsigned int hrLinearIndex;

      ContinuousIndexType nbIndex;

      PointType lrPoint;            //PSF center in world coordinates (PointType = worldcoordinate for ITK)
      PointType nbPoint;            //PSF point in world coordinate
      PointType transformedPoint;   //after applying current transform (registration)

      unsigned int offset = 0;
      for(unsigned int im2 = 0; im2 < im; im2++)
      {
          offset += m_Regions[im2].GetNumberOfPixels();
      }
      unsigned int i=inputIndex[2];

      // Iteration over slices
      for ( i=inputIndex[2]; i < inputIndex[2] + inputSize[2]; i++ )
      {

        RegionType wholeSliceRegion;
        wholeSliceRegion = m_Regions[im];

        IndexType  wholeSliceRegionIndex = wholeSliceRegion.GetIndex();
        SizeType   wholeSliceRegionSize  = wholeSliceRegion.GetSize();

        wholeSliceRegionIndex[2]= i;
        wholeSliceRegionSize[2] = 1;

        wholeSliceRegion.SetIndex(wholeSliceRegionIndex);
        wholeSliceRegion.SetSize(wholeSliceRegionSize);

        ConstIteratorType fixedIt( m_Images[im], wholeSliceRegion);

        double lrValue;
        double hrValue;

        for(fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt)
        {
          //Current index in the LR image
          lrIndex = fixedIt.GetIndex();

          //World coordinates of lrIndex using the image header
          m_Images[im] -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );

          if ( m_Masks.size() > 0)
            if ( ! m_Masks[im] -> IsInside(lrPoint) )
              continue;

          //Compute the coordinates in the SR using the estimated registration
          transformedPoint = m_Transforms[im]->TransformPoint( lrPoint );

          //check if this point is in the SR image (m_ReferenceImage)
          if ( ! interpolator -> IsInsideBuffer( transformedPoint ) )
            continue;

          //From the LR image coordinates to the LR ROI coordinates
          lrDiffIndex[0] = lrIndex[0] - inputIndex[0];
          lrDiffIndex[1] = lrIndex[1] - inputIndex[1];
          lrDiffIndex[2] = lrIndex[2] - inputIndex[2];

          //Compute the corresponding linear index of lrDiffIndex
          lrLinearIndex = lrDiffIndex[0] + lrDiffIndex[1]*inputSize[0] +
                          lrDiffIndex[2]*inputSize[0]*inputSize[1];

          //Get the intensity value in the LR image
          m_Y[lrLinearIndex + offset] = fixedIt.Get();

  //        if(m_UseOutliers)
  //        {
  //            if(m_Outliers[im][i])
  //            {
  //                if(i>0 && i<inputIndex[2] + inputSize[2])
  //                {
  //                    double p = 0.0;
  //                    lrIndex[2] = i-1;
  //                    p+=m_Images[im]->GetPixel(lrIndex);
  //                    lrIndex[2] = i+1;
  //                    p+=m_Images[im]->GetPixel(lrIndex);
  //                    lrIndex[2]= i;
  //                    p = p/2.0;
  //                    Y[lrLinearIndex + offset] = p;
  //                }
  //                else if(i==0)
  //                {
  //                    double p = 0.0;
  //                    lrIndex[2] = i+1;
  //                    p+=m_Images[im]->GetPixel(lrIndex);
  //                    lrIndex[2]= i;
  //                    Y[lrLinearIndex + offset] = p;
  //                }
  //                else if(i ==inputIndex[2] + inputSize[2] -1 )
  //                {
  //                    double p = 0.0;
  //                    lrIndex[2] = i-1;
  //                    p+=m_Images[im]->GetPixel(lrIndex);
  //                    lrIndex[2]= i;
  //                    Y[lrLinearIndex + offset] = p;
                 // }
  //            }
  //        }

          //Set the center point of the PSF
          //function -> SetCenter( lrPoint );
          m_PSFs[im]->SetCenter(lrPoint);

          typename PsfImageType::Pointer PSF = PsfImageType::New();
          PSF = m_PSFs[im]->GetPsfImage();

          //btk::ImageHelper< PsfImageType >::WriteImage(PSF, "GaussianPSF.nii.gz");
          ImageRegionConstIteratorWithIndex< itk::Image< float, 3> > itPsf(PSF, PSF->GetLargestPossibleRegion());



          //Loop over points of the PSF
          //for(unsigned int k=0; k<deltaIndexes.size(); k++)
          for(itPsf.GoToBegin(); !itPsf.IsAtEnd(); ++itPsf)
          {
//            //Coordinates in the LR image
//            nbIndex[0] = deltaIndexes[k][0] + lrIndex[0];//0 + ...
//            nbIndex[1] = deltaIndexes[k][1] + lrIndex[1];//0 + ...
//            nbIndex[2] = deltaIndexes[k][2] + lrIndex[2];//-0.5, 0, +0.5

              nbIndex = itPsf.GetIndex();



            //World coordinates using LR image header
            //m_Images[im] -> TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );
            PSF->TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );
            //Compute the PSF value at this point
           // lrValue = function -> Evaluate(nbPoint);
            //lrValue = m_PSFs[im]->Evaluate(nbPoint);
            lrValue = itPsf.Get();

            //std::cout<<lrValue<<std::endl;

            if ( lrValue > 0)
            {
              //Compute the world coordinate of this point in the SR image
              transformedPoint = m_Transforms[im]->TransformPoint( nbPoint );

              //Set this coordinate in continuous index in SR image space
              m_ReferenceImage -> TransformPhysicalPointToContinuousIndex(
                                  transformedPoint, hrContIndex );

              bool isInsideHR = true;

              // FIXME This checking should be done for all points first, and
              // discard the point if al least one point is out of the reference
              // image

//              if ( (hrContIndex[0] < start_hr[0]) || (hrContIndex[0] > end_hr[0]) ||
//                   (hrContIndex[1] < start_hr[1]) || (hrContIndex[1] > end_hr[1]) ||
//                   (hrContIndex[2] < start_hr[2]) || (hrContIndex[2] > end_hr[2]) )
//                 isInsideHR = false;




//              if ( isInsideHR )
//              {
              if(m_OutputImageRegion.IsInside(hrContIndex))
              {

                //Compute the corresponding value in the SR image -> useless
                //Allows to compute the set of contributing neighbors
                hrValue = interpolator -> Evaluate( transformedPoint );

                //Loop over points affected using the interpolation
                for(unsigned int n=0; n<interpolator -> GetContributingNeighbors();
                    n++)
                {
                  //Index in the SR image
                  hrIndex = interpolator -> GetIndex(n);



                  //Index in the ROI of the SR index
                  hrDiffIndex[0] = hrIndex[0] - start_hr[0];
                  hrDiffIndex[1] = hrIndex[1] - start_hr[1];
                  hrDiffIndex[2] = hrIndex[2] - start_hr[2];

                  //Compute the corresponding linear index
                  hrLinearIndex = hrDiffIndex[0] + hrDiffIndex[1]*size_hr[0] +
                                  hrDiffIndex[2]*size_hr[0]*size_hr[1];

                  //Add the correct value in H !
                  m_H(lrLinearIndex + offset, hrLinearIndex) += interpolator -> GetOverlap(n)* lrValue;


                }

              }

            }

          }

        }

      }

    }

    // Normalize H
    std::cout<<"H rows : "<<m_H.rows()<<std::endl;
    std::cout<<"H cols : "<<m_H.cols()<<std::endl;
    for (unsigned int i = 0; i < m_H.rows(); i++)
    {
      double sum = m_H.sum_row(i);

      //VnlSparseMatrixType::row & r = m_H.get_row(i);
      //VnlSparseMatrixType::row::iterator col_iter = r.begin();

      //vnl_sparse_matrix< float >::row & r = m_H.get_row(i);
      //vnl_sparse_matrix< float >::row::iterator col_iter = r.begin();

      vnl_sparse_matrix< PrecisionType >::row & r = m_H.get_row(i);
      vnl_sparse_matrix< PrecisionType >::row::iterator col_iter = r.begin();

      for ( ;col_iter != r.end(); ++col_iter)
      {
        (*col_iter).second = (*col_iter).second / sum;
          //std::cout<<"H: "<<(*col_iter).second <<std::endl;

      }


    }
    //std::cout<<"size : "<<s<<std::endl;
    //std::cout<<" size Y : "<<m_Y.size()<<std::endl;

    // Precalcule Ht*Y. Note that this is calculated as Y*H since
    // Ht*Y = (Yt*H)t and for Vnl (Yt*H)t = (Yt*H) = Y*H because
    // the vnl_vector doesn't have a 2nd dimension. This allows us
    // to save a lot of memory because we don't need to store Ht.

    m_H.pre_mult(m_Y,m_HtY);




    m_IsHComputed = true;
}

//-------------------------------------------------------------------------------------------------
template< class TImage >
void
SRHMatrixComputation< TImage >::ComputeSimulatedLRImages()
{
    if(!m_IsHComputed )
    {
        this->Update();
    }
    // H should previoulsy be computed
    m_H.mult(m_X,m_SimY);

    unsigned int offset = 0;

    for (unsigned int im = 0; im < m_Images.size(); im++)
    {
      IndexType absIndex;

      IndexType start = m_Regions[im].GetIndex();
      SizeType  size  = m_Regions[im].GetSize();
      unsigned int nvoxels = m_Regions[im].GetNumberOfPixels();
      IndexType diffIndex;

      for( unsigned int i=0; i<nvoxels; i++)
      {

        diffIndex[2] = i / (size[0]*size[1]);

        diffIndex[1] = i - diffIndex[2]*size[0]*size[1];
        diffIndex[1] = diffIndex[1] / size[0];

        diffIndex[0] = i - diffIndex[2]*size[0]*size[1] - diffIndex[1]*size[0];

        absIndex[0] = diffIndex[0] + start[0];
        absIndex[1] = diffIndex[1] + start[1];
        absIndex[2] = diffIndex[2] + start[2];

        m_SimulatedImages[im]->SetPixel(absIndex,m_SimY[i + offset]);

      }

      offset = offset + nvoxels;

    }

    m_X.clear();
}
//-------------------------------------------------------------------------------------------------
}//namespace
