/*==========================================================================
  
  Â© UniversitÃ© de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 28/08/2013
  Author(s): Frederic Champ (champ(at)unistra.fr)
  
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


#include "btkDwiReconstructionFilter.h"
namespace btk
{

DwiReconstructionFilter::DwiReconstructionFilter()
{
    m_Dataset = DatasetType::New();
    m_VolumeRegistration = true;
    m_SliceRegistration  = false;
    m_NumberOfIterations = 100;

}

//----------------------------------------------------------------------------------------

DwiReconstructionFilter::~DwiReconstructionFilter()
{
    // ----
}

//----------------------------------------------------------------------------------------

void DwiReconstructionFilter::PrintSelf(std::ostream &os, itk::Indent indent) const
{
    //Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------------------

void DwiReconstructionFilter::Initialize()
{
    ////////////////////////////////////////////////////////////////////////////
    // Used in case you whant to use the mean image as reference for slice by slice
    ////////////////////////////////////////////////////////////////////////////
    //    m_Region4D = m_InputSequence -> GetLargestPossibleRegion();
    //    m_Size4D = m_InputSequence -> GetLargestPossibleRegion().GetSize();
    //    m_GradientImages.resize(m_Size4D[3]-1);

    //    ////////////////////////////////////////////////////////////////////////////
    //    // Extract gradient images
    //    ////////////////////////////////////////////////////////////////////////////
    //    for(unsigned int i=1; i< m_Size4D[3]; i++)
    //    {
    //        SequenceType::RegionType  ImageRegion = m_InputSequence -> GetLargestPossibleRegion();
    //        ImageRegion.SetSize(3,0);
    //        ImageRegion.SetIndex(3,i);

    //        ExtractorPointer extractor = ExtractorType::New();
    //        extractor -> SetInput(m_InputSequence);
    //        extractor -> SetExtractionRegion(ImageRegion);
    //        extractor -> SetDirectionCollapseToSubmatrix( );
    //        extractor -> Update();

    //        m_GradientImages[i-1] = extractor -> GetOutput();
    //    }

}

void DwiReconstructionFilter::Update()
{

    SequenceType::RegionType  region4D = m_InputSequence -> GetLargestPossibleRegion();
    SequenceType::SizeType size4D = region4D.GetSize();
    std::vector<TransformPointer> initialTransforms(size4D[3]-1, TransformType::New());

    ////////////////////////////////////////////////////////////////////////////
    // Check reference image
    ////////////////////////////////////////////////////////////////////////////
    //
    // If reference image hasn't been set, set the B0 as reference
    //
    if(m_ReferenceImage.IsNull())
    {
        SequenceType::RegionType  B0region = region4D;
        B0region.SetSize(3,0);
        B0region.SetIndex(3,0);

        ExtractorPointer B0extractor = ExtractorType::New();
        B0extractor -> SetInput(m_InputSequence);
        B0extractor -> SetExtractionRegion(B0region);
        B0extractor -> SetDirectionCollapseToSubmatrix( );
        B0extractor -> Update();

        m_ReferenceImage = B0extractor -> GetOutput();

    }

    ////////////////////////////////////////////////////////////////////////////
    // Volume registration
    ////////////////////////////////////////////////////////////////////////////
    //
    // If m_VolumeRegistration = true, start with a volume to volume registration
    // in order to avoid mis overlapping between the reference and Dwi volumes.
    //
    JoinerPointer JoinerR = JoinerType::New();

    if(m_VolumeRegistration)
    {
        btkCoutMacro("Start Registration:");

        JoinerR -> SetOrigin( 0.0);
        JoinerR -> SetSpacing(m_InputSequence->GetSpacing()[3] );
        //Joiner -> SetInput(0,   m_ResampledImages[0]);

        for(unsigned int i=0; i< size4D[3]-1; i++)
        {

            ////////////////////////////////////////////////////////////////////////////
            // Extract the Nth Dwi
            ////////////////////////////////////////////////////////////////////////////
            SequenceType::RegionType  DwiRegion =region4D; // m_InputSequence -> GetLargestPossibleRegion();
            DwiRegion.SetSize(3,0);
            DwiRegion.SetIndex(3,i+1);

            ExtractorPointer DwiExtractor = ExtractorType::New();
            DwiExtractor -> SetInput(m_InputSequence);
            DwiExtractor -> SetExtractionRegion(DwiRegion);
            DwiExtractor -> SetDirectionCollapseToSubmatrix();
            DwiExtractor -> Update();

            ImagePointer DiffusionImage = DwiExtractor->GetOutput();

            ImageType::IndexType centerIndex;
            centerIndex[0] = DwiRegion.GetIndex(0) + DwiRegion.GetSize(0) / 2.0;
            centerIndex[1] = DwiRegion.GetIndex(1) + DwiRegion.GetSize(1) / 2.0;
            centerIndex[2] = DwiRegion.GetIndex(2) + DwiRegion.GetSize(2) / 2.0;

            ImageType::PointType  centerPoint;
            DiffusionImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);

            initialTransforms[i]=TransformType::New();
            initialTransforms[i] -> SetIdentity();
            initialTransforms[i] -> SetCenter(centerPoint);


            ////////////////////////////////////////////////////////////////////////////
            // Initialize loop through the image resolution
            ////////////////////////////////////////////////////////////////////////////
            unsigned int startResolution = 8; // Set the start resolution (mm per voxel's sides)
            unsigned int lastResolution = 2;  // Set the last resolution
            unsigned int numberOfStep = ((startResolution -lastResolution) /2) +1;

            unsigned int resolution = startResolution;
            for(unsigned int j=0;j<numberOfStep; j++)
            {

                std::cout<<" \r -> Volume "<<i<<" - Resolution "<<resolution<<" ... "<<std::flush;

                ////////////////////////////////////////////////////////////////////////////
                //Resample gradient images with the given resolution
                ////////////////////////////////////////////////////////////////////////////
                ResolutionFilterType::Pointer ResolutionSampler = ResolutionFilterType::New();
                ResolutionSampler -> SetInputImage(DiffusionImage.GetPointer());
                ResolutionSampler -> SetResolution(resolution);
                ResolutionSampler -> Update();

                ImagePointer resampledDwi = ResolutionSampler->GetOutput();

                ////////////////////////////////////////////////////////////////////////////
                // Resample the reference image on the same grid of resampled gradients
                ////////////////////////////////////////////////////////////////////////////
                TransformPointer id = TransformType::New();
                id -> SetIdentity();

                ResamplerType::Pointer resampler = ResamplerType::New();
                resampler->SetTransform(id);
                resampler->SetInput(m_ReferenceImage);
                resampler->SetUseReferenceImage(true);
                resampler->SetReferenceImage(resampledDwi);
                resampler->SetDefaultPixelValue(0);
                resampler->Update();

                ImagePointer resampledRef  = resampler -> GetOutput();
                ImageType::RegionType region  = resampledRef -> GetLargestPossibleRegion();

                ////////////////////////////////////////////////////////////////////////////
                //Initialize optimizer
                ////////////////////////////////////////////////////////////////////////////
                PowellOptimizerPointer optimizer = PowellOptimizerType::New();
                optimizer->SetMaximize(false);
                optimizer->SetStepLength(0.1);
                optimizer->SetStepTolerance( 1e-4);
                optimizer->SetValueTolerance(1e-4 );
                optimizer->SetMaximumIteration( m_NumberOfIterations );


                ////////////////////////////////////////////////////////////////////////////
                // Initialize metric
                ////////////////////////////////////////////////////////////////////////////
                MattesMetricPointer metric = MattesMetricType::New();
                metric -> SetNumberOfHistogramBins( 64/resolution);
                metric -> SetNumberOfSpatialSamples(0.5*region.GetNumberOfPixels());

                ////////////////////////////////////////////////////////////////////////////
                // Initialize interpolator
                ////////////////////////////////////////////////////////////////////////////
                InterPolatorPointer interpolator = InterpolatorType::New();

                ////////////////////////////////////////////////////////////////////////////
                // Connect components
                ////////////////////////////////////////////////////////////////////////////
                RegistrationPointer registration = RegistrationType::New();
                registration -> SetInitialTransformParameters(initialTransforms[i] -> GetParameters() );
                registration -> SetTransform(   initialTransforms[i] );
                registration -> SetMetric(       metric );
                registration -> SetOptimizer(    optimizer );
                registration -> SetInterpolator( interpolator );
                registration -> SetFixedImage(   resampledRef );
                registration -> SetMovingImage(  resampledDwi);
                registration -> SetFixedImageRegion(resampledDwi -> GetLargestPossibleRegion());
                registration -> Initialize();
                registration -> Update();

                initialTransforms[i] -> SetParameters(registration -> GetLastTransformParameters());

                resolution-=2;
                if(resolution <= 0)
                {
                    resolution =1;
                }
            }

            ResamplerType::Pointer resampler = ResamplerType::New();
            resampler->SetTransform(initialTransforms[i]);
            resampler->SetInput(DiffusionImage);
            resampler->SetUseReferenceImage(true);
            resampler->SetReferenceImage(DiffusionImage);
            resampler->SetDefaultPixelValue(0);
            resampler->Update();
            JoinerR -> SetInput(i,resampler->GetOutput());

        } // end for loop


        JoinerR -> Update();
        std::ostringstream osR;
        osR << "./Registred_sequence.nii.gz";
        btk::ImageHelper<  SequenceType>::WriteImage(JoinerR->GetOutput(), osR.str());

        btkCoutMacro(" Done.");



    } // end if fetal



    //////////////////////////////////////////////////////////////////////////
    //Initializations of the data structure with previous transforms
    //////////////////////////////////////////////////////////////////////////
    m_Dataset->SetInputSequence(m_InputSequence.GetPointer());
    m_Dataset -> SetInitialTransforms(initialTransforms);
    m_Dataset -> Initialize();
    m_Dataset -> RemoveOutliers();


    ////////////////////////////////////////////////////////////////////////////
    // Slice registration
    ////////////////////////////////////////////////////////////////////////////
    //
    // If m_SliceRegistration = true, start with a slice by slice registration
    //
    if(m_SliceRegistration)
    {

        //////////////////////////////////////////////////////////////////////////
        //Initializations of the data structure with previous transforms
        //////////////////////////////////////////////////////////////////////////
        //          this->Initialize()
        //        std::vector< itk::Image<double,3 >::Pointer >  RegistredImages(m_Size4D[3]-1);
        //        for(unsigned int r =0; r <RegistredImages.size(); r++)
        //        {
        //            typedef itk::CastImageFilter< ImageType,  itk::Image<double,3 > > CastFilterType;
        //            CastFilterType::Pointer castFilter = CastFilterType::New();
        //            castFilter->SetInput(m_GradientImages[r]);
        //            castFilter -> Update();
        //            RegistredImages[r]=castFilter -> GetOutput();
        //        }

        //        std::vector< double > weights(RegistredImages.size(),1.0f/(RegistredImages.size()));

        //        WeightedSumPointer sumFilter = WeightedSumType::New();
        //        sumFilter->SetInputs(RegistredImages);
        //        sumFilter->SetWeights(weights);
        //        sumFilter->Update();


        //        typedef itk::CastImageFilter< itk::Image<double,3 >,  ImageType> CastFilterType;
        //        CastFilterType::Pointer castFilter = CastFilterType::New();
        //        castFilter->SetInput(sumFilter->GetOutput());
        //        castFilter -> Update();


        TransformPointer id = TransformType::New();
        id -> SetIdentity();

        ResamplerType::Pointer resampler = ResamplerType::New();
        resampler->SetTransform(id);
        resampler->SetInput(m_ReferenceImage);
        resampler->SetUseReferenceImage(true);
        resampler->SetReferenceImage(m_Dataset->GetB0());
        resampler->SetDefaultPixelValue(0);
        resampler->Update();


        ImagePointer  meanImage= resampler -> GetOutput();
        btk::ImageHelper<  ImageType >::WriteImage(meanImage, "./mean_image.nii.gz");
        unsigned int loop =0;

        ImagePointer resampledRef = meanImage; //resampler->GetOutput();
        for(DatasetType::DataVectorType::iterator it = m_Dataset->begin(); it != m_Dataset->end(); it++)
        {
            DiffusionSlice::Pointer slice = (*it);

            ImageType::RegionType sliceRegion = slice -> GetLargestPossibleRegion();
            ImageType::RegionType refRegion = resampledRef-> GetLargestPossibleRegion();

            ImageType::IndexType centerIndex;
            centerIndex[0] = refRegion.GetIndex(0) + refRegion.GetSize(0) / 2.0;
            centerIndex[1] = refRegion.GetIndex(1) + refRegion.GetSize(1) / 2.0;
            centerIndex[2] = refRegion.GetIndex(2) + refRegion.GetSize(2) / 2.0;

            ImageType::PointType  centerPoint;
            slice->TransformIndexToPhysicalPoint(centerIndex, centerPoint);

            TransformPointer initialTransform=TransformType::New();
            initialTransform -> SetIdentity();
            initialTransform -> SetCenter(centerPoint);

            ///////////////////////////////////////////////////////////////////////////
            //Initialize optimizer
            ////////////////////////////////////////////////////////////////////////////
            PowellOptimizerPointer optimizer = PowellOptimizerType::New();
            optimizer->SetMaximize(false);
            optimizer->SetStepLength(0.01);
            optimizer->SetStepTolerance( 1e-4);
            optimizer->SetValueTolerance(1e-4 );
            optimizer->SetMaximumIteration( 100 );

            //        GradOptimizerPointer optimizer = GradOptimizerType::New();

            //        optimizer->MinimizeOn();
            //        optimizer->SetMaximumStepLength( 0.01 );
            //        optimizer->SetMinimumStepLength( 0.0001 );
            //        optimizer->SetNumberOfIterations( 100 );
            //        optimizer->SetRelaxationFactor( 0.8 );

            //        GradOptimizerType::ScalesType optimizerScales(12 );

            //        optimizerScales[0] =  1.0;
            //        optimizerScales[1] =  1.0;
            //        optimizerScales[2] =  1.0;
            //        optimizerScales[3] =  1.0;
            //        optimizerScales[4] =  1.0;
            //        optimizerScales[5] =  1.0;
            //        optimizerScales[6] =  1.0;
            //        optimizerScales[7] =  1.0;
            //        optimizerScales[8] =  1.0;
            //        optimizerScales[9] =  1.0 / 1000.0;
            //        optimizerScales[10] =  1.0 / 1000.0;
            //        optimizerScales[11] =  1.0 / 1000.0;

            //        optimizer->SetScales( optimizerScales );

            ////////////////////////////////////////////////////////////////////////////
            // Initialize metric
            ////////////////////////////////////////////////////////////////////////////

            MattesMetricPointer metric = MattesMetricType::New();
            metric -> SetNumberOfHistogramBins( 32);
            metric->UseAllPixelsOn();
            //        metric -> SetNumberOfSpatialSamples(0.8*sliceRegion.GetNumberOfPixels());
            //        NormalizedMetricPointer metric = NormalizedMetricType::New();


            ////////////////////////////////////////////////////////////////////////////
            // Initialize interpolator
            ////////////////////////////////////////////////////////////////////////////
            InterPolatorPointer interpolator = InterpolatorType::New();

            ////////////////////////////////////////////////////////////////////////////
            // Connect components
            ////////////////////////////////////////////////////////////////////////////
            RegistrationPointer registration = RegistrationType::New();
            registration -> SetInitialTransformParameters(initialTransform -> GetParameters() );
            registration -> SetTransform(   initialTransform );
            registration -> SetMetric(       metric );
            registration -> SetOptimizer(    optimizer );
            registration -> SetInterpolator( interpolator );
            registration -> SetFixedImage(  slice  );
            registration -> SetMovingImage( resampledRef );
            registration -> SetFixedImageRegion(sliceRegion);
            registration -> Initialize();

            try{
                registration -> Update();
                initialTransform -> SetParameters(registration -> GetLastTransformParameters());
                initialTransform->GetInverse(initialTransform);
                slice->Transform(initialTransform);

            }
            catch (itk::ExceptionObject)
            {
                // --- ?
            }



            loop++;
            btkCoutMacro("  -> Number of slices registred = "<<loop<<" "<<slice->GetOrigin()<<" "<<registration->GetLastTransformParameters()<<" (total: "<<m_Dataset->GetSize()<<" )");
            //        std::cout<<"  -> Number of slices registred = "<<loop<<" "<<slice->GetOrigin<<" "<<initialTransform->GetParameters()<<" (total: "<<m_Dataset->GetSize()<<" )"<<std::endl; //<<std::flush;

            //        }
        }
        btkCoutMacro(" Done.");
    }




}


} // end namespace btk
