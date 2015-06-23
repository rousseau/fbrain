/*
 Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 
 15 december 2014
 rousseau@unistra.fr
 
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
 */

/* Standard includes */
#include <tclap/CmdLine.h>
#include "vector"
#include "sstream"
//#include "numeric"      // std::accumulate

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <limits>

/* Itk includes */
#include "itkImage.h"
#include "itkVector.h"

#include "itkMatrixOffsetTransformBase.h"
//#include "itkEuler3DTransform.h"
#include "btkIOTransformHelper.h"
//#include "itkTransformFileReader.h"
//#include "itkResampleImageFilter.h"

//#include "vnl/vnl_sparse_matrix.h"
#include <omp.h>

/*Btk includes*/
#include "btkImageHelper.h"
#include "btkPandoraBoxImageFilters.h"
#include "btkPandoraBoxRegistrationFilters.h"
#include "btkPandoraBoxCostFunction.h"

inline bool less_than_second( const std::pair< vnl_vector< double >, double > & p1, const std::pair< vnl_vector< double >, double > & p2 ){
   return p1.second < p2.second;
}

int main(int argc, char** argv)
{
  try {
    
    TCLAP::CmdLine cmd("3D image registration", ' ', "0.1", true);
    
    TCLAP::ValueArg<std::string> refImageArg        ("r","reference","reference image",true,"","string",cmd);
    TCLAP::ValueArg<std::string> movImageArg        ("m","moving","moving image",true,"","string",cmd);
    TCLAP::ValueArg<std::string> refMaskArg         ("","refmask","reference mask",false,"","string",cmd);
    TCLAP::ValueArg<std::string> movMaskArg         ("","movmask","moving mask",false,"","string",cmd);
    TCLAP::ValueArg<std::string> registeredImageArg ("o","output","registered image file",false,"","string",cmd);
    TCLAP::ValueArg<std::string> transformArg       ("t","transf","estimated transform",false,"","string",cmd);
    TCLAP::ValueArg< int >       dofArg             ("","dof","Number of degrees of freedom (6,7,9,12)",false,6,"int",cmd);
    TCLAP::MultiArg< float>      scaleArg           ("","scale","scales for registration (default: 8,4,2)",false,"float",cmd);
    TCLAP::ValueArg< int >       similarityArg      ("s","sim","similarity measure (MSE:0 (default), MI:1, NMI:2, openMP MSE: 3, openMP MI: 4, openMP NMI: 5)",false,0,"int",cmd);
    TCLAP::ValueArg< int >       binArg             ("b","bin","number of bin of the joint histogram (default: 64)",false,64,"int",cmd);
    TCLAP::ValueArg< int >       bestCandidatesArg  ("","best","number of best candidates kept (default: 2)",false,2,"int",cmd);
    TCLAP::ValueArg< int >       numberOfPerturbationsArg  ("","perturb","number of perturbations (default: 2)",false,2,"int",cmd);
    TCLAP::ValueArg< float >     rotationXRangeArg  ("","rx","range of rotation (x axis) (default: 10)",false,10,"float",cmd);
    TCLAP::ValueArg< float >     rotationYRangeArg  ("","ry","range of rotation (y axis) (default: 10)",false,10,"float",cmd);
    TCLAP::ValueArg< float >     rotationZRangeArg  ("","rz","range of rotation (z axis) (default: 10)",false,10,"float",cmd);
    TCLAP::ValueArg< int >       numberOfStartingEstimateArg  ("","start","number of starting estimates (per rotation axis) (default: 3)",false,3,"int",cmd);
    TCLAP::ValueArg< int >       orderInterpolationArg ("","order","order of the interpolation spline",false,1,"int",cmd);

    
    //TODO

    // Parse the args.
    cmd.parse( argc, argv );
    
    // Get the value parsed by each arg.

    std::string                reference_image_filename   = refImageArg.getValue();
    std::string                moving_image_filename      = movImageArg.getValue();
    std::string                reference_mask_filename    = refMaskArg.getValue();
    std::string                moving_mask_filename       = movMaskArg.getValue();

    std::string                registered_image_filename  = registeredImageArg.getValue();
    std::string                output_transform_filename  = transformArg.getValue();

    //DOF MAX -> need to adjust code for that **************************************************
    int       dof              = dofArg.getValue();
    int       similarity       = similarityArg.getValue();
    int       bin              = binArg.getValue();
    std::vector<float> scales  = scaleArg.getValue();
    int numberOfBestCandidates = bestCandidatesArg.getValue();
    int numberOfPerturbations  = numberOfPerturbationsArg.getValue();

    vnl_vector< double > rotationRange(3,0);
    rotationRange[0] = rotationXRangeArg.getValue();
    rotationRange[1] = rotationYRangeArg.getValue();
    rotationRange[2] = rotationZRangeArg.getValue();
    int samplingRate = numberOfStartingEstimateArg.getValue();
    int interpolationOrder = orderInterpolationArg.getValue();



    if(scales.size() < 1)
    {
        scales.push_back(8);
        scales.push_back(4);
        scales.push_back(2);
    }
    //Ordering scales if necessary
    std::sort(scales.begin(), scales.end());
    std::reverse(scales.begin(), scales.end());
    
    typedef itk::Image< float, 3 >                                   itkFloatImage;
    typedef itk::MatrixOffsetTransformBase<double,3,3>               itkTransformType;
    //itk::TransformFactory<itkTransformType>::RegisterTransform();
    typedef itk::Vector<double, 3>                              itkVector;

    
    //*******************************************************************************************************
    // READING IMAGE DATA
    //*******************************************************************************************************
    
    //Read images---------------------------------------------------------------------------------------------------
    itkFloatImage::Pointer referenceImage = btk::ImageHelper< itkFloatImage > ::ReadImage(reference_image_filename);
    itkFloatImage::Pointer movingImage    = btk::ImageHelper< itkFloatImage > ::ReadImage(moving_image_filename);
    
    std::cout<<"Reference image : \n";
    btk::PandoraBoxImageFilters::DisplayImageInfo(referenceImage);
    std::cout<<"Moving image : \n";
    btk::PandoraBoxImageFilters::DisplayImageInfo(movingImage);
    
    //Read masks----------------------------------------------------------------------------------------------------
    itkFloatImage::Pointer referenceMask = btk::ImageHelper< itkFloatImage >::ReadOrCreateImage(reference_mask_filename,referenceImage,1);
    itkFloatImage::Pointer movingMask    = btk::ImageHelper< itkFloatImage >::ReadOrCreateImage(moving_mask_filename,movingImage,1);
        
    //Computer center of the reference image
    //It will be used as center of the transform

    itk::Vector<double, 3>  center;

    itkFloatImage::SizeType    size    = referenceImage->GetLargestPossibleRegion().GetSize();
    //itkFloatImage::IndexType centerIndex;
    itk::ContinuousIndex<double,3> centerIndex;

    centerIndex[0] = (size[0]-1)/2.0;
    centerIndex[1] = (size[1]-1)/2.0;
    centerIndex[2] = (size[2]-1)/2.0;
    itkFloatImage::PointType refPoint; //physical point location of reference image

    //referenceImage->TransformIndexToPhysicalPoint(centerIndex,refPoint);
    referenceImage->TransformContinuousIndexToPhysicalPoint(centerIndex,refPoint);
    center[0] = refPoint[0];
    center[1] = refPoint[1];
    center[2] = refPoint[2];

    std::cout<<"Image center of the estimated transform:"<<centerIndex<<std::endl;
    std::cout<<"Physical center of the estimated transform:"<<center<<std::endl;

    //Initialisation of the registration parameters : by default, we assume that the headers provide relevant information
    //Otherwise, we should initialize the parameters using image moments, or center etc.
    vnl_vector< double > inputParam(dof,0);
    if(dof > 6)
        for(unsigned int i=6; i < dof; i++)
            inputParam(i) = 1.0;
    if(dof > 9)
        for(unsigned int i=9; i < dof; i++)
            inputParam(i) = 0.0;

    //define center of transformation
    vnl_vector< double > outputParam(inputParam.size(),0);


    itkFloatImage::SpacingType spacingReferenceImage = referenceImage->GetSpacing();

    //parameterRange should be filled with cautious (for scaling for instance). It defines the size of the simplex.
    vnl_vector< double > finalParameterRange(inputParam.size(),0.05);
    //Rotation
    finalParameterRange(0) = 5;
    finalParameterRange(1) = 5;
    finalParameterRange(2) = 5;
    //Translations (could be computed from voxel size)
    finalParameterRange(3) = 1 * spacingReferenceImage[0];
    finalParameterRange(4) = 1 * spacingReferenceImage[1];
    finalParameterRange(5) = 1 * spacingReferenceImage[2];

    //The final tolerance defines the accuracy we want for each registration parameter
    vnl_vector< double > finalTolerance(inputParam.size(),0.01);
    //Rotation
    finalTolerance(0) = 0.1;
    finalTolerance(1) = 0.1;
    finalTolerance(2) = 0.1;
    //Translations (could be computed from voxel size)
    finalTolerance(3) = 0.1 * spacingReferenceImage[0];
    finalTolerance(4) = 0.1 * spacingReferenceImage[1];
    finalTolerance(5) = 0.1 * spacingReferenceImage[2];

    //Scaling and skew : 0.01 by default


    //Core of the algorithm :
    //1- Choose the similarity measure then loop over :
    //2- Choose the resolution
    //3- Setup the set of possible parameters
    //4- Possibly rerun the optimization using perturbations

    btk::PandoraBoxCostFunction * myCostFunction;
    switch(similarity)
    {
        case 0:
            myCostFunction = new btk::PandoraBoxCostFunctionMSE();
            break;
        case 1:
            myCostFunction = new btk::PandoraBoxCostFunctionMI();
            break;
        case 2:
            myCostFunction = new btk::PandoraBoxCostFunctionNMI();
            break;
        case 3:
            myCostFunction = new btk::PandoraBoxCostFunctionMSEOpenMP();
            break;
        case 4:
            myCostFunction = new btk::PandoraBoxCostFunctionMIOpenMP();
            break;
        case 5:
            myCostFunction = new btk::PandoraBoxCostFunctionNMIOpenMP();
            break;
        default:
            myCostFunction = new btk::PandoraBoxCostFunctionMSE();
    }

    itkFloatImage::Pointer tmpReferenceImage;
    itkFloatImage::Pointer tmpMovingImage;
    itkFloatImage::Pointer tmpReferenceMask;
    itkFloatImage::Pointer tmpMovingMask;

    //FIRST SCALE *********************************************************************************
    //Fast approximate search of good sets of parameters
    float factor = scales[0];
    std::cout<<"Factor : "<<factor<<std::endl;

    //INITIALIZATION (same of all scales) ----------------------------------------------------------
    //Downscale with respect to the finest scaling factor of the original images

    btk::PandoraBoxImageFilters::DownscaleImage(referenceImage, tmpReferenceImage,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(movingImage, tmpMovingImage,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(referenceMask, tmpReferenceMask,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(movingMask, tmpMovingMask,factor);

    myCostFunction->SetReferenceImage(tmpReferenceImage);
    myCostFunction->SetMovingImage(tmpMovingImage);
    myCostFunction->SetMovingMask(tmpMovingMask);
    myCostFunction->SetReferenceMask(tmpReferenceMask);
    myCostFunction->SetCenter(center);

    if( (similarity==1) || (similarity==2) || (similarity==4) || (similarity==5) )
        myCostFunction->InitializeJointHistogram(bin,bin);

    vnl_vector< double > tmpParameterRange = finalParameterRange;
    vnl_vector< double > tmpTolerance      = finalTolerance;

    //update accordingly to image resolution
    for(unsigned int i=3; i < 6; i++)
    {
        tmpParameterRange(i) = finalParameterRange(i) * factor / 2;
    }
    for(unsigned int i=0; i < 6; i++)
    {
        tmpTolerance(i) = finalTolerance(i) * factor;
    }
    //END OF INITIALIZATION ------------------------------------------------------------------------

    //Generate parameters (using uniform sampling over rotation domain)
    std::vector< vnl_vector< double > > candidates;
    btk::PandoraBoxRegistrationFilters::GenerateUniformlyDistributedParameters(inputParam, candidates, rotationRange, samplingRate);

    std::vector< std::pair< vnl_vector< double >, double > > pairOfCandidates(candidates.size());

    std::cout<<"Run registration for each set of candidate parameters\n";
    for(unsigned int i=0; i < candidates.size(); i++)
    {
        btk::PandoraBoxRegistrationFilters::Register3DImages(*myCostFunction, candidates[i], outputParam, tmpParameterRange, tmpTolerance);
        //std::cout<<"Estimated parameters : "<<outputParam<<std::endl;
        pairOfCandidates[i].first = outputParam;
        pairOfCandidates[i].second = (*myCostFunction)(outputParam);
        //A cost function value equal to 0 is not possible (meaning that there is no more overlap between the two images.
        if(pairOfCandidates[i].second == 0)
            pairOfCandidates[i].second = std::numeric_limits<double>::max();
        //std::cout<<"Cost function : "<< pairOfCandidates[i].second <<std::endl;
    }

    std::cout<<"Sort results by cost function values\n";
    std::sort( pairOfCandidates.begin(), pairOfCandidates.end(), less_than_second );
    for(unsigned int i=0; i < candidates.size(); i++)
    {
        std::cout<<"param : "<<pairOfCandidates[i].first<<", cost : "<<pairOfCandidates[i].second<<std::endl;
    }
    std::vector< vnl_vector< double > > bestCandidates(numberOfBestCandidates);
    for(unsigned int i=0; i < numberOfBestCandidates; i++)
        bestCandidates[i] = pairOfCandidates[i].first;


    //SECOND SCALE *********************************************************************************
    //Keep the best candidates of the previous scale, add perturbations and optimize
    if(scales.size() > 1)
        factor = scales[1];
    else
        factor = scales[0];

    std::cout<<"Factor : "<<factor<<std::endl;

    //INITIALIZATION (same of all scales) ----------------------------------------------------------
    //Downscale with respect to the finest scaling factor of the original images
    btk::PandoraBoxImageFilters::DownscaleImage(referenceImage, tmpReferenceImage,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(movingImage, tmpMovingImage,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(referenceMask, tmpReferenceMask,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(movingMask, tmpMovingMask,factor);

    myCostFunction->SetReferenceImage(tmpReferenceImage);
    myCostFunction->SetMovingImage(tmpMovingImage);
    myCostFunction->SetMovingMask(tmpMovingMask);
    myCostFunction->SetReferenceMask(tmpReferenceMask);
    myCostFunction->SetCenter(center);

    if( (similarity==1) || (similarity==2) || (similarity==4) || (similarity==5) )
        myCostFunction->InitializeJointHistogram(bin,bin);

    tmpParameterRange = finalParameterRange;
    tmpTolerance      = finalTolerance;

    //update accordingly to image resolution
    for(unsigned int i=3; i < 6; i++)
    {
        tmpParameterRange(i) = finalParameterRange(i) * factor / 2;
    }
    for(unsigned int i=0; i < 6; i++)
    {
        tmpTolerance(i) = finalTolerance(i) * factor;
    }
    //END OF INITIALIZATION ------------------------------------------------------------------------

    std::vector< std::pair< vnl_vector< double >, double > > pairOfParamsAndCostFunctionValues(numberOfBestCandidates * numberOfPerturbations);

    //Generate random registration parameters
    srand (time(NULL));
    int index = 0;
    for(unsigned int i=0; i < numberOfBestCandidates; i++)
        for(unsigned int j=0; j < numberOfPerturbations; j++)
        {
            btk::PandoraBoxRegistrationFilters::GenerateRandomParameters(bestCandidates[i], pairOfParamsAndCostFunctionValues[index].first, tmpParameterRange);
            index++;
        }

    //Run registration for each set of parameters
    for(unsigned int i=0; i < numberOfBestCandidates * numberOfPerturbations; i++)
    {
        btk::PandoraBoxRegistrationFilters::Register3DImages(*myCostFunction, pairOfParamsAndCostFunctionValues[i].first, outputParam, tmpParameterRange, tmpTolerance);
        std::cout<<"Estimated parameters : "<<outputParam<<std::endl;
        pairOfParamsAndCostFunctionValues[i].first = outputParam;
        pairOfParamsAndCostFunctionValues[i].second = (*myCostFunction)(outputParam);
        //A cost function value equal to 0 is not possible (meaning that there is no more overlap between the two images.
        if(pairOfParamsAndCostFunctionValues[i].second == 0)
            pairOfParamsAndCostFunctionValues[i].second = std::numeric_limits<double>::max();
        std::cout<<"Cost function : "<< pairOfParamsAndCostFunctionValues[i].second <<std::endl;
    }

    //Sort results by cost function values
    std::sort( pairOfParamsAndCostFunctionValues.begin(), pairOfParamsAndCostFunctionValues.end(), less_than_second );
    for(unsigned int i=0; i < numberOfPerturbations; i++)
        std::cout<<"param : "<<pairOfParamsAndCostFunctionValues[i].first<<", cost : "<<pairOfParamsAndCostFunctionValues[i].second<<std::endl;

    //Assign best estimations
    vnl_vector< double > bestParamSoFar(inputParam.size(),0);

    bestParamSoFar = pairOfParamsAndCostFunctionValues[0].first;


    //THE OTHER SCALES ****************************************************************************************
    //Optimize only one set of parameters
    if(scales.size() > 2){
        for(unsigned int scaleIndex=2; scaleIndex < scales.size(); scaleIndex++)
        {
            factor = scales[scaleIndex];
            std::cout<<"Factor : "<<factor<<std::endl;

            //Downscale with respect to the finest scaling factor of the original images
            btk::PandoraBoxImageFilters::DownscaleImage(referenceImage, tmpReferenceImage,factor);
            btk::PandoraBoxImageFilters::DownscaleImage(movingImage, tmpMovingImage,factor);
            btk::PandoraBoxImageFilters::DownscaleImage(referenceMask, tmpReferenceMask,factor);
            btk::PandoraBoxImageFilters::DownscaleImage(movingMask, tmpMovingMask,factor);

            myCostFunction->SetReferenceImage(tmpReferenceImage);
            myCostFunction->SetMovingImage(tmpMovingImage);
            myCostFunction->SetMovingMask(tmpMovingMask);
            myCostFunction->SetReferenceMask(tmpReferenceMask);
            myCostFunction->SetCenter(center);

            if( (similarity==1) || (similarity==2) || (similarity==4) || (similarity==5) )
                myCostFunction->InitializeJointHistogram(bin,bin);

            tmpParameterRange = finalParameterRange;
            tmpTolerance      = finalTolerance;

            //update accordingly to image resolution
            for(unsigned int i=3; i < 6; i++)
            {
                tmpParameterRange(i) = finalParameterRange(i) * factor / 2;
            }
            for(unsigned int i=0; i < 6; i++)
            {
                tmpTolerance(i) = finalTolerance(i) * factor;
            }
            //END OF INITIALIZATION ------------------------------------------------------------------------

            btk::PandoraBoxRegistrationFilters::Register3DImages(*myCostFunction, bestParamSoFar, outputParam, tmpParameterRange, tmpTolerance);
            std::cout<<"Estimated parameters : "<<outputParam<<std::endl;
            bestParamSoFar = outputParam;

        }
    }

/*
    // FIRST LEVEL : FACTOR 8
    std::cout<<"Factor : "<<factor<<std::endl;
    int numberOfBestCandidates = 2;
    int numberOfPerturbations = 2;

    //Downscale with respect to the finest scaling factor of the original images
    btk::PandoraBoxImageFilters::DownscaleImage(referenceImage, tmpReferenceImage,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(movingImage, tmpMovingImage,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(referenceMask, tmpReferenceMask,factor);
    btk::PandoraBoxImageFilters::DownscaleImage(movingMask, tmpMovingMask,factor);

    btk::ImageHelper<itkFloatImage>::WriteImage(tmpReferenceImage, "tmpref.nii.gz");
    btk::ImageHelper<itkFloatImage>::WriteImage(tmpMovingImage, "tmpmov.nii.gz");
    btk::ImageHelper<itkFloatImage>::WriteImage(tmpReferenceMask, "tmprefmask.nii.gz");
    btk::ImageHelper<itkFloatImage>::WriteImage(tmpMovingMask, "tmpmovmask.nii.gz");



    myCostFunction->SetReferenceImage(tmpReferenceImage);
    myCostFunction->SetMovingImage(tmpMovingImage);
    myCostFunction->SetMovingMask(tmpMovingMask);
    myCostFunction->SetReferenceMask(tmpReferenceMask);
    myCostFunction->SetCenter(center);

    if( (similarity==1) || (similarity==2) || (similarity==4) || (similarity==5) )
        myCostFunction->InitializeJointHistogram(bin,bin);

    vnl_vector< double > tmpParameterRange = finalParameterRange;
    vnl_vector< double > tmpTolerance      = finalTolerance;

    //update accordingly to image resolution
    for(unsigned int i=3; i < 6; i++)
    {
        tmpParameterRange(i) = finalParameterRange(i) * factor / 2;
    }
    for(unsigned int i=0; i < 6; i++)
    {
        tmpTolerance(i) = finalTolerance(i) * factor;
    }


    //Generate parameters (using uniform sampling over rotation domain)
    std::vector< vnl_vector< double > > candidates;
    vnl_vector< double > rotationRange(3,10);
    int samplingRate = 3;
    btk::PandoraBoxRegistrationFilters::GenerateUniformlyDistributedParameters(inputParam, candidates, rotationRange, samplingRate);

    std::vector< std::pair< vnl_vector< double >, double > > pairOfCandidates(candidates.size());

    std::cout<<"Run registration for each set of candidate parameters\n";
    for(unsigned int i=0; i < candidates.size(); i++)
    {
        btk::PandoraBoxRegistrationFilters::Register3DImages(*myCostFunction, candidates[i], outputParam, tmpParameterRange, tmpTolerance);
        //std::cout<<"Estimated parameters : "<<outputParam<<std::endl;
        pairOfCandidates[i].first = outputParam;
        pairOfCandidates[i].second = (*myCostFunction)(outputParam);
        //A cost function value equal to 0 is not possible (meaning that there is no more overlap between the two images.
        if(pairOfCandidates[i].second == 0)
            pairOfCandidates[i].second = std::numeric_limits<double>::max();
        //std::cout<<"Cost function : "<< pairOfCandidates[i].second <<std::endl;
    }

    std::cout<<"Sort results by cost function values\n";
    std::sort( pairOfCandidates.begin(), pairOfCandidates.end(), less_than_second );
    for(unsigned int i=0; i < candidates.size(); i++)
    {
        std::cout<<"param : "<<pairOfCandidates[i].first<<", cost : "<<pairOfCandidates[i].second<<std::endl;
    }
    std::vector< vnl_vector< double > > bestCandidates(numberOfBestCandidates);
    for(unsigned int i=0; i < numberOfBestCandidates; i++)
        bestCandidates[i] = pairOfCandidates[i].first;




    std::vector< std::pair< vnl_vector< double >, double > > pairOfParamsAndCostFunctionValues(numberOfBestCandidates * numberOfPerturbations);

    //Generate random registration parameters
    srand (time(NULL));
    int index = 0;
    for(unsigned int i=0; i < numberOfBestCandidates; i++)
        for(unsigned int j=0; j < numberOfPerturbations; j++)
        {
            btk::PandoraBoxRegistrationFilters::GenerateRandomParameters(bestCandidates[i], pairOfParamsAndCostFunctionValues[index].first, tmpParameterRange);
            index++;
        }

    //Run registration for each set of parameters
    for(unsigned int i=0; i < numberOfBestCandidates * numberOfPerturbations; i++)
    {
        btk::PandoraBoxRegistrationFilters::Register3DImages(*myCostFunction, pairOfParamsAndCostFunctionValues[i].first, outputParam, tmpParameterRange, tmpTolerance);
        std::cout<<"Estimated parameters : "<<outputParam<<std::endl;
        pairOfParamsAndCostFunctionValues[i].first = outputParam;
        pairOfParamsAndCostFunctionValues[i].second = (*myCostFunction)(outputParam);
        //A cost function value equal to 0 is not possible (meaning that there is no more overlap between the two images.
        if(pairOfParamsAndCostFunctionValues[i].second == 0)
            pairOfParamsAndCostFunctionValues[i].second = std::numeric_limits<double>::max();
        std::cout<<"Cost function : "<< pairOfParamsAndCostFunctionValues[i].second <<std::endl;
    }

    //Sort results by cost function values
    std::sort( pairOfParamsAndCostFunctionValues.begin(), pairOfParamsAndCostFunctionValues.end(), less_than_second );
    for(unsigned int i=0; i < numberOfPerturbations; i++)
        std::cout<<"param : "<<pairOfParamsAndCostFunctionValues[i].first<<", cost : "<<pairOfParamsAndCostFunctionValues[i].second<<std::endl;

    //Assign best estimations
    outputParam = pairOfParamsAndCostFunctionValues[0].first;
*/
    std::cout<<"*********************************************************************\n";
    std::cout<<"Estimated parameters : "<<outputParam<<std::endl;
    std::cout<<"*********************************************************************\n";

    std::cout<<"Convert parameters to ITK matrix\n";
    itkTransformType::Pointer estimatedTransform = itkTransformType::New();
    
    btk::PandoraBoxTransform::ConvertParametersToMatrix(estimatedTransform, outputParam, center);
    std::cout<<"Estimated parameters : "<<outputParam<<std::endl;

    std::cout<<"Resampling\n";
    itkFloatImage::Pointer registeredImage = itkFloatImage::New();
    btk::PandoraBoxImageFilters::ResampleImageUsingTransform(movingImage, registeredImage, referenceImage, interpolationOrder, estimatedTransform);

    std::cout<<"Writing\n";
    if(registered_image_filename != "")
        btk::ImageHelper<itkFloatImage>::WriteImage(registeredImage, registered_image_filename);
    if(output_transform_filename != "")
        btk::IOTransformHelper< itkTransformType >::WriteTransform(estimatedTransform,output_transform_filename);

    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    
    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }
  
}
