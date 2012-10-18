/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

16 march 2011
rousseau @ unistra . fr

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
#include "string"
#include "iomanip"

/* Itk includes */
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkCropImageFilter.h"
#include "itkImageIOBase.h"

/* Btk includes */
#include "btkImageHelper.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkIOImageHelper.h"
#include "btkResampleImagesToBiggestImageFilter.h"
#include "btkWeightedSumOfImagesFilter.h"
#include "btkCropImageUsingMaskFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"


template< class TImage,class TMask,int Dimension >
void Process(std::vector< std::string > &inputFileNames, std::vector< std::string > &outputFileNames, std::vector< std::string > &maskFileNames)
{
    // Read image files
    std::vector< typename TImage::Pointer > inputs = /*btk::DiffusionSequenceHelper::ReadSequenceArray(inputFileNames):*/ btk::ImageHelper< TImage >::ReadImage(inputFileNames);
    std::vector< typename TMask::Pointer >   masks = btk::ImageHelper< TMask >::ReadImage(maskFileNames);


    // Compute the cropping mask by using the union of all masks
    typename TMask::Pointer croppingMask = NULL;

    if(masks.size() > 1)
    {
        // First, resample masks if needed
        if(!btk::ImageHelper< TMask >::IsInSamePhysicalSpace(masks))
        {
            btkCoutMacro("Resampling masks in same physical space... ");

            typename btk::ResampleImagesToBiggestImageFilter< TMask >::Pointer resampleFilter = btk::ResampleImagesToBiggestImageFilter< TMask >::New();
            resampleFilter->SetInputs(masks);
            resampleFilter->SetInterpolator(itk::NearestNeighborInterpolateImageFunction< TMask >::New());
            resampleFilter->Update();
            masks = resampleFilter->GetOutputs();

            btkCoutMacro("done.");
        }

        // Combine the masks
        btkCoutMacro("Combining masks... ");

        std::vector< float > weights(masks.size(), 1.0f);

        typename btk::WeightedSumOfImagesFilter< TMask >::Pointer sumFilter = btk::WeightedSumOfImagesFilter< TMask >::New();
        sumFilter->SetInputs(masks);
        sumFilter->SetWeights(weights);
        sumFilter->Update();

        croppingMask = sumFilter->GetOutput();
//btk::ImageHelper< TMask >::WriteImage(croppingMask, "test.nii.gz");
        btkCoutMacro("done.");
    }
    else // masks.size() == 1
    {
        croppingMask = masks[0];
    }


    // Crop input images
    btkCoutMacro("Cropping images... ");

    typename btk::CropImageUsingMaskFilter< TImage,TMask >::Pointer cropFilter = btk::CropImageUsingMaskFilter< TImage,TMask >::New();
    cropFilter->SetMask(croppingMask);
    cropFilter->SetInputs(inputs);
    cropFilter->Update();

    std::vector< typename TImage::Pointer > outputs = cropFilter->GetOutputs();

    btkCoutMacro("done.");


    // Write output files
    btk::ImageHelper< TImage >::WriteImage(outputs, outputFileNames);
}


int main(int argc, char** argv)
{
    try
    {
        //
        // Command line parser
        //

        TCLAP::CmdLine cmd("btkCropImageUsingMask: Crop an image (3D - anatomical, or 4D - diffusion weighted) using a 3D mask (non-zero values). Input image and mask image must have the same pixel type.", ' ', "2.0", true);

        TCLAP::MultiArg< std::string > inputFileNamesArg("i", "input", "Input image file", true, "string", cmd);
        TCLAP::MultiArg< std::string > maskFileNamesArg("m", "mask", "Mask image file", true, "string", cmd);
        TCLAP::MultiArg< std::string > outputFileNamesArg("o", "output", "Output image file", true, "string", cmd);

        // Parse the args.
        cmd.parse( argc, argv );

        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();
        std::vector< std::string > maskFileNames = maskFileNamesArg.getValue();
        std::vector< std::string > outputFileNames = outputFileNamesArg.getValue();


        //
        // Processing
        //

        // Check arguments
        if(inputFileNames.size() != maskFileNames.size() || inputFileNames.size() != outputFileNames.size())
        {
            btkException("The number of inputs/masks/outputs is not coherent ! There should be a mask and an output for each input !");
        }

        // Check image informations
        itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFileNames[0].c_str(), itk::ImageIOFactory::ReadMode);
        imageIO->SetFileName(inputFileNames[0]);
        imageIO->ReadImageInformation();

        if(imageIO->GetPixelType() != itk::ImageIOBase::SCALAR)
        {
            btkException("Unsupported image pixel type (only scalar are supported) !");
        }

        if(imageIO->GetNumberOfDimensions() == 3)
        {
            switch(imageIO->GetComponentType())
            {
                case itk::ImageIOBase::SHORT:
                    Process< itk::Image< short,3 >,itk::Image< unsigned char,3 >,3 >(inputFileNames, outputFileNames, maskFileNames);
                    break;

                case itk::ImageIOBase::FLOAT:
                    Process< itk::Image< float,3 >,itk::Image< unsigned char,3 >,3 >(inputFileNames, outputFileNames, maskFileNames);
                    break;

                case itk::ImageIOBase::DOUBLE:
                    Process< itk::Image< double,3 >,itk::Image< unsigned char,3 >,3 >(inputFileNames, outputFileNames, maskFileNames);
                    break;

                default:
                    btkException("Unsupported component type (only short, float or double types are supported) !");
            }
        }
        else if(imageIO->GetNumberOfDimensions() == 4)
        {
//            switch(imageIO->GetComponentType())
//            {
//                case itk::ImageIOBase::SHORT:
//                    Process< btk::DiffusionSequence,itk::Image< unsigned char,3 >,4 >(inputFileNames, outputFileNames, maskFileNames);
//                    break;

//                case itk::ImageIOBase::FLOAT:
//                    Process< btk::DiffusionSequence,itk::Image< unsigned char,3 >,4 >(inputFileNames, outputFileNames, maskFileNames);
//                    break;

//                case itk::ImageIOBase::DOUBLE:
//                    Process< btk::DiffusionSequence,itk::Image< unsigned char,3 >,4 >(inputFileNames, outputFileNames, maskFileNames);
//                    break;

//                default:
//                    btkException("Unsupported component type (only short, float or double types are supported) !");
//            }
            btkException("Dimension 4 is not implemented yet !");
        }
        else // imageIO->GetNumberOfDimensions() != 3 && imageIO->GetNumberOfDimensions() != 4
        {
            btkException("Unsupported dimension (only anatomical - 3D - or diffusion weighted - 4D - images supported) !");
        }
    }
    catch(TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch(std::string &e)  // catch any exceptions
    {
        std::cerr << "error: " << e <<std::endl;
    }
    catch(itk::ExceptionObject &e)
    {
        btkCerrMacro(e.GetDescription());
    }

    return EXIT_SUCCESS;
}




