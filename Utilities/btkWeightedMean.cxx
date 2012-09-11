/*==========================================================================
 
 © Université de Strasbourg - Centre National de la Recherche Scientifique
 
 Date: 05/10/2011
 Author(s): François Rousseau (rousseau@unistra.fr)
            Julien Pontabry (pontabry@unistra.fr)
 
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

// STL includes
#include "string"
#include "vector"

// TCLAP includes
#include <tclap/CmdLine.h>

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkDisplacementFieldTransform.h"

// Local includes
#include "btkMacro.h"
#include "btkImageHelper.h"


// ITK definitions
typedef itk::Image< float,3 >                                             Image;
typedef itk::DisplacementFieldTransform< float,3 >::DisplacementFieldType DeformationField;


// Weight sum over template images
template< typename TImage >
void ComputeAndWriteWeightedMean(std::vector< std::string > &inputFileNames, std::vector< float > &weights, const std::string &outputFileName)
{
    //
    // Read & check
    //

    // Verify number of images and number of weights
    if(weights.size() > 0 && inputFileNames.size() != weights.size())
        throw(std::string("The number of input images is different than the number of weights !"));

    // Read images
    std::vector< typename TImage::Pointer > inputImages = btk::ImageHelper< TImage >::ReadImage(inputFileNames);

    // Verify sizes of images
    if(!btk::ImageHelper< TImage >::IsInSamePhysicalSpace(inputImages))
        throw(std::string("Input images are not in the same physical space !"));


    //
    // Compute the weighted mean
    //

    typename TImage::Pointer outputImage = btk::ImageHelper< TImage >::CreateNewImageFromPhysicalSpaceOf(inputImages[0]);

    itk::ImageRegionIterator< TImage > outputIt(outputImage, outputImage->GetLargestPossibleRegion());

    if(weights.size() > 0)
    {
        for(unsigned int i = 0; i < inputImages.size(); i++)
        {
            itk::ImageRegionIterator< TImage > inputIt(inputImages[i], inputImages[i]->GetLargestPossibleRegion());

            for(inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd() && !outputIt.IsAtEnd(); ++inputIt, ++outputIt)
            {
                outputIt.Set( outputIt.Get() + (inputIt.Get() * weights[i]) );
            } // for each voxel
        } // for each image
    }
    else // weights.size() == 0
    {
        float defaultValue = 1.0 / static_cast< float >(inputImages.size());

        for(unsigned int i = 0; i < inputImages.size(); i++)
        {
            itk::ImageRegionIterator< TImage > inputIt(inputImages[i], inputImages[i]->GetLargestPossibleRegion());

            for(inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd() && !outputIt.IsAtEnd(); ++inputIt, ++outputIt)
            {
                outputIt.Set( outputIt.Get() + (inputIt.Get() * defaultValue) );
            } // for each voxel
        } // for each image
    }


    //
    // Write output
    //

    btk::ImageHelper< TImage >::WriteImage(outputImage, outputFileName);
}


// Main function
int main(int argc, char *argv[])
{ 
    try
    {

        //
        // Command line parser
        //

        // Command line
        TCLAP::CmdLine cmd("Compute a weighted mean of 3D images", ' ', "2.0", true);
    
        // Arguments
        TCLAP::MultiArg< std::string > inputFileNamesArg("i", "input", "Input image filenames", true, "string", cmd);
        TCLAP::MultiArg< float >              weightsArg("w", "weight", "Image weight (default: 1/N)", false, "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output image filename (default: \"out.nii.gz\")", false, "out.nii.gz", "string", cmd);

        TCLAP::SwitchArg deformationFieldArg("f", "deformation_field", "Operate on deformation field images", cmd);
    
        // Parse the args.
        cmd.parse(argc, argv);

        // Get the value parsed by each arg.
        std::vector< std::string > inputFileNames = inputFileNamesArg.getValue();
        std::vector< float >              weights = weightsArg.getValue();
        std::string                outputFileName = outputFileNameArg.getValue();

        bool deformationField = deformationFieldArg.getValue();
 

        //
        // Processing
        //

        if(deformationField)
        {
            ComputeAndWriteWeightedMean< DeformationField >(inputFileNames, weights, outputFileName);
        }
        else
        {
            ComputeAndWriteWeightedMean< Image >(inputFileNames, weights, outputFileName);
        }
    }
    catch(TCLAP::ArgException &e)
    {
        btkCoutMacro("Exception: " << e.error() << " for arg " << e.argId());
    }
    catch(std::string &message)
    {
        btkCoutMacro("Exception: " << message);
    }
  
  return EXIT_SUCCESS;
}




