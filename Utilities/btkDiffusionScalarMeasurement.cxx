/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 15/05/2013
  Author(s): Julien Pontabry (pontabry@unistra.fr)

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

// TCLAP includes
#include <tclap/CmdLine.h>

// STL includes
#include "string"

// ITK includes
#include "itkImage.h"
#include "itkImageRegionIterator.h"

// BTK includes
#include "btkMacro.h"
#include "btkDiffusionSequence.h"
#include "btkDiffusionSequenceHelper.h"
#include "btkImageHelper.h"
#include "btkSphericalHarmonicsDiffusionDecompositionFilter.h"
#include "btkOrientationDiffusionFunctionModel.h"


const unsigned int Dimension = 3;
typedef double Scalar;
typedef itk::Image< Scalar,Dimension >          ScalarImage;
typedef itk::ImageRegionIterator< ScalarImage > ScalarImageIterator;
typedef btk::SphericalHarmonicsDiffusionDecompositionFilter::OutputImageType CoefficientImage;
typedef itk::ImageRegionIterator< CoefficientImage > CoefficientImageIterator;
typedef itk::Image< short,Dimension > MaskImage;
typedef itk::ImageRegionIterator< MaskImage > MaskImageIterator;


double gentr(std::vector< float > &response)
{
    double coefficient = 3.0 / (2.0 * M_PI);

    double sum = 0.0;

    for(unsigned int i = 0; i < response.size(); i++)
    {
        sum += response[i];
    }

    return coefficient * sum;
}


int main(int argc, char * argv[])
{
    try
    {

        //
        // Command line parser
        //

        // Command line
        TCLAP::CmdLine cmd("Computes diffusion scalar measurements", ' ', "Unversioned");

        // Arguments
        TCLAP::ValueArg< std::string > inputFileNameArg("d", "input", "DWI sequence.", true, "", "string", cmd);
        TCLAP::ValueArg< std::string > outputFileNameArg("o", "output", "Output scalar volume", true, "", "string", cmd);

        // Options
        TCLAP::ValueArg< std::string > maskFileNameArg("m", "mask", "Mask filename", false, "", "string", cmd);
        TCLAP::ValueArg< std::string > scalarArg("", "scalar", "Scalar measurement to compute (FMI, R0, R2, Rm, GA, MSD, GFA) (default: GFA)", false, "GFA", "string", cmd);
        TCLAP::ValueArg< unsigned int > shOrderArg("", "sh_order", "Spherical harmonics order (default: 4)", false, 4, "positive even integer", cmd);

        // Parse arguments
        cmd.parse(argc, argv);

        // Get arguments'values back
        std::string inputFileName  = inputFileNameArg.getValue();
        std::string outputFileName = outputFileNameArg.getValue();

        std::string maskFileName = maskFileNameArg.getValue();
        std::string scalar       = scalarArg.getValue();
        unsigned int shOrder     = shOrderArg.getValue();


        //
        // Testing parameters
        //

        // Test scalar
        if(scalar != "FMI" && scalar != "R0" && scalar != "R2" && scalar != "Rm" && scalar != "GA" && scalar != "MSD" && scalar != "GFA")
        {
            throw(std::string("Scalar is unknown !"));
        }


        //
        // Read sequence
        //

        // Read sequence
        btk::DiffusionSequence::Pointer sequence = btk::DiffusionSequenceHelper::ReadSequence(inputFileName);

        // Read mask
        MaskImage::Pointer mask = NULL;

        if(!maskFileName.empty())
        {
            mask = btk::ImageHelper< MaskImage >::ReadImage(maskFileName);
        }


        //
        // Preprocessing
        //

        std::cout << "Preprocessing..." << std::flush;

        // Compute spherical harmonics coefficients
        btk::SphericalHarmonicsDiffusionDecompositionFilter::Pointer shFilter = btk::SphericalHarmonicsDiffusionDecompositionFilter::New();
        shFilter->SetInput(sequence);
        shFilter->SetSphericalHarmonicsOrder(shOrder);

        if(scalar == "GFA")
        {
            shFilter->SetEstimationType(btk::SphericalHarmonicsDiffusionDecompositionFilter::DIFFUSION_SIGNAL);
        }
        else // scalar != "GFA"
        {
            shFilter->SetEstimationType(btk::SphericalHarmonicsDiffusionDecompositionFilter::APPARENT_DIFFUSION_PROFILE);
        }

        shFilter->Update();

        CoefficientImage::Pointer adc = shFilter->GetOutput();

        // Compute the model function
        btk::OrientationDiffusionFunctionModel::Pointer adcModel = btk::OrientationDiffusionFunctionModel::New();
        adcModel->SetInputModelImage(adc);
        adcModel->Update();

        // Get the number of coefficients
        unsigned int numberOfShCoefficients = adc->GetNumberOfComponentsPerPixel();

        // Clean memory
        sequence = NULL;

        std::cout << "done." << std::endl;


        //
        // Processing
        //

        std::cout << "Processing..." << std::flush;

        // Create new image
        ScalarImage::Pointer output = btk::ImageHelper< CoefficientImage,ScalarImage >::CreateNewImageFromPhysicalSpaceOf(adc);

        // Create mask if needed
        if(maskFileName.empty())
        {
            mask = btk::ImageHelper< ScalarImage,MaskImage >::CreateNewImageFromPhysicalSpaceOf(output, itk::NumericTraits< MaskImage::PixelType >::OneValue());
        }

        // Iterate over image to compute scalar measurements
        ScalarImageIterator outIt(output, output->GetLargestPossibleRegion());
        CoefficientImageIterator inIt(adc, adc->GetLargestPossibleRegion());
        MaskImageIterator mIt(mask, mask->GetLargestPossibleRegion());

        if(scalar == "FMI")
        {
            for(inIt.GoToBegin(), outIt.GoToBegin(), mIt.GoToBegin(); !inIt.IsAtEnd() && !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++inIt, ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get sh coefficients
                    CoefficientImage::PixelType inPixel = inIt.Get();

                    // Sum of the squared sh coefficient with order greater or equal to 4
                    double numerator = 0.0;

                    for(unsigned int i = 6; i < numberOfShCoefficients; i++)
                    {
                        numerator += inPixel[i]*inPixel[i];
                    }

                    // Sum of the squared sh coefficient with order equal to 2
                    double denominator = 0.0;

                    for(unsigned int i = 1; i < 6; i++)
                    {
                        denominator += inPixel[i]*inPixel[i];
                    }

                    // Compute the FMI scalar measurement
                    double FMI = numerator / denominator;

                    // Set output value
                    outIt.Set(FMI);
                }
                else
                {
                    outIt.Set(0.0);
                }
            }
        }
        else if(scalar == "R0")
        {
            for(inIt.GoToBegin(), outIt.GoToBegin(), mIt.GoToBegin(); !inIt.IsAtEnd() && !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++inIt, ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get sh coefficients
                    CoefficientImage::PixelType inPixel = inIt.Get();

                    // Sum of the absolute sh coefficient
                    double denominator = 0.0;

                    for(unsigned int i = 0; i < numberOfShCoefficients; i++)
                    {
                        denominator += std::abs(inPixel[i]);
                    }

                    // Compute the FMI scalar measurement
                    double R0 = std::abs(inPixel[0]) / denominator;

                    // Set output value
                    outIt.Set(R0);
                }
                else
                {
                    outIt.Set(0.0);
                }
            }
        }
        else if(scalar == "R2")
        {
            for(inIt.GoToBegin(), outIt.GoToBegin(), mIt.GoToBegin(); !inIt.IsAtEnd() && !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++inIt, ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get sh coefficients
                    CoefficientImage::PixelType inPixel = inIt.Get();

                    // Sum of the absolute sh coefficient with order equal to 2
                    double numerator = 0.0;

                    for(unsigned int i = 1; i < 6; i++)
                    {
                        numerator += std::abs(inPixel[i]);
                    }

                    // Sum of the absolute sh coefficient
                    double denominator = 0.0;

                    for(unsigned int i = 0; i < numberOfShCoefficients; i++)
                    {
                        denominator += std::abs(inPixel[i]);
                    }

                    // Compute the FMI scalar measurement
                    double R2 = numerator / denominator;

                    // Set output value
                    outIt.Set(R2);
                }
                else
                {
                    outIt.Set(0.0);
                }
            }
        }
        else if(scalar == "Rm")
        {
            for(inIt.GoToBegin(), outIt.GoToBegin(), mIt.GoToBegin(); !inIt.IsAtEnd() && !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++inIt, ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get sh coefficients
                    CoefficientImage::PixelType inPixel = inIt.Get();

                    // Sum of the absolute sh coefficient with order equal to 2
                    double numerator = 0.0;

                    for(unsigned int i = 6; i < numberOfShCoefficients; i++)
                    {
                        numerator += std::abs(inPixel[i]);
                    }

                    // Sum of the absolute sh coefficient
                    double denominator = 0.0;

                    for(unsigned int i = 0; i < numberOfShCoefficients; i++)
                    {
                        denominator += std::abs(inPixel[i]);
                    }

                    // Compute the FMI scalar measurement
                    double Rm = numerator / denominator;

                    // Set output value
                    outIt.Set(Rm);
                }
                else
                {
                    outIt.Set(0.0);
                }
            }
        }
        else if(scalar == "GA")
        {
            double min = DBL_MAX, max = DBL_MIN;

            for(outIt.GoToBegin(), mIt.GoToBegin(); !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get current index
                    ScalarImage::IndexType index = outIt.GetIndex();

                    btk::OrientationDiffusionFunctionModel::ContinuousIndex cindex;
                    cindex[0] = index[0]; cindex[1] = index[1]; cindex[2] = index[2];

                    // Get ADC profile
                    std::vector< float > adcResponse = adcModel->SignalAt(cindex);

                    // Normalize coefficients
                    double normCoefficient = gentr(adcResponse);

                    for(unsigned int i = 0; i < adcResponse.size(); i++)
                    {
                        adcResponse[i] /= normCoefficient;
                        adcResponse[i] *= adcResponse[i];
                    }

                    // Compute variance (approximation of GA)
                    double GA = (1.0/3.0) * (gentr(adcResponse) - (1.0/3.0));

                    if(std::isnan(GA))
                    {
                        GA = 0.0;
                    }

                    // Set output value
                    outIt.Set(GA);

                    if(GA < min)
                    {
                        min = GA;
                    }

                    if(GA > max)
                    {
                        max = GA;
                    }
                }
                else
                {
                    outIt.Set(0.0);
                }
            }

            double maxMinusMin = max - min;

            for(outIt.GoToBegin(), mIt.GoToBegin(); !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    outIt.Set(((outIt.Get() - min) / maxMinusMin)*1000.0);
                }
            }
        }
        else if(scalar == "MSD")
        {
            for(outIt.GoToBegin(), mIt.GoToBegin(); !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get current index
                    ScalarImage::IndexType index = outIt.GetIndex();

                    btk::OrientationDiffusionFunctionModel::ContinuousIndex cindex;
                    cindex[0] = index[0]; cindex[1] = index[1]; cindex[2] = index[2];

                    // Get ADC profile
                    std::vector< float > adcResponse = adcModel->SignalAt(cindex);

                    // Compute variance
                    double MSD = (1.0/3.0) * (gentr(adcResponse));

                    if(std::isnan(MSD))
                    {
                        MSD = 0.0;
                    }

                    // Set output value
                    outIt.Set(MSD);
                }
                else
                {
                    outIt.Set(0.0);
                }
            }
        }
        else if(scalar == "GFA")
        {
            for(outIt.GoToBegin(), mIt.GoToBegin(); !outIt.IsAtEnd() && !mIt.IsAtEnd(); ++outIt, ++mIt)
            {
                if(mIt.Get() != 0)
                {
                    // Get current index
                    ScalarImage::IndexType index = outIt.GetIndex();

                    btk::OrientationDiffusionFunctionModel::ContinuousIndex cindex;
                    cindex[0] = index[0]; cindex[1] = index[1]; cindex[2] = index[2];

                    // Get ODF
                    std::vector< float > response = adcModel->ModelAt(cindex);

                    // Compute the mean and the mean squarred
                    double mean = 0.0, meanSq = 0.0;

                    for(unsigned int i = 0; i < response.size(); i++)
                    {
                        mean   += response[i];
                        meanSq += response[i] * response[i];
                    }

                    mean /= response.size();

                    // Compute the variance
                    double variance = 0.0;

                    for(unsigned int i = 0; i < response.size(); i++)
                    {
                        double deviation = response[i] - mean;
                        variance       += deviation * deviation;
                    }

                    // Compute the GFA
                    double GFA = std::sqrt( (response.size() * variance) / ((response.size()-1) * meanSq) );

                    if(std::isnan(GFA))
                    {
                        GFA = 0.0;
                    }

                    // Set output value
                    outIt.Set(GFA);
                }
                else
                {
                    outIt.Set(0.0);
                }
            }
        }
        else // scalar == anything else
        {
            throw(std::string("Scalar is unknown !"));
        }

        std::cout << "done." << std::endl;


        //
        // Write output
        //

        // Write modified sequence
        btk::ImageHelper< ScalarImage >::WriteImage(output, outputFileName);
    }
    catch(TCLAP::ArgException &e)
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch(std::string &e)
    {
        std::cerr << "Error: " << e << std::endl;
    }

    return EXIT_SUCCESS;
}
