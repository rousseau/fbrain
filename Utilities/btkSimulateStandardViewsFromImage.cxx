/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 20/12/2010
  Author(s): Estanislao Oubel (oubel@unistra.fr)

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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "btkOrientedSpatialFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkContinuousIndex.h"


#include "tclap/CmdLine.h"


int main( int argc, char *argv[] )
{

  try {

  std::vector< std::string > input;
  std::vector< std::string > output;

  float sx, sy, sz;

  const char *inputName = NULL;
  const char *outputName = NULL;

  // Parse arguments

  TCLAP::CmdLine cmd("Compute a low resolution image from a high resolution image by using a generative model", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","High resolution image",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Low resolution image",true,"","string",cmd);

  TCLAP::ValueArg<float> sxArg  ("","sx","Spacing in x-direction",false,1,"float",cmd);
  TCLAP::ValueArg<float> syArg  ("","sy","Spacing in y-direction",false,1,"float",cmd);
  TCLAP::ValueArg<float> szArg  ("","sz","Spacing in z-direction",false,1,"float",cmd);

  TCLAP::SwitchArg axlSwitch("","axl","Lower resolution in x", false);
  TCLAP::SwitchArg corSwitch("","cor","Lower resolution in y", false);
  TCLAP::SwitchArg sagSwitch("","sag","Lower resolution in z", false);

  std::vector<TCLAP::Arg*>  xorlist;
  xorlist.push_back(&axlSwitch);
  xorlist.push_back(&corSwitch);
  xorlist.push_back(&sagSwitch);

  cmd.xorAdd( xorlist );

  // Parse the argv array.
  cmd.parse( argc, argv );

  inputName  = inputArg.getValue().c_str();
  outputName = outputArg.getValue().c_str();

  sx = sxArg.getValue();
  sy = syArg.getValue();
  sz = szArg.getValue();

  bool axl = axlSwitch.getValue();
  bool cor = corSwitch.getValue();
  bool sag = sagSwitch.getValue();


  typedef short PixelType;
  const   unsigned int    Dimension = 3;

  typedef itk::Image< PixelType, Dimension >  ImageType;

  typedef itk::ImageFileReader< ImageType >   ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >   WriterType;


  // Read reference image

  ImageReaderType::Pointer hrReader = ImageReaderType::New();
  hrReader -> SetFileName( inputName );
  hrReader -> Update();

  ImageType::Pointer hrImage = hrReader -> GetOutput();


  // Modify direction if necessary (cor or sag = true)

  vnl_matrix_fixed<double,3,3> vnlDirection  = hrImage -> GetDirection().GetVnlMatrix();
  vnl_matrix_fixed<double,3,3> vnlDirectionMod = hrImage -> GetDirection().GetVnlMatrix();

  if ( cor )
  {
    vnlDirectionMod.set_column(1,vnlDirection.get_column(2));
    vnlDirectionMod.set_column(2,vnlDirection.get_column(1));
  } else if ( sag )
    {
      vnlDirectionMod.set_column(0,vnlDirection.get_column(2));
      vnlDirectionMod.set_column(2,vnlDirection.get_column(0));
    }

  ImageType::DirectionType newDirection;
  newDirection = vnlDirectionMod;

//  hrImage -> SetDirection ( newDirection );


  // Create low resolution image

  ImageType::SpacingType hrSpacing = hrImage -> GetSpacing();
  ImageType::SizeType    hrSize    = hrImage -> GetLargestPossibleRegion().GetSize();

  ImageType::SpacingType lrSpacing;
  lrSpacing[0] = sx; lrSpacing[1] = sy; lrSpacing[2] = sz;

  ImageType::Pointer lrImage = ImageType::New();

  ImageType::RegionType lrRegion;

  ImageType::IndexType  lrStart;
  lrStart[0] = 0; lrStart[1] = 0; lrStart[2] = 0;

  ImageType::SizeType   lrSize;

  if ( axl )
  {
    lrSize[0] = floor(hrSize[0] * hrSpacing[0] / lrSpacing[0]);
    lrSize[1] = floor(hrSize[1] * hrSpacing[1] / lrSpacing[1]);
    lrSize[2] = floor(hrSize[2] * hrSpacing[2] / lrSpacing[2]);
  } else if ( cor )
    {
      lrSize[0] = floor(hrSize[0] * hrSpacing[0] / lrSpacing[0]);
      lrSize[1] = floor(hrSize[2] * hrSpacing[2] / lrSpacing[1]);
      lrSize[2] = floor(hrSize[1] * hrSpacing[1] / lrSpacing[2]);
    } else if ( sag )
      {
        lrSize[0] = floor(hrSize[2] * hrSpacing[2] / lrSpacing[0]);
        lrSize[1] = floor(hrSize[1] * hrSpacing[1] / lrSpacing[1]);
        lrSize[2] = floor(hrSize[0] * hrSpacing[0] / lrSpacing[2]);
      }


  lrRegion.SetIndex( lrStart );
  lrRegion.SetSize(  lrSize );

  lrImage -> SetRegions( lrRegion );
  lrImage -> Allocate();

  lrImage -> SetSpacing( lrSpacing );
  lrImage -> SetDirection( newDirection );
  lrImage -> SetOrigin( hrImage -> GetOrigin() + (lrSpacing - hrSpacing)/ 2.0 );


  typedef btk::OrientedSpatialFunction<double, Dimension, ImageType::PointType> FunctionType;
  FunctionType::Pointer function = FunctionType::New();
  function -> SetDirection( lrImage -> GetDirection() );
  function -> SetSpacing(   lrImage -> GetSpacing() );

  typedef itk::ContinuousIndex<double, Dimension> ContinuousIndexType;

  std::vector<ContinuousIndexType> deltaIndexes;
  int npoints =  lrSpacing[2] / (2.0 * hrSpacing[2]) ;
  ContinuousIndexType delta;
  delta[0] = 0.0; delta[1] = 0.0;

  for (int i = -npoints ; i <= npoints; i++ )
  {
    delta[2] = i * 0.5 / (double) npoints;
    deltaIndexes.push_back(delta);

//    std::cout << delta << std::endl;
  }

  typedef itk::LinearInterpolateImageFunction< ImageType, double> InterpolatorType;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator -> SetInputImage( hrImage );

  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  IteratorType lrIt( lrImage, lrImage -> GetLargestPossibleRegion() );

  ImageType::IndexType lrIndex;
  ImageType::PointType lrPoint;
  ContinuousIndexType  nbIndex;
  ImageType::PointType nbPoint;

  unsigned int counter;
  double value;

  for ( lrIt.GoToBegin(); ! lrIt.IsAtEnd(); ++lrIt)
  {
    lrIndex = lrIt.GetIndex();
    lrImage -> TransformIndexToPhysicalPoint( lrIndex, lrPoint );
    function -> SetCenter( lrPoint );

    counter = 0;
    value   = 0.0;

    for(unsigned int k=0; k<deltaIndexes.size(); k++)
    {
      nbIndex[0] = deltaIndexes[k][0] + lrIndex[0];
      nbIndex[1] = deltaIndexes[k][1] + lrIndex[1];
      nbIndex[2] = deltaIndexes[k][2] + lrIndex[2];

      lrImage -> TransformContinuousIndexToPhysicalPoint( nbIndex, nbPoint );

      if ( function -> Evaluate( nbPoint ) > 0)
      {
        // I assume an identity transform here ...

        if ( interpolator -> IsInsideBuffer( nbPoint ) )
        {
          value += interpolator -> Evaluate( nbPoint );
          counter++;

        }

      }

    }

    lrIt.Set(value / counter );
//    lrIt.Set( 1 );

  }


  // Write image

  WriterType::Pointer writer =  WriterType::New();
  writer -> SetFileName( outputName );
  writer -> SetInput( lrImage );

  if ( strcmp(outputName,"") != 0)
  {
    std::cout << "Writing " << outputName << " ... ";
    writer -> Update();
    std::cout << "done." << std::endl;
  }

  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

  return EXIT_SUCCESS;
}

