/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/11/2010
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

#include <iostream>
#include <string>

#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "CmdLine.h"

int main( int argc, char * argv[] )
{
  try{

  const char *outputFile = NULL;

  // Parse arguments

  TCLAP::CmdLine cmd("Converts a dwi sequence in nifti format to nrrd format", ' ', "Unversioned");

  TCLAP::ValueArg<std::string> inputArg("i","input","Diffusion-weighted sequence",true,"","string",cmd);
  TCLAP::ValueArg<std::string> outputArg("o","output","Sequence in NRRD format",true,"","string",cmd);

  TCLAP::SwitchArg  lsssSwitchArg("v","vectorImage","The output is volume interleaved (an image of "
  "vectors of diffusion values). By default, the output is voxel interleaved (you go through the "
  "DWI values for each voxel before moving on to the values for the next voxel). ",cmd,false);

  TCLAP::SwitchArg  lpsSwitchArg("","lps","Word coordinates expressed in LPS (Left-Posterior-Superior). "
  "By default RAS (Right-Anterior-Superior) is used.",cmd,false);

  cmd.parse( argc, argv );

  outputFile = outputArg.getValue().c_str();

  std::string rawFile = outputArg.getValue();
  rawFile.replace(rawFile.size()-4,7,"raw.gz");

  char inputFile[255];
  strcpy( inputFile, (char*)inputArg.getValue().c_str() );
  strcat ( inputFile,".nii" );

  char bvecFile[255];
  strcpy( bvecFile, (char*)inputArg.getValue().c_str() );
  strcat ( bvecFile,".bvec" );

  char bvalFile[255];
  strcpy( bvalFile, (char*)inputArg.getValue().c_str() );
  strcat ( bvalFile,".bval" );


  // Read dwi sequence

  typedef short         PixelType;
  const   unsigned int  Dimension = 4;

  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >    ReaderType;

  ReaderType::Pointer reader = ReaderType::New();

  reader -> SetFileName( inputFile );
  reader -> Update();

  ImageType::Pointer image = reader -> GetOutput();

  ImageType::SizeType       imageSize       = image -> GetLargestPossibleRegion().GetSize();
  ImageType::PointType      imageOrigin     = image -> GetOrigin();
  ImageType::SpacingType    imageSpacing    = image -> GetSpacing();
  ImageType::DirectionType  imageDirection  = image -> GetDirection();

// Read b-values

  std::vector< int > bValues;
  std :: cout << "Reading bvalues from "  << bvalFile << std::endl; std::cout.flush();

  FILE* f;

  f = fopen( bvalFile, "r" );
  int bvalue;
  while ( !feof(f) )
  {
    fscanf( f, "%d ", &bvalue);
    bValues.push_back(bvalue);
  }
  fclose (f);

// Read gradients

  std::vector< std::vector< float > > gValues(3);

  std :: cout << "Reading gradients from " << bvecFile << std::endl; std::cout.flush();

  f = fopen( bvecFile, "r" );
  float gvalue;
  int coor = 0;
  for (unsigned int i=1; i<=3*bValues.size(); i++)
  {
    if ( ( i % bValues.size() ) == 0 )
      fscanf( f, "%f\n", &gvalue);
    else
      fscanf( f, "%f ", &gvalue);

    gValues[coor].push_back(gvalue);

    if ( ( i % bValues.size() ) == 0 )
      coor++;
  }
  fclose (f);


  if (lsssSwitchArg.isSet())
  {

  // Create vector image

    typedef itk::VectorImage< PixelType, 3 >   VectorImageType;
    VectorImageType::Pointer nrrdImage = VectorImageType::New();
    nrrdImage -> SetVectorLength( bValues.size() );

    VectorImageType::IndexType  start;
    VectorImageType::SizeType   size;
    VectorImageType::RegionType region;

    start[0] = 0; start[1] = 0; start[2] = 0;
    size[0]  = imageSize[0];  size[1]  = imageSize[1];  size[2]  = imageSize[2];

    region.SetSize( size );
    region.SetIndex( start );

    nrrdImage -> SetRegions( region );
    nrrdImage -> Allocate();


    // Fill vector image

    typedef itk::ImageRegionIteratorWithIndex< VectorImageType >  IteratorType;
    IteratorType nrrdIt( nrrdImage, nrrdImage -> GetLargestPossibleRegion() );

    VectorImageType::IndexType nrrdIndex;
    VectorImageType::PixelType nrrdPixel;
    ImageType::IndexType niftiIndex;

    for( nrrdIt.GoToBegin(); !nrrdIt.IsAtEnd(); ++nrrdIt)
    {
      nrrdIndex = nrrdIt.GetIndex();
      nrrdPixel = nrrdIt.Get();

      niftiIndex[0] = nrrdIndex[0]; niftiIndex[1] = nrrdIndex[1]; niftiIndex[2] = nrrdIndex[2];

      for ( unsigned int k=0; k < nrrdImage -> GetVectorLength(); k++)
      {
        niftiIndex[3] = k;
        nrrdPixel[k] = image -> GetPixel( niftiIndex );

      }

      nrrdIt.Set( nrrdPixel );

    }

    // Write binary data
    // TODO How to compress data? NrrdIO option?

    typedef itk::ImageFileWriter< VectorImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();

    writer -> SetInput( nrrdImage );
    writer -> SetFileName( outputFile );
    writer -> UseCompressionOn();
    writer -> Update();

  }
  else
  {
    // Write binary data
    // TODO How to compress data? NrrdIO option?

    typedef itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();

    writer -> SetInput( image );
    writer -> SetFileName( outputFile );
    writer -> UseCompressionOn();
    writer -> Update();

  }

  // Write header

  FILE* fw;
  fw = fopen( outputFile, "w" );
  fprintf( fw, "NRRD0005\n" );

  char* pch;
  char* rawFile_cstr = (char*)rawFile.c_str();

  pch = strtok(rawFile_cstr,"/");
  while (pch != NULL)
  {
    rawFile_cstr = pch;
    pch = strtok (NULL, "/");
  }

  fprintf( fw, "content: exists(%s,0)\n",rawFile_cstr);
  fprintf( fw, "type: short\n");
  fprintf( fw, "dimension: 4\n");
  //TODO How to change space? NrrdIO option?
  if ( lpsSwitchArg.isSet() )
    fprintf( fw, "space: left-posterior-superior\n");
  else
    fprintf( fw, "space: right-anterior-superior\n");

  if (lsssSwitchArg.isSet())
    fprintf( fw, "sizes: %ld %ld %ld %ld\n", imageSize[3], imageSize[0], imageSize[1], imageSize[2]);
  else
    fprintf( fw, "sizes: %ld %ld %ld %ld\n", imageSize[0], imageSize[1], imageSize[2], imageSize[3]);

  fprintf( fw, "thicknesses: NaN NaN %lf NaN\n", imageSpacing[2]);

  vnl_matrix<double> niftiDir;
  niftiDir = imageDirection.GetVnlMatrix();

  vnl_matrix<double> nrrdDir;
  nrrdDir = niftiDir.extract(3,3);

  if ( !lpsSwitchArg.isSet() )
  {
    nrrdDir.scale_row(0,-1);
    nrrdDir.scale_row(1,-1);
  }

  vnl_matrix<double> nrrdSpa(3,3);
  nrrdSpa.set_identity();
  nrrdSpa(0,0) = imageSpacing[0]; nrrdSpa(1,1) = imageSpacing[1]; nrrdSpa(2,2) = imageSpacing[2];

  vnl_matrix<double> nrrdSpaDir;
  nrrdSpaDir = nrrdDir * nrrdSpa;

  if (lsssSwitchArg.isSet())
    fprintf( fw, "space directions: none (%lf,%lf,%lf) (%lf,%lf,%lf) (%lf,%lf,%lf)\n",
           nrrdSpaDir(0,0), nrrdSpaDir(1,0), nrrdSpaDir(2,0),
           nrrdSpaDir(0,1), nrrdSpaDir(1,1), nrrdSpaDir(2,1),
           nrrdSpaDir(0,2), nrrdSpaDir(1,2), nrrdSpaDir(2,2) );
  else
    fprintf( fw, "space directions: (%lf,%lf,%lf) (%lf,%lf,%lf) (%lf,%lf,%lf) none\n",
           nrrdSpaDir(0,0), nrrdSpaDir(1,0), nrrdSpaDir(2,0),
           nrrdSpaDir(0,1), nrrdSpaDir(1,1), nrrdSpaDir(2,1),
           nrrdSpaDir(0,2), nrrdSpaDir(1,2), nrrdSpaDir(2,2) );


  fprintf( fw, "centerings: cell cell cell ???\n");

  if (lsssSwitchArg.isSet())
    fprintf( fw, "kinds: list space space space\n");
  else
    fprintf( fw, "kinds: space space space list\n");

  fprintf( fw, "endian: little\n");
  fprintf( fw, "encoding: gzip\n");
  fprintf( fw, "space units: \"mm\" \"mm\" \"mm\"\n");

  if ( lpsSwitchArg.isSet() )
    fprintf( fw, "space origin: (%lf,%lf,%lf)\n", imageOrigin[0], imageOrigin[1], imageOrigin[2]);
  else
    fprintf( fw, "space origin: (%lf,%lf,%lf)\n", -imageOrigin[0], -imageOrigin[1], imageOrigin[2]);

  fprintf( fw, "data file: %s\n",rawFile_cstr);
  fprintf( fw, "measurement frame: (1,0,0) (0,1,0) (0,0,1)\n");
  fprintf( fw, "modality:=DWMRI\n");

  // Express gradient directions in WC

  vnl_vector< double > g(3);
  vnl_vector< double > gwc(3);

  for ( unsigned int k=0; k < imageSize[3]; k++ )
  {
    g(0) = gValues[0][k]; g(1) = gValues[1][k]; g(2) = gValues[2][k];
    gwc = nrrdDir * g;
    gValues[0][k] = gwc(0); gValues[1][k] = gwc(1); gValues[2][k] = gwc(2);
  }

  fprintf( fw, "DWMRI_b-value:= %d\n",bValues[bValues.size()-1]);

  for ( unsigned int k=0; k < imageSize[3]; k++)
  {
    if ( lpsSwitchArg.isSet() )
      fprintf( fw, "DWMRI_gradient_%.4d:= %.14lf %.14lf %.14lf\n",k,-gValues[0][k],-gValues[1][k],gValues[2][k]);
    else
      fprintf( fw, "DWMRI_gradient_%.4d:= %.14lf %.14lf %.14lf\n",k,gValues[0][k],gValues[1][k],gValues[2][k]);
  }

  fclose (fw);


  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }


  return EXIT_SUCCESS;
}
