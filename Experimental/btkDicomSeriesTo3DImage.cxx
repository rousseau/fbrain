/*
Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique

31 january 2014
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

/* Itk includes */
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"

int main(int argc, char** argv)
{
  try {
  
    TCLAP::CmdLine cmd("Dummy Dicom converter (GDCM based - ITK4.5 example)", ' ', "0.1", true);

    TCLAP::ValueArg<std::string> inputDirArg ("i","image_dir","input image directory",true,"","string", cmd);

  	// Parse the args.
  	cmd.parse( argc, argv );

  	// Get the value parsed by each arg. 
    std::string input_dir                      = inputDirArg.getValue();



    typedef signed short    PixelType;
    const unsigned int      Dimension = 3;

    typedef itk::Image< PixelType, Dimension >         ImageType;
    typedef itk::ImageSeriesReader< ImageType >        ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    typedef itk::GDCMImageIO       ImageIOType;
    ImageIOType::Pointer dicomIO = ImageIOType::New();

    reader->SetImageIO( dicomIO );

    typedef itk::GDCMSeriesFileNames NamesGeneratorType;
    NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

    nameGenerator->SetUseSeriesDetails( true );
    nameGenerator->AddSeriesRestriction("0008|0021" );

    nameGenerator->SetDirectory( input_dir );

    try
      {
      std::cout << std::endl << "The directory: " << std::endl;
      std::cout << std::endl << input_dir << std::endl << std::endl;
      std::cout << "Contains the following DICOM Series: ";
      std::cout << std::endl << std::endl;


      typedef std::vector< std::string >    SeriesIdContainer;

      const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();

      SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
      SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
      while( seriesItr != seriesEnd )
        {
        std::cout << seriesItr->c_str() << std::endl;
        ++seriesItr;
        }

      std::string seriesIdentifier;

      seriesIdentifier = seriesUID.begin()->c_str();

      std::cout << std::endl << std::endl;
      std::cout << "Now reading series: " << std::endl << std::endl;
      std::cout << seriesIdentifier << std::endl;
      std::cout << std::endl << std::endl;

      typedef std::vector< std::string >   FileNamesContainer;
      FileNamesContainer fileNames;

      fileNames = nameGenerator->GetFileNames( seriesIdentifier );

      reader->SetFileNames( fileNames );

      try
        {
        reader->Update();
        }
      catch (itk::ExceptionObject &ex)
        {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
        }

      typedef itk::MetaDataDictionary   DictionaryType;
      const  DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
      typedef itk::MetaDataObject< std::string > MetaDataStringType;
      DictionaryType::ConstIterator itr = dictionary.Begin();
      DictionaryType::ConstIterator end = dictionary.End();

      std::string acquisition_date =  dynamic_cast<const MetaDataStringType *>( dictionary.Find( "0008|0022" )->second.GetPointer() )->GetMetaDataObjectValue();
      std::string acquisition_time =  dynamic_cast<const MetaDataStringType *>( dictionary.Find( "0008|0032" )->second.GetPointer() )->GetMetaDataObjectValue();
      std::string protocol = dynamic_cast<const MetaDataStringType *>( dictionary.Find( "0018|1030" )->second.GetPointer() )->GetMetaDataObjectValue();

      std::cout << "Acquisition Date: " << acquisition_date << std::endl;
      std::cout << "Acquisition Time: " << acquisition_time << std::endl;
      std::cout << "Protocol Name: " << protocol << std::endl;

      std::cout << "Slice Thickness: " << dynamic_cast<const MetaDataStringType *>( dictionary.Find( "0018|0050" )->second.GetPointer() )->GetMetaDataObjectValue() << std::endl;
      std::cout << "Pixel Spacing: " << dynamic_cast<const MetaDataStringType *>( dictionary.Find( "0028|0030" )->second.GetPointer() )->GetMetaDataObjectValue() << std::endl;

      std::cout << "Image Orientation: " << dynamic_cast<const MetaDataStringType *>( dictionary.Find( "0020|0037" )->second.GetPointer() )->GetMetaDataObjectValue() << std::endl;

      typedef itk::ImageFileWriter< ImageType > WriterType;
      WriterType::Pointer writer = WriterType::New();

      std::string output_filename = acquisition_date+"_"+acquisition_time+protocol;

      //Remove spaces if any in the output_filename
      std::replace(output_filename.begin(), output_filename.end(), ' ', '_');
      std::replace(output_filename.begin(), output_filename.end(), '.', '_');

      writer->SetFileName( output_filename+".nii.gz" );

      writer->SetInput( reader->GetOutput() );

      std::cout  << "Writing the image as " << std::endl << std::endl;
      std::cout  << output_filename << std::endl << std::endl;

      try
        {
        writer->Update();
        }
      catch (itk::ExceptionObject &ex)
        {
        std::cout << ex << std::endl;
        return EXIT_FAILURE;
        }
      }
    catch (itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
      }

    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------------

    return 1;
  } catch (TCLAP::ArgException &e)  // catch any exceptions
  { std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }

}
