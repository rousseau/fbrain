/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 30/11/2011
  Author(s): François Rousseau (rousseau@unistra.fr)

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

#ifndef __btkSuperResolutionDataManager_h
#define __btkSuperResolutionDataManager_h

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageDuplicator.h"

#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"


class SuperResolutionDataManager
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 3>          itkImage;
  typedef itkImage::Pointer                 itkPointer;  
  typedef itk::ImageFileReader< itkImage >  itkReader;
  typedef itk::ImageFileWriter< itkImage >  itkWriter;
  typedef itk::ImageDuplicator< itkImage >  itkDuplicator;
  typedef itk::Image< unsigned char, 3 >     itkImageMask;
  typedef itk::ImageFileReader< itkImageMask > itkMaskReader;
  typedef itk::ImageMaskSpatialObject< 3 >   itkMask;
  typedef itkMask::Pointer                   itkMaskPointer;

  typedef itk::AffineTransform<double,3>     itkAffineDeformation;
  typedef itk::TransformFileReader           itkTransformReader;
  
  itkPointer                  m_inputHRImage;
  itkPointer                  m_outputHRImage;
  itkPointer                  m_maskHRImage; //not yet used
  std::vector<itkPointer>     m_inputLRImages;
  std::vector<itkMaskPointer>  m_maskLRImages;
  std::vector<itkImage::RegionType> m_regionLRImages;
  std::vector<itkAffineDeformation::Pointer> m_affineTransform;

  
  
  //Image information (size, spacing etc.)
  itkImage::SpacingType m_spacing;
  itkImage::SizeType    m_size;
  itkImage::RegionType  m_region;
  
  typedef itk::ImageFileReader< itkImage >   ReaderType;
  typedef itk::ImageFileWriter< itkImage >   WriterType;
  
  typedef itk::TransformFileReader     TransformReaderType;
  typedef TransformReaderType::TransformListType * TransformListType;

  
  void ReadHRImage(std::string input_file);
  void ReadLRImages(std::vector<std::string> & input_file);
  void ReadMaskLRImages(std::vector<std::string> & input_file);
  void ReadAffineTransform(std::vector<std::string> & input_file);
  
};

void SuperResolutionDataManager::ReadHRImage(std::string input_file)
{
  std::cout<<"All images are casted into float type.\n";
  std::cout<<"Reading the initial estimate of the high resolution image : "<<input_file<<std::endl;
  itkReader::Pointer reader = itkReader::New();
  reader->SetFileName( input_file );
  reader->Update();
  m_inputHRImage = reader->GetOutput();
  
  //compute characteristics of the input image
  m_region  = m_inputHRImage->GetLargestPossibleRegion();
  m_size    = m_region.GetSize();
  m_spacing = m_inputHRImage->GetSpacing();
  
  //duplicate the input image into the output images to keep all header information
  itkDuplicator::Pointer duplicator = itkDuplicator::New();
  duplicator->SetInputImage( m_inputHRImage );
  duplicator->Update();
  m_outputHRImage = duplicator->GetOutput();
  m_outputHRImage->FillBuffer(0); 
}

void SuperResolutionDataManager::ReadLRImages(std::vector<std::string> & input_file)
{
  m_inputLRImages.resize(input_file.size());
  
  for(unsigned int i=0;i<input_file.size();i++){
    std::cout<<"Reading Low-Resolution Input Image : "<<input_file[i]<<"\n";
    itkReader::Pointer reader = itkReader::New();
    reader->SetFileName( input_file[i]  );
    reader->Update();
    m_inputLRImages[i] = reader->GetOutput();
  }
}

void SuperResolutionDataManager::ReadMaskLRImages(std::vector<std::string> & input_file)
{
  
  m_maskLRImages.resize(input_file.size());
  m_regionLRImages.resize(input_file.size());
  
  for(unsigned int i=0;i<input_file.size();i++){
    std::cout<<"Reading Low-Resolution Input Mask : "<<input_file[i]<<"\n";
    itkMaskReader::Pointer reader = itkMaskReader::New();
    reader->SetFileName(input_file[i]);
    try
    {
      reader->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "Error while reading the mask file !" << std::endl;
      std::cerr << err << std::endl;
      std::cerr << "[FAILED]" << std::endl;
    }
    //m_maskLRImages[i]->SetImage(reader->GetOutput()); //buggy!!!!!!!!!!!!!!!!!!!!!!!!!!    
    //m_regionLRImages[i] = m_maskLRImages[i]->GetAxisAlignedBoundingBoxRegion();
  }
}

void SuperResolutionDataManager::ReadAffineTransform(std::vector<std::string> & input_file)
{
  m_affineTransform.resize(input_file.size());
  
  for(unsigned int i=0;i<input_file.size();i++){
    std::cout<<"Reading Affine Transform : "<<input_file[i]<<"\n";
    itkTransformReader::Pointer reader = itkTransformReader::New();
    reader->SetFileName(input_file[i]);
    try
    {
      reader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Error while reading the transform file" << std::endl;
      std::cerr << excp << std::endl;
      std::cerr << "[FAILED]" << std::endl;
    } 
    TransformListType transforms = reader->GetTransformList();
    TransformReaderType::TransformListType::const_iterator titr = transforms->begin();

    m_affineTransform[i] = dynamic_cast<itkAffineDeformation *>( titr->GetPointer() ) ; 
  }
}
#endif
