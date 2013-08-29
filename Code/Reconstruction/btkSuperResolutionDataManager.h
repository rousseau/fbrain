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
#include "itkTransformFactory.h"

#include "itkAffineTransform.h"
#include "itkRigid3DTransform.h"
#include "itkEuler3DTransform.h"
#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"
#include "itkMatrixOffsetTransformBase.h"

#include "btkIOTransformHelper.h"

class SuperResolutionDataManager
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 3>          itkImage;
  typedef itkImage::Pointer                  itkPointer;  
  typedef itk::ImageFileReader< itkImage >   itkReader;
  typedef itk::ImageFileWriter< itkImage >   itkWriter;
  typedef itk::ImageDuplicator< itkImage >   itkDuplicator;
  typedef itk::Image< unsigned char, 3 >     itkImageMask;
  typedef itk::ImageFileReader< itkImageMask > itkMaskReader;
  typedef itk::ImageMaskSpatialObject< 3 >   itkMask;
  typedef itkMask::Pointer                   itkMaskPointer;

  //typedef itk::Rigid3DTransform<double>     itkAffineDeformation;
  //typedef itk::AffineTransform<double,3>     itkAffineDeformation;
  typedef itk::TransformFileReader           itkTransformReader;
  typedef itk::MatrixOffsetTransformBase<double,3,3> itkTransformType;
  
  itkPointer                  m_inputHRImage;
  itkPointer                  m_currentHRImage;
  itkPointer                  m_outputHRImage;
  itkPointer                  m_maskHRImage; //should be in short or bool to save memory
  std::vector<itkPointer>     m_inputLRImages;
  std::vector<itkPointer>     m_simulatedInputLRImages;
  std::vector<itkMaskPointer>  m_maskLRImages;
  std::vector<itkImage::RegionType> m_regionLRImages;
  std::vector<itkTransformType::Pointer> m_affineTransform;
  std::vector<itkTransformType::Pointer> m_inverseAffineTransform;

  
  
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
  void WriteSimulatedLRImages(std::vector<std::string> & input_file);
  void WriteOutputHRImage(std::string & input_file);
  void WriteOneImage(itkPointer image, std::string & filename);
  
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
  
  itkDuplicator::Pointer duplicator2 = itkDuplicator::New();
  duplicator2->SetInputImage( m_inputHRImage );
  duplicator2->Update();
  m_currentHRImage = duplicator2->GetOutput();
  
  itkDuplicator::Pointer duplicator3 = itkDuplicator::New();
  duplicator3->SetInputImage( m_inputHRImage );
  duplicator3->Update();
  m_maskHRImage = duplicator3->GetOutput();
}

void SuperResolutionDataManager::ReadLRImages(std::vector<std::string> & input_file)
{
  m_inputLRImages.resize(input_file.size());
  m_simulatedInputLRImages.resize(input_file.size());
  
  for(unsigned int i=0;i<input_file.size();i++){
    std::cout<<"Reading Low-Resolution Input Image : "<<input_file[i]<<"\n";
    itkReader::Pointer reader = itkReader::New();
    reader->SetFileName( input_file[i]  );
    reader->Update();
    m_inputLRImages[i] = reader->GetOutput();
    
    //duplicate the input LR image into the output simulated image to keep all header information
    itkDuplicator::Pointer duplicator = itkDuplicator::New();
    duplicator->SetInputImage( m_inputLRImages[i] );
    duplicator->Update();
    m_simulatedInputLRImages[i] = duplicator->GetOutput();
    m_simulatedInputLRImages[i]->FillBuffer(0);
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
  m_affineTransform.resize(m_inputLRImages.size());
  m_inverseAffineTransform.resize(m_inputLRImages.size());

  if(input_file.size() > 0){
  
    if(m_inputLRImages.size() != input_file.size())
      std::cout<<"WARNING : the number of input transforms is different to the number of LR images!!\n";
    
    itk::TransformFactory<itkTransformType>::RegisterTransform();
  
    for(unsigned int i=0;i<input_file.size();i++){
     std::cout<<"Reading Affine Transform : "<<input_file[i]<<"\n";
     /*
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
      */
      m_affineTransform[i] = btk::IOTransformHelper< itkTransformType >::ReadTransform( input_file[i] );
      //m_affineTransform[i] = static_cast<itkAffineDeformation *>( reader->GetTransformList()->front().GetPointer());
      m_inverseAffineTransform[i] = itkTransformType::New();
      m_inverseAffineTransform[i]->SetCenter( m_affineTransform[i]->GetCenter() );   
      m_affineTransform[i]->GetInverse(m_inverseAffineTransform[i]); 
      //std::cout<<"affine transform:"<<m_affineTransform[i]<<"\n---------------\n";
      //--------
      std::cout<<"Transform :"<<m_affineTransform[i]->GetNameOfClass()<<std::endl;
      std::cout<<"Inverse Transform :"<<m_inverseAffineTransform[i]->GetNameOfClass()<<std::endl;

      //--------

    }
  }
  //else set the transforms with identity 
  else
  {
    std::cout<<"Affine transforms are set to identity\n";
    for(unsigned int i=0;i<m_inputLRImages.size();i++){
      m_affineTransform[i] = itkTransformType::New();
      m_inverseAffineTransform[i] = itkTransformType::New();
      std::cout<<"affine transform:"<<m_affineTransform[i]<<"\n";
    }
  
  }  
}


void SuperResolutionDataManager::WriteSimulatedLRImages(std::vector<std::string> & input_file)
{
  for(unsigned int i=0;i<input_file.size();i++){    
    std::string output_file = "simulated_"+input_file[i];
    WriteOneImage(m_simulatedInputLRImages[i], output_file);
  }
}

void SuperResolutionDataManager::WriteOutputHRImage(std::string & input_file)
{
  //Copy current estimate to the output HR image
  itkDuplicator::Pointer duplicator = itkDuplicator::New();
  duplicator->SetInputImage( m_currentHRImage );
  duplicator->Update();
  m_outputHRImage = duplicator->GetOutput();
  
  WriteOneImage(m_outputHRImage, input_file);
}

void SuperResolutionDataManager::WriteOneImage(itkPointer image, std::string & filename)
{
  std::cout<<"Writing : "<<filename<<"\n";
  WriterType::Pointer writer =  WriterType::New();
  writer -> SetFileName( filename );
  writer -> SetInput( image );
  writer -> Update();
}
 
#endif
