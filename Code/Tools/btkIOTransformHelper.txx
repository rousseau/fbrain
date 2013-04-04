/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 17/04/2012
  Author(s): Marc Schweitzer (marc.schweitzer(at)unistra.fr)

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

#ifndef __BTK_IOTRANSFORMHELPER_TXX__
#define __BTK_IOTRANSFORMHELPER_TXX__

#include "btkIOTransformHelper.h"

#include "itkTransformFactory.h"
#include "iostream"

namespace btk
{
//---------------------------------------------------------------------------------
template <class TTransform >
typename IOTransformHelper< TTransform >::TransformPointerType
IOTransformHelper< TTransform >::
ReadTransform(std::string _fileName)
{
    itk::TransformFactory<TTransform>::RegisterTransform();
    itkTReader::Pointer reader = itkTReader::New();
    typedef itkTReader::TransformListType * TransformListType;

    reader->SetFileName(_fileName);
    try
    {
        std::cout << "Reading " << _fileName.c_str() << " ... " ; std::cout.flush();
        reader->Update();

    }
    catch(itk::ExceptionObject & excpt)
    {
        std::cerr << "Error while reading transform" << std::endl;
        std::cerr << excpt << std::endl;
        std::cout << "[FAILED]" << std::endl;
        throw excpt;
    }
    TransformListType transforms = reader->GetTransformList();
    //    itkTReader::TransformListType::const_iterator titr = transforms->begin();
    //    TransformPointerType trans = dynamic_cast< TTransform * >( titr->GetPointer() );
    TransformPointerType trans = static_cast< TTransform * >( transforms->front().GetPointer() );
    std::cout << " done! " << std::endl;

    return trans;



}
//---------------------------------------------------------------------------------
template <class TTransform>
std::vector< typename IOTransformHelper< TTransform >::TransformPointerType >&
IOTransformHelper< TTransform >::
ReadTransform(std::vector< std::string >& _fileNames)
{
    itk::TransformFactory<TTransform>::RegisterTransform();

    std::vector< TransformPointerType > *ptrTransforms = new std::vector< TransformPointerType >;
    std::vector< TransformPointerType > &transforms = *ptrTransforms;
    transforms.resize(_fileNames.size());

    for(int i = 0; i< _fileNames.size(); i++)
    {
        itkTReader::Pointer reader = itkTReader::New();
        typedef itkTReader::TransformListType * TransformListType;

        reader->SetFileName(_fileNames[i]);
        try
        {
            std::cout << "Reading " << _fileNames[i].c_str() << " ... " ; std::cout.flush();
            reader->Update();

        }
        catch(itk::ExceptionObject & excpt)
        {
            std::cerr << "Error while reading transform" << std::endl;
            std::cerr << excpt << std::endl;
            std::cout << "[FAILED]" << std::endl;
            throw excpt;
        }
        TransformListType tlist = reader->GetTransformList();
        //        itkTReader::TransformListType::const_iterator titr = tlist->begin();
        //transforms[i] = TransformType::New();
        //        TransformPointerType trans = dynamic_cast< TTransform * >( titr->GetPointer() );
        TransformPointerType trans = static_cast< TTransform * >( tlist->front().GetPointer() );
        transforms[i] = trans;
        std::cout << " done! " << std::endl;
    }

    return transforms;
}


//---------------------------------------------------------------------------------
template <class TTransform>
void IOTransformHelper< TTransform >::WriteTransform(TransformPointerType _transform, const std::string _fileName)
{

    itkTWriter::Pointer writer = itkTWriter::New();
    writer->SetInput(_transform);
    writer->SetFileName(_fileName.c_str());
    try
    {
        std::cout << "Writing " << _fileName.c_str() << " ... " ; std::cout.flush();
        writer->Update();
        std::cout << " done! " << std::endl;
    }
    catch(itk::ExceptionObject & excpt)
    {
        std::cerr << "Error while saving transform" << std::endl;
        std::cerr << excpt << std::endl;
        std::cout << "[FAILED]" << std::endl;
        throw excpt;
    }


}

//---------------------------------------------------------------------------------
template <class TTransform>
void IOTransformHelper< TTransform >::WriteTransform(std::vector< TransformPointerType > &_transforms, std::vector< std::string > &_fileNames)
{
    if(_transforms.size() > 0 && _transforms.size() == _fileNames.size())
    {
        for(int i = 0; i< _transforms.size(); i++)
        {


            WriteTransform(_transforms[i],_fileNames[i]);
            //            itkTWriter::Pointer writer = itkTWriter::New();
            //            writer->SetInput(_transforms[i]);
            //            writer->SetFileName(_fileNames[i].c_str());
            //            try
            //            {
            //                std::cout << "Writing " << _fileNames[i].c_str() << " ... " ; std::cout.flush();
            //                writer->Update();
            //                std::cout << " done! " << std::endl;
            //            }
            //            catch(itk::ExceptionObject & excpt)
            //            {
            //                std::cerr << "Error while saving transform" << std::endl;
            //                std::cerr << excpt << std::endl;
            //                std::cout << "[FAILED]" << std::endl;
            //                throw excpt;
            //            }
        }
    }
    else
    {
        itk::ExceptionObject  excpt;
        excpt.SetDescription("Transforms are empty, or does not have the same size as fileNames !");
        throw excpt;
    }
}

//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
}
#endif
