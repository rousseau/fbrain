/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 05/07/2012
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

#include "btkDiffusionSequenceHelper.h"


// STL includes
#include "iostream"

// Local includes
#include "btkDiffusionSequenceFileReader.h"
#include "btkDiffusionSequenceFileWriter.h"


namespace btk
{

void DiffusionSequenceHelper::WriteSequence(btk::DiffusionSequence::Pointer sequence, const std::string &fileName)
{
    std::cout << "Writing \"" << fileName << "\"... " << std::flush;

    btk::DiffusionSequenceFileWriter::Pointer writer = btk::DiffusionSequenceFileWriter::New();
    writer->SetFileName(fileName);
    writer->SetInput(sequence);
    writer->Update();

    std::cout << "done." << std::endl;
}

//----------------------------------------------------------------------------------------

void DiffusionSequenceHelper::WriteSequenceArray(std::vector< btk::DiffusionSequence::Pointer > &sequences, std::vector< std::string > &fileNames)
{
    if(sequences.size() == fileNames.size())
    {
        for(int i = 0; i < sequences.size(); i++)
        {
            WriteSequence(sequences[i], fileNames[i]);
        }
    }
    else
    {
        std::string err("vector of images and vector of names have not the same size !");
        throw err;

    }
}

//----------------------------------------------------------------------------------------

btk::DiffusionSequence::Pointer DiffusionSequenceHelper::ReadSequence(const std::string &fileName)
{
    btk::DiffusionSequenceFileReader::Pointer reader = btk::DiffusionSequenceFileReader::New();
    reader->SetFileName(fileName);
    reader->Update();
    std::cout << "Reading image \"" << fileName << "\"... done." << std::endl;
    return reader->GetOutput();
}

//----------------------------------------------------------------------------------------

std::vector< btk::DiffusionSequence::Pointer > &DiffusionSequenceHelper::ReadSequenceArray(std::vector<std::string> &fileNames)
{
    std::vector< btk::DiffusionSequence::Pointer > *ptrSequences = new std::vector< btk::DiffusionSequence::Pointer >;
    std::vector< btk::DiffusionSequence::Pointer > &sequences = *ptrSequences;
    sequences.resize(fileNames.size());

    for(int i = 0; i < fileNames.size(); i++)
    {
        sequences[i] = ReadSequence(fileNames[i]);
    }

    return sequences;
}

//----------------------------------------------------------------------------------------

btk::DiffusionSequence::Pointer DiffusionSequenceHelper::CreateNewSequenceFromPhysicalSpaceOf(btk::DiffusionSequence::Pointer sequence)
{
//    // Define header
//    btk::DiffusionSequence::SizeType           sequenceSize = sequence->GetLargestPossibleRegion().GetSize();
//    btk::DiffusionSequence::PointType        sequenceOrigin = sequence->GetOrigin();
//    btk::DiffusionSequence::SpacingType     sequenceSpacing = sequence->GetSpacing();
//    btk::DiffusionSequence::DirectionType sequenceDirection = sequence->GetDirection();

//    typename TSequenceOutput::SizeType newImageSize;
//    newImageSize[0] = sequenceSize[0]; newImageSize[1] = sequenceSize[1]; newImageSize[2] = sequenceSize[2];

//    typename TSequenceOutput::PointType newImageOrigin;
//    newImageOrigin[0] = sequenceOrigin[0]; newImageOrigin[1] = sequenceOrigin[1]; newImageOrigin[2] = sequenceOrigin[2];

//    typename TSequenceOutput::SpacingType newImageSpacing;
//    newImageSpacing[0] = sequenceSpacing[0]; newImageSpacing[1] = sequenceSpacing[1]; newImageSpacing[2] = sequenceSpacing[2];

//    typename TSequenceOutput::DirectionType newImageDirection;
//    newImageDirection(0,0) = sequenceDirection(0,0); newImageDirection(0,1) = sequenceDirection(0,1); newImageDirection(0,2) = sequenceDirection(0,2);
//    newImageDirection(1,0) = sequenceDirection(1,0); newImageDirection(1,1) = sequenceDirection(1,1); newImageDirection(1,2) = sequenceDirection(1,2);
//    newImageDirection(2,0) = sequenceDirection(2,0); newImageDirection(2,1) = sequenceDirection(2,1); newImageDirection(2,2) = sequenceDirection(2,2);

//    // Create new image
//    typename TSequenceOutput::Pointer newImage = TSequenceOutput::New();

//    // If the new image is a diffusion sequence, additionnal header
//    if(typeid(TSequenceOutput) == typeid(btk::DiffusionSequence))
//    {
//        newImageSize[3] = sequenceSize[3];
//        newImageOrigin[3] = sequenceOrigin[3];
//        newImageSpacing[3] = sequenceSpacing[3];

//        newImageDirection(0,3) = sequenceDirection(0,3);
//        newImageDirection(1,3) = sequenceDirection(1,3);
//        newImageDirection(2,3) = sequenceDirection(2,3);
//        newImageDirection(3,3) = sequenceDirection(3,3);
//        newImageDirection(3,0) = sequenceDirection(3,0);
//        newImageDirection(3,1) = sequenceDirection(3,1);
//        newImageDirection(3,2) = sequenceDirection(3,2);

//        newImage->SetGradientTable(sequence->GetGradientTable);
//        newImage->SetBValues(sequence->GetBValues);
//    }

    btk::DiffusionSequence::Pointer newSequence = btk::DiffusionSequence::New();
    newSequence->SetRegions(sequence->GetLargestPossibleRegion());
    newSequence->SetOrigin(sequence->GetOrigin());
    newSequence->SetSpacing(sequence->GetSpacing());
    newSequence->SetDirection(sequence->GetDirection());
    newSequence->SetGradientTable(sequence->GetGradientTable());
    newSequence->SetBValues(sequence->GetBValues());

    newSequence->Allocate();
    newSequence->FillBuffer(0);

    return newSequence;
}

//----------------------------------------------------------------------------------------

std::vector< btk::DiffusionSequence::Pointer > &DiffusionSequenceHelper::CreateNewSequenceFromPhysicalSpaceOf(std::vector< btk::DiffusionSequence::Pointer > &sequences)
{
    std::vector< btk::DiffusionSequence::Pointer > *ptrNewImages = new std::vector< btk::DiffusionSequence::Pointer >;
    std::vector< btk::DiffusionSequence::Pointer > &newImages = *ptrNewImages;

    for(std::vector< btk::DiffusionSequence::Pointer >::iterator it = sequences.begin(); it != sequences.end(); it++)
    {
        newImages.push_back(CreateNewSequenceFromPhysicalSpaceOf(*it));
    }

    return newImages;
}

} // namespace btk
