/*
 * Copyright or © or Copr. Université de Strasbourg - Centre National de la Recherche Scientifique
 *
 * 06 july 2011
 * < pontabry at unistra dot fr >
 *
 * This software is governed by the CeCILL-B license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL-B
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-B license and that you accept its terms.
 */

// TCLAP : Templatized C++ Command Line Parser includes
#include <tclap/CmdLine.h>

// STL includes
#include "string"
#include "iostream"

// VTK includes
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkAppendPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyLine.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"

// Local includes
#include "btkTypes.h"
#include "btkPoint.h"
#include "btkDirection.h"
#include "btkSignal.h"
#include "btkSignalExtractor.h"
#include "btkSHModel.h"
#include "btkSHModelEstimator.h"
#include "btkFileNameTools.h"


using namespace btk;


bool IsInsideMask(Point p, Mask::Pointer mask)
{
    Mask::IndexType index;
    index[0] = std::floor(p.x() + 0.5);
    index[1] = std::floor(p.y() + 0.5);
    index[2] = std::floor(p.z() + 0.5);

    MaskRegion region = mask->GetLargestPossibleRegion();

    bool isIn = false;

    if(index[0] >= 0 && index[0] < (unsigned int)region.GetSize(0) &&
       index[1] >= 0 && index[1] < (unsigned int)region.GetSize(1) &&
       index[2] >= 0 && index[2] < (unsigned int)region.GetSize(2))
        isIn = (mask->GetPixel(index) > 0) ? true : false;

    return isIn;
}

Direction MeanDirection(SHModel *modelFun, Point currentPoint, Direction lastDirection, Real angleThreshold)
{
    std::vector<Direction> maxima = modelFun->getMaxDirectionsAt(currentPoint);

    if(maxima.size() > 1)
    {
        Vector lastVector = lastDirection.toVector();
        std::vector<Direction> mus;
        std::vector<Real> values;

        // we keep only maxima in the solid angle
        for(std::vector<Direction>::iterator it=maxima.begin(); it!=maxima.end(); it++)
        {
            Vector v = it->toVector();

            Real scalProd = v.x()*lastVector.x() + v.y()*lastVector.y() + v.z()*lastVector.z();
            Real normU    = v.x()*v.x() + v.y()*v.y() + v.z()*v.z();
            Real normV    = lastVector.x()*lastVector.x() + lastVector.y()*lastVector.y() + lastVector.z()*lastVector.z();
            Real alpha    = std::acos( scalProd / std::sqrt(normU * normV) );

            if(alpha <= angleThreshold)
            {
                mus.push_back(*it);
                values.push_back(modelFun->odfAt(*it, currentPoint));
            }
        } // for each maxima

        if(mus.size() > 1)
        {
            // Find maximal direction
            Direction maxDir;
            Real maxValue = 0;

            for(unsigned int i=0; i<mus.size(); i++)
            {
                Real tmp = values[i];

                if(maxValue < tmp)
                {
                    maxValue = tmp;
                    maxDir   = mus[i];
                }
            }

            return maxDir;
        }
        else if(mus.size() == 1)
            return mus.front();
        else // mus.size() == 0
        {
            Direction nullDir;
            nullDir.setNull();

            return nullDir;
        }
    }
    else if(maxima.size() == 1)
        return maxima.front();
    else // maxima.size() == 0
    {
        Direction nullDir;
        nullDir.setNull();

        return nullDir;
    }
}

Real ComputeGFA(SHModel *model, Signal* signal, Point p)
{
    Matrix S = signal->signalAt(p);
    Matrix M = model->signalAt(p);

    Real psi_average = 0;
    for(unsigned int i=0; i<S.Rows(); i++)
        psi_average += S(i,0);
    psi_average /= S.Rows();

    Real psi_std = 0;
    for(unsigned int i=0; i<S.Rows(); i++)
        psi_std += (S(i,0)-psi_average) * (S(i,0)-psi_average);
    psi_std /= S.Rows();
    psi_std = std::sqrt(psi_std);

    Real psi_rms = 0;
    for(unsigned int i=0; i<S.Rows(); i++)
        psi_rms += S(i,0) * S(i,0);
    psi_rms /= S.Rows();
    psi_rms = std::sqrt(psi_rms);

    Real GFA = psi_std / psi_rms;

    return GFA;
}


int main(int argc, char *argv[])
{
    // Command line variables
    std::string dwiFileName;
    std::string vecFileName;
    std::string maskFileName;
    std::string labelFilename;
    std::vector<int> select_labels;

    std::string outFibersFileName;

    bool verboseMode;
    bool quietMode;
    bool saveTmpFiles;
    bool lps;

    unsigned int modelOrder;
    Real lambda;
    Real stepSize;
    Real angleThreshold;
    Float seedSpacing;

    std::string inRadix;


    try
    {
        //
        // Program's command line parser definition
        //

            // Defines command line parser
            TCLAP::CmdLine cmd("BTK Tractography", ' ', "0.2");

            // Defines arguments
            TCLAP::ValueArg<std::string>   dwiArg("d", "dwi", "Dwi sequence", true, "", "string", cmd);
            TCLAP::ValueArg<std::string>  maskArg("m", "mask", "White matter mask", true, "", "string", cmd);
            TCLAP::ValueArg<std::string> labelArg("l", "label", "Label volume of seeds", true, "", "string", cmd);
            TCLAP::MultiArg<int>        labelsArg("", "labels", "Label select for processing", false, "int", cmd);

            TCLAP::ValueArg<std::string> outFibersArg("", "fibers", "Output fibers filename (default \"fibers-label\")", false, "fibers", "string", cmd);

            TCLAP::SwitchArg verboseSwitchArg("", "verbose", "Display more informations on standard output", cmd, false);
            TCLAP::SwitchArg quietSwitchArg("", "quiet", "Display no information on either standard and error outputs", cmd, false);
            TCLAP::SwitchArg saveTmpSwitchArg("", "save_temporary_files", "Save diffusion signal, model coefficients, variance and spherical coordinates of gradient directions in files", cmd, false);
            TCLAP::SwitchArg lpsSwitchArg("", "lps", "Word coordinates expressed in LPS (Left-Posterior-Superior). By default RAS (Right-Anterior-Superior) is used.", cmd, false);

            TCLAP::ValueArg<unsigned int> orderArg("", "model_order", "Order of the model (i.e. of spherical harmonics) (default 4)", false, 4, "unsigned int", cmd);
            TCLAP::ValueArg<Real>    lambdArg("", "model_regularization", "Regularization coefficient of the model (default 0.006)", false, 0.006, "Real", cmd);
            TCLAP::ValueArg<Real>    stepSizeArg("", "step_size", "Step size of particles displacement (default 0.5)", false, 0.5, "Real", cmd);
            TCLAP::ValueArg<Real>    angleThreshArg("", "angular_threshold", "Angular threshold between successive displacement vector of a particle's trajectory (default PI/3)", false, M_PI/2., "Real", cmd);
            TCLAP::ValueArg<Real>    seedSpacingArg("", "seed_spacing", "Spacing in mm between seeds (default 1)", false, 1.0, "float", cmd);

            // Parsing arguments
            cmd.parse(argc, argv);

            // Get back arguments' values
            dwiFileName    = dwiArg.getValue();
            maskFileName   = maskArg.getValue();
            labelFilename  = labelArg.getValue();
            select_labels  = labelsArg.getValue();

            outFibersFileName = outFibersArg.getValue();

            verboseMode  = verboseSwitchArg.getValue();
            quietMode    = quietSwitchArg.getValue();
            saveTmpFiles = saveTmpSwitchArg.getValue();
            lps          = lpsSwitchArg.getValue();

            modelOrder     = orderArg.getValue();
            lambda         = lambdArg.getValue();
            stepSize       = stepSizeArg.getValue();
            angleThreshold = angleThreshArg.getValue();
            seedSpacing    = seedSpacingArg.getValue();


            inRadix     = GetRadixOf(dwiFileName);
            vecFileName = inRadix + ".bvec";
    }
    catch(TCLAP::ArgException &e)
    {
        std::cout << "TCLAP error: " << e.error() << " for argument " << e.argId() << std::endl;
        std::exit(EXIT_FAILURE);
    }


    char displayMode = quietMode ? 0 : (verboseMode ? 2 : 1);


    //
    // Diffusion signal extraction
    //

        SignalExtractor *extractor = new SignalExtractor(vecFileName, dwiFileName, maskFileName, displayMode);
        extractor->extract();

        std::vector<Direction> *directions = extractor->GetDirections();
        Sequence::Pointer       signal     = extractor->GetSignal();
        std::vector<Real>      *sigmas     = extractor->GetSigmas();
        Mask::Pointer           mask       = extractor->GetMask();

        if(saveTmpFiles)
            extractor->save();

        delete extractor;


    //
    // Model estimation
    //

        SHModelEstimator *estimator = new SHModelEstimator(signal, directions, mask, modelOrder, lambda, displayMode);
        estimator->estimate();

        Sequence::Pointer model = estimator->GetModel();

        if(saveTmpFiles)
            estimator->save();

        delete estimator;




    //
    // Tractography (streamline propagation)
    //


        // Get spherical harmonics model
        SHModel *modelFun = new SHModel(model, directions, displayMode);

        // Get signal
        Signal *signalFun = new Signal(signal, sigmas, directions, displayMode);

        // Read label image if any and verify sizes
        LabelMapReader::Pointer labelReader = LabelMapReader::New();
        labelReader->SetFileName(labelFilename);
        labelReader->Update();
        LabelMap::Pointer labelVolume = labelReader->GetOutput();


        // Initialize output structures
        std::vector<int> labels;
        std::vector<vtkSmartPointer<vtkAppendPolyData> > fibers;


        // Resample label map
        LabelResampler::Pointer labelResampler = LabelResampler::New();
        labelResampler->SetInput(labelVolume);

        LabelInterpolator::Pointer labelInterpolator = LabelInterpolator::New();
        labelResampler->SetInterpolator(labelInterpolator);

        LabelMap::SpacingType spacing = labelVolume->GetSpacing();
        LabelMap::SizeType size = labelVolume->GetLargestPossibleRegion().GetSize();
        size[0] *= spacing[0]/seedSpacing; size[1] *= spacing[1]/seedSpacing; size[2] *= spacing[2]/seedSpacing;
        labelResampler->SetSize(size);

        LabelMap::DirectionType direction = labelVolume->GetDirection();
        labelResampler->SetOutputDirection(direction);

        LabelMap::PointType origin = labelVolume->GetOrigin();
        labelResampler->SetOutputOrigin(origin);

        spacing[0] = seedSpacing; spacing[1] = seedSpacing; spacing[2] = seedSpacing;
        labelResampler->SetOutputSpacing(spacing);

        labelResampler->Update();
        LabelMap::Pointer resampledLabelVolume = labelResampler->GetOutput();

        // Apply filter on each labeled voxels
        LabelMapIterator labelIt(resampledLabelVolume, resampledLabelVolume->GetLargestPossibleRegion());

        unsigned int numOfSeeds = 0;

        for(labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt)
        {
            if(labelIt.Get() > 0 &&
               (select_labels.size() == 0 ||
                *std::find(select_labels.begin(),select_labels.end(),labelIt.Get()) == labelIt.Get()) )
                numOfSeeds++;
        }

        unsigned int numOfSeedsDone = 0;

        std::cout << std::fixed;
        std::cout << std::setprecision(2);
        Display1(displayMode, std::cout << "Running tractography..." << std::endl);
        for(labelIt.GoToBegin(); !labelIt.IsAtEnd(); ++labelIt)
        {
            int label = (int)labelIt.Get();

            if(label > 0 &&
               (select_labels.size() == 0 ||
                *std::find(select_labels.begin(),select_labels.end(),label) == label) )
            {
                LabelMap::IndexType index = labelIt.GetIndex();

                itk::Point<Real,3> worldPoint;
                resampledLabelVolume->TransformIndexToPhysicalPoint(index, worldPoint);

                Mask::IndexType maskIndex;
                mask->TransformPhysicalPointToIndex(worldPoint, maskIndex);

                // If the seed is not in the mask, there is no need to continue this one
                if(mask->GetPixel(maskIndex) != 0)
                {
                    Point begin(worldPoint[0], worldPoint[1], worldPoint[2]);

                    // Let's start filtering
                    Display1(displayMode, std::cout << "\tProgress: " << (Real)numOfSeedsDone / (Real)numOfSeeds * 100 << "%\r" << std::flush);
                    Display2(displayMode, std::cout << "\tSeed's world coordinates: (" << worldPoint[0] << "," << worldPoint[1] << "," << worldPoint[2] << ")" << std::endl);
                    Display2(displayMode, std::cout << "\tSeed's image coordinates: (" << maskIndex[0] << "," << maskIndex[1] << "," << maskIndex[2] << ")" << std::endl);


                    //
                    // Streamline propagation
                    //

                    // Entrées : modelFun, mask, begin, stepSize, angleThreshold
                    // Sortie : fiber (vtkPolyData)

                    // Principe : propagation streamline (suivi direction max de chaque ODF -à chaque point de l'image) en avançant dans la direction d'un pas égal à stepSize
                    // Condition d'arrêt : Pas de direction max dans l'angle demandé ou point hors du masque

                    // Algorithme
                    // chercher direction max du point courant et sa symétrique
                    // pour chaque direction, exécuter les instructions suivantes :
                    // initialiser dernière direction à la direction de départ
                    // initialiser le point courant à begin + déplacement dans la dernière direction
                    // initialiser continuer à vrai
                    // tant qu'on peut continuer
                    // faire
                    // |    si point courant pas dans le masque
                    // |    |   alors stop
                    // |    |   sinon si recherche direction max dans seuil ne donne rien
                    // |    |   |   alors stop
                    // |    |   |   sinon propager et mettre à jour point courant et dernière direction
                    // |    |   fsi
                    // |    fsi
                    // ffaire

                    // Variables
                    Direction lastDirection;
                    Point currentPoint;
                    std::vector<Point> part1, part2;
                    std::vector<Real> gfa1, gfa2;
                    bool goOn;

                    itk::Point<Real,3> worldPoint;
                    worldPoint[0] = begin.x(); worldPoint[1] = begin.y(); worldPoint[2] = begin.z();

                    itk::ContinuousIndex<Real,3> continuousIndex;
                    if(!mask->TransformPhysicalPointToContinuousIndex(worldPoint, continuousIndex))
                    {
                        std::cout << "Error: continuous index not in image !" << std::endl;
                        std::cerr << "(" << worldPoint[0] << "," << worldPoint[1] << "," << worldPoint[2] << ") --> (" << continuousIndex[0] << "," << continuousIndex[1] << "," << continuousIndex[2] << ")" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    Point start(continuousIndex[0], continuousIndex[1], continuousIndex[2]);

                    // Get maximal and symetric directions at start
                    Direction maxDirection = modelFun->getMaxDirectionAt(start);
                    Direction symDirection(
                            M_PI - maxDirection.theta(),
                            (maxDirection.phi() < M_PI) ? maxDirection.phi() + M_PI : maxDirection.phi() - M_PI
                    );

                    // Treat maximal direction
                    lastDirection = maxDirection;
                    currentPoint = start + lastDirection.toVector()*stepSize;
                    part1.push_back(currentPoint);
                    gfa1.push_back(ComputeGFA(modelFun,signalFun,currentPoint));
                    goOn = true;

                    while(goOn)
                    {
                        if(!IsInsideMask(currentPoint,mask))
                        {
                            goOn = false;
                        }
                        else // IsInsideMask(currentPoint,mask)
                        {
                            Direction mean = MeanDirection(modelFun, currentPoint, lastDirection, angleThreshold);

                            if(mean.isNull())
                            {
                                goOn = false;
                            }
                            else // !mean.isNull()
                            {
                                lastDirection = mean;
                                currentPoint = currentPoint + mean.toVector()*stepSize;
                                part1.push_back(currentPoint);
                                gfa1.push_back(ComputeGFA(modelFun,signalFun,currentPoint));
                            }
                        }
                    }

                    // Treat symetric direction
                    lastDirection = symDirection;
                    currentPoint = start + lastDirection.toVector()*stepSize;
                    part2.push_back(currentPoint);
                    gfa2.push_back(ComputeGFA(modelFun,signalFun,currentPoint));
                    goOn = true;

                    while(goOn)
                    {
                        if(!IsInsideMask(currentPoint,mask))
                        {
                            goOn = false;
                        }
                        else // IsInsideMask(currentPoint,mask)
                        {
                            Direction mean = MeanDirection(modelFun, currentPoint, lastDirection, angleThreshold);

                            if(mean.isNull())
                            {
                                goOn = false;
                            }
                            else // !mean.isNull()
                            {
                                lastDirection = mean;
                                currentPoint = currentPoint + mean.toVector()*stepSize;
                                part2.push_back(currentPoint);
                                gfa2.push_back(ComputeGFA(modelFun,signalFun,currentPoint));
                            }
                        }
                    }


                    //
                    // Compute polydata
                    //

                    // VTK structures
                    vtkSmartPointer<vtkPolyData> fiber = vtkSmartPointer<vtkPolyData>::New();

                    vtkSmartPointer<vtkPoints>   points = vtkSmartPointer<vtkPoints>::New();
                    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
                    vtkSmartPointer<vtkPolyLine>   line = vtkSmartPointer<vtkPolyLine>::New();
                    vtkSmartPointer<vtkDoubleArray> gfa = vtkSmartPointer<vtkDoubleArray>::New();

                    // Data
                    vtkIdType pid;
                    unsigned int id = 0;
                    double valueGFA = 0;

                    ImageContinuousIndex ci;
                    Image::PointType wp;


                    // Build fiber with the fiber in world coordinates
                    gfa->SetNumberOfComponents(1);
                    line->GetPointIds()->SetNumberOfIds(part1.size() + part2.size());

                    // First point
                    Point p = part1[part1.size()-1];
                    ci[0] = p.x(); ci[1] = p.y(); ci[2] = p.z();
                    mask->TransformContinuousIndexToPhysicalPoint(ci,wp);

                    if(lps)
                        pid = points->InsertNextPoint(wp[0], wp[1], wp[2]);
                    else // ras
                        pid = points->InsertNextPoint(-wp[0], -wp[1], wp[2]);

                    valueGFA = gfa1[gfa1.size()-1];
                    gfa->InsertNextTuple(&valueGFA);
                    line->GetPointIds()->SetId(id++,pid);

                    // First part
                    for(int k=part1.size()-2; k>=0; k--)
                    {
                        Point p = part1[k];

                        Image::PointType wp;
                        ImageContinuousIndex cip;

                        cip[0] = p.x(); cip[1] = p.y(); cip[2] = p.z();
                        mask->TransformContinuousIndexToPhysicalPoint(cip, wp);

                        if(lps)
                            pid = points->InsertNextPoint(wp[0], wp[1], wp[2]);
                        else // ras
                            pid = points->InsertNextPoint(-wp[0], -wp[1], wp[2]);

                        valueGFA = gfa1[k];
                        gfa->InsertNextTuple(&valueGFA);
                        line->GetPointIds()->SetId(id++,pid);
                    }


                    // Last part
                    for(unsigned int k=0; k<part2.size(); k++)
                    {
                        Point p = part2[k];

                        Image::PointType wp;
                        ImageContinuousIndex cip;

                        cip[0] = p.x(); cip[1] = p.y(); cip[2] = p.z();
                        mask->TransformContinuousIndexToPhysicalPoint(cip, wp);

                        if(lps)
                            pid = points->InsertNextPoint(wp[0], wp[1], wp[2]);
                        else // ras
                            pid = points->InsertNextPoint(-wp[0], -wp[1], wp[2]);

                        valueGFA = gfa2[k];
                        gfa->InsertNextTuple(&valueGFA);
                        line->GetPointIds()->SetId(id++,pid);
                    }

                    lines->InsertNextCell(line);
                    fiber->SetPoints(points);
                    fiber->SetLines(lines);
                    fiber->GetPointData()->SetScalars(gfa);


                    // Save data

                    for(int l=labels.size(); l<label; l++)
                        labels.push_back(-1);

                    if(labels[label-1] == -1)
                    {
                        fibers.push_back(vtkSmartPointer<vtkAppendPolyData>::New());
                        labels[label-1] = fibers.size()-1;
                    }

                    fibers[labels[label-1]]->AddInput(fiber);

                    Display2(displayMode, std::cout << "done." << std::endl);

                    numOfSeedsDone++;
                }
            }
        } // for each labeled voxels
        Display1(displayMode, std::cout << "\tProgress: 100.00%" << std::endl << "done." << std::endl);

        delete signalFun;
        delete modelFun;
        delete directions;
        delete sigmas;


    //
    // Writing files
    //

        Display1(displayMode, std::cout << "Writing results..." << std::endl);

        // Write fiber polydata
        Display2(displayMode, std::cout << "\tWriting fibers..." << std::flush);
        for(unsigned int label = 0; label < labels.size(); label++)
        {
            if(labels[label] != -1)
            {
                unsigned int nbOfInputs = fibers[labels[label]]->GetNumberOfInputPorts();
                bool saveFibers = true;
                unsigned int i = 0;

                do
                {
                    if(fibers[labels[label]]->GetNumberOfInputConnections(i++) < 1)
                        saveFibers = false;
                } while(saveFibers && i<nbOfInputs);

                if(saveFibers)
                {
                    std::stringstream filename;
                    filename << outFibersFileName << "-" << label+1 << ".vtk";

                    fibers[labels[label]]->Update();
                    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
                    writer->SetInput(fibers[labels[label]]->GetOutput());
                    writer->SetFileName(filename.str().c_str());
                    writer->Write();
                }
            }
        }
        Display2(displayMode, std::cout << "done." << std::endl);

    Display1(displayMode, std::cout << "done." << std::endl);


    return EXIT_SUCCESS;
}
