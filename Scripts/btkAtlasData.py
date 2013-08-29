#!/usr/bin/python
# -*- coding: utf-8 -*-
#############################################################################
#
#  © Université de Strasbourg - Centre National de la Recherche Scientifique
#
#  Date: 01/12/2011
#  Author(s): Julien Pontabry (pontabry@unistra.fr)
#
#  This software is governed by the CeCILL-B license under French law and
#  abiding by the rules of distribution of free software.  You can  use,
#  modify and/ or redistribute the software under the terms of the CeCILL-B
#  license as circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and,  more generally, to use and operate it in the
#  same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL-B license and that you accept its terms.
#
#############################################################################


# Path to data directory
dataPath = '/path/to/data'

# Path to output directory
outputPath = '/path/to/output'

# Activate or desactivate script (display commands only)
scriptOn = False

# Number of processes used
nbOfProcesses = 1


#############################################################################
#                             Modalities                                    #
#############################################################################

modalities = {}
UseInRegistration = 'UseInRegistration'
UseInRegression   = 'UseInRegression'
ModalityWeight    = 'Weight'
ModalityDataPath  = 'DataPath'
IsTissueMap       = 'IsTissueMap'
TissueLabel       = 'TissueLabel'

## There you can put all the modalities you need, following the example below.

## Modality example
# ModalityImage = 'ModalityID'                                            # Identifier of the modality
# modalities[ModalityImage] = {}                                          # Initialization
# modalities[ModalityImage][UseInRegistration] = True                     # Is the modality used during registration ?
# modalities[ModalityImage][UseInRegression]   = True                     # Is the modality used during regression ?
# modalities[ModalityImage][ModalityWeight]    = 1                        # The weight of the modality (normalized after)
# modalities[ModalityImage][ModalityDataPath]  = dataPath + '/ModalityID' # The path to modality
# modalities[ModalityImage][IsTissueMap]       = False                    # Is the modality a label
# modalities[ModalityImage][TissueLabel]       = 0                        # If the modality is a label, which number in tissue map ?



#############################################################################
#                               Patients                                    #
#############################################################################

## There you can fill the patients information, following the example below.

# patients = [ ('Identifier1', age1), ('Identifier2', age2), ... , ('IdentifierN', ageN) ]



#############################################################################
#                          Template creation                                #
#############################################################################

# Reference patient for template creation
patientReference = 'Identifier1'

# Working directory of template creation
templatePath = outputPath + '/template'

# Number of iterations for each step (ANTS syntax)
registrationSteps = '1000x500x200x150x150x150'

# Gradient step (ANTS syntax)
gradientStep = '0.25'


#############################################################################
#                      Longitudinal Atlas creation                          #
#############################################################################

# Time step for sampling
timeStep = 1.0

# Bandwith parameter
bandwith = 1.0

# Working directory of atlas creation
atlasPath = outputPath + '/atlas'


#############################################################################
#                             External programs                             #
#############################################################################

BtkBinaryDir         = ''
BinarizeLabels       = 'btkBinarizeLabels'
GaussianFilter       = 'btkImageGaussianFilter'
ProbMapNormalization = 'btkProbabilityMapNormalization'
WeightedSum          = 'btkWeightedSumOfImages'
WeightedSumAffine    = 'btkWeightedSumOfAffineTransforms'
BinarizeMaps         = 'btkBinarizeTissueProbabilityMaps'
InverseField         = 'btkInverseDisplacementField'
CropUsingMask        = 'btkCropImageUsingMask'
HistogramMatching    = 'btkHistogramMatching'
ComputeDistAffine    = 'btkComputeDistanceBetweenAffineTransforms'
ComputeDistImage     = 'btkComputeDistanceBetweenImages'

AntsBinaryDir    = ''
ANTS             = 'ANTS'
Warp             = 'WarpImageMultiTransform'
ComposeTransform = 'ComposeMultiTransform'


if len(BtkBinaryDir) > 0:
	BtkBinaryDir += '/'

if len(AntsBinaryDir) > 0:
	AntsBinaryDir += '/'

