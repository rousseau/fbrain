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

# T2 modality
T2Image = 'T2'
modalities[T2Image] = {}
modalities[T2Image][UseInRegistration] = True
modalities[T2Image][UseInRegression]   = True
modalities[T2Image][ModalityWeight]    = 1
modalities[T2Image][ModalityDataPath]  = dataPath + '/T2'
modalities[T2Image][IsTissueMap]       = False
modalities[T2Image][TissueLabel]       = 0

# GM modality
GMImage = 'GM'
modalities[GMImage] = {}
modalities[GMImage][UseInRegistration] = True
modalities[GMImage][UseInRegression]   = True
modalities[GMImage][ModalityWeight]    = 1
modalities[GMImage][ModalityDataPath]  = dataPath + '/GM'
modalities[GMImage][IsTissueMap]       = True
modalities[GMImage][TissueLabel]       = 1

# WM modality
WMImage = 'WM'
modalities[WMImage] = {}
modalities[WMImage][UseInRegistration] = True
modalities[WMImage][UseInRegression]   = True
modalities[WMImage][ModalityWeight]    = 1
modalities[WMImage][ModalityDataPath]  = dataPath + '/WM'
modalities[WMImage][IsTissueMap]       = True
modalities[WMImage][TissueLabel]       = 2

# Cervelet modality
CerveletImage = 'Cervelet'
modalities[CerveletImage] = {}
modalities[CerveletImage][UseInRegistration] = True
modalities[CerveletImage][UseInRegression]   = True
modalities[CerveletImage][ModalityWeight]    = 1
modalities[CerveletImage][ModalityDataPath]  = dataPath + '/Cervelet'
modalities[CerveletImage][IsTissueMap]       = True
modalities[CerveletImage][TissueLabel]       = 3

# Brainstem modality
BrainstemImage = 'Brainstem'
modalities[BrainstemImage] = {}
modalities[BrainstemImage][UseInRegistration] = True
modalities[BrainstemImage][UseInRegression]   = True
modalities[BrainstemImage][ModalityWeight]    = 1
modalities[BrainstemImage][ModalityDataPath]  = dataPath + '/Brainstem'
modalities[BrainstemImage][IsTissueMap]       = True
modalities[BrainstemImage][TissueLabel]       = 4

# CSF modality
CSFImage = 'CSF'
modalities[CSFImage] = {}
modalities[CSFImage][UseInRegistration] = False
modalities[CSFImage][UseInRegression]   = True
modalities[CSFImage][ModalityWeight]    = 0
modalities[CSFImage][ModalityDataPath]  = dataPath + '/CSF'
modalities[CSFImage][IsTissueMap]       = True
modalities[CSFImage][TissueLabel]       = 5

# Other modality
OtherImage = 'Other'
modalities[OtherImage] = {}
modalities[OtherImage][UseInRegistration] = False
modalities[OtherImage][UseInRegression]   = True
modalities[OtherImage][ModalityWeight]    = 0
modalities[OtherImage][ModalityDataPath]  = dataPath + '/Other'
modalities[OtherImage][IsTissueMap]       = True
modalities[OtherImage][TissueLabel]       = 0


#############################################################################
#                               Patients                                    #
#############################################################################

# patients = [ ('Identifier1', age1), ('Identifier2', age2), ... , ('IdentifierN', ageN) ]
patients = [ ('ARS_Hu', 33), ('AYD_Na', 28), ('BAL_In', 32.5), ('BER_Sl', 30), ('BLA_El', 32), ('CHA_Sa', 28), ('DAH_Au', 30), ('DER_An', 26), ('ELO_Ha', 28), ('ESC_Mi', 34), ('FRE_St', 28), ('HER_Au', 32), ('HIE_Au', 28), ('KOG_Fa', 29), ('KRA_Na', 32), ('LIP_La', 33), ('MAG_Ai', 30), ('NEF_Dr', 32), ('RYC_Ca', 27), ('STR_Re', 31), ('TAB_Ar_01', 27), ('TAB_Ar_02', 30), ('TRO_Sa', 32) ]



#############################################################################
#                          Template creation                                #
#############################################################################

# Reference patient for template creation
patientReference = 'KOG_Fa'

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

