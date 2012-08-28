########################################################
#        Name: data.py                                 #      
# Description: Configuration script for Atlas creation #
#      Author: Julien Pontabry                         #
#        Date: December 2011                           #
########################################################


# Path to data directory
dataPath = '/home/miv/pontabry/Recherches/Atlas/Atlas-23F/data'

# Path to working directory
workPath = '/home/miv/pontabry/Recherches/Atlas/Atlas-23F'

# Activate or desactivate script (display commands only)
scriptOn = False

# Number of processes used
nbOfProcesses = 6


########################################################
#                    Modalities                        #
########################################################

modalities = {}
UseInRegistration = 'UseInRegistration'
UseInRegression   = 'UseInRegression'
ModalityWeight    = 'Weight'
ModalityDataPath  = 'DataPath'

# T2 modality
T2Image = 'T2'
modalities[T2Image] = {}
modalities[T2Image][UseInRegistration] = True
modalities[T2Image][UseInRegression]   = True
modalities[T2Image][ModalityWeight]    = 1
modalities[T2Image][ModalityDataPath]  = dataPath + '/T2'

# GM modality
GMImage = 'GM'
modalities[GMImage] = {}
modalities[GMImage][UseInRegistration] = True
modalities[GMImage][UseInRegression]   = True
modalities[GMImage][ModalityWeight]    = 1
modalities[GMImage][ModalityDataPath]  = dataPath + '/GM'

# WM modality
WMImage = 'WM'
modalities[WMImage] = {}
modalities[WMImage][UseInRegistration] = True
modalities[WMImage][UseInRegression]   = True
modalities[WMImage][ModalityWeight]    = 1
modalities[WMImage][ModalityDataPath]  = dataPath + '/WM'

# Cervelet modality
CerveletImage = 'Cervelet'
modalities[CerveletImage] = {}
modalities[CerveletImage][UseInRegistration] = True
modalities[CerveletImage][UseInRegression]   = True
modalities[CerveletImage][ModalityWeight]    = 1
modalities[CerveletImage][ModalityDataPath]  = dataPath + '/Cervelet'

# Brainstem modality
BrainstemImage = 'Brainstem'
modalities[BrainstemImage] = {}
modalities[BrainstemImage][UseInRegistration] = True
modalities[BrainstemImage][UseInRegression]   = True
modalities[BrainstemImage][ModalityWeight]    = 1
modalities[BrainstemImage][ModalityDataPath]  = dataPath + '/Brainstem'

# CSF modality
CSFImage = 'CSF'
modalities[CSFImage] = {}
modalities[CSFImage][UseInRegistration] = False
modalities[CSFImage][UseInRegression]   = True
modalities[CSFImage][ModalityWeight]    = 0
modalities[CSFImage][ModalityDataPath]  = dataPath + '/CSF'

# Other modality
OtherImage = 'Other'
modalities[OtherImage] = {}
modalities[OtherImage][UseInRegistration] = False
modalities[OtherImage][UseInRegression]   = True
modalities[OtherImage][ModalityWeight]    = 0
modalities[OtherImage][ModalityDataPath]  = dataPath + '/Other'


########################################################
#                     Patients                         #
########################################################

# patients = [ ('Identifier1', age1), ('Identifier2', age2), ... , ('IdentifierN', ageN) ]
patients = [ ('ARS_Hu', 33), ('AYD_Na', 28), ('BAL_In', 32.5), ('BER_Sl', 30), ('BLA_El', 32), ('CHA_Sa', 28), ('DAH_Au', 30), ('DER_An', 26), ('ELO_Ha', 28), ('ESC_Mi', 34), ('FRE_St', 28), ('HER_Au', 32), ('HIE_Au', 28), ('KOG_Fa', 29), ('KRA_Na', 32), ('LIP_La', 33), ('MAG_Ai', 30), ('NEF_Dr', 32), ('RYC_Ca', 27), ('STR_Re', 31), ('TAB_Ar_01', 27), ('TAB_Ar_02', 30), ('TRO_Sa', 32) ]



########################################################
#                 Template creation                    #
########################################################

# Reference patient for template creation
patientReference = 'KOG_Fa'

# Working directory of template creation
templatePath = workPath + '/template'


########################################################
#                   Atlas creation                     #
########################################################

# Time step for sampling
timeStep = 1.0

# Bandwith parameter
bandwith = 1.0

# Working directory of atlas creation
atlasPath = workPath + '/atlas'


########################################################
#                 External programs                    #
########################################################

ANTS             = 'ANTS'
Warp             = 'WarpImageMultiTransform'
AverageImages    = 'AverageImages'
AverageAffines   = 'AverageAffineTransform'
AverageFields    = 'AverageDeformationFields'
ShapeRegression  = 'DeformationKernelRegression'
WeightedSum      = 'pxbinaryimageoperator'
CreateNewImage   = 'pxcreatezeroimage'
ComposeTransform = 'ComposeMultiTransform'
