#!/usr/bin/python
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


import btkAtlasData
import tools
import numpy
import os
import multiprocessing


# Pool used for multiprocessing
pool = multiprocessing.Pool(btkAtlasData.nbOfProcesses)

if btkAtlasData.scriptOn and not(os.path.isdir(btkAtlasData.templatePath)):
	os.mkdir(btkAtlasData.templatePath)


print 'Computing population average...'

#############################################################################
#               1.1. Estimate deformation to reference                      #
#############################################################################

print '\tEstimating deformations to reference...'

jobs = []

for patient in btkAtlasData.patients:
	if patient[0] != btkAtlasData.patientReference:
		outputImage = '{0}/{1}to{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference)
		goANTS      = '{0} 3 '.format(btkAtlasData.ANTS)

		for modality in btkAtlasData.modalities.keys():
			if btkAtlasData.modalities[modality][btkAtlasData.UseInRegistration]:
				fixedImage  = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], btkAtlasData.patientReference, modality)
				movingImage = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], patient[0], modality)
				weight      = btkAtlasData.modalities[modality][btkAtlasData.ModalityWeight]
			
				goANTS += '-m CC[{0},{1},{2},5] '.format(fixedImage, movingImage, weight)
		
		goANTS += '-o {0} -i {1} -t SyN[{2}] --use-Histogram-Matching > {3}.log 2> {3}.errlog'.format(outputImage, btkAtlasData.registrationSteps, btkAtlasData.gradientStep, outputImage+btkAtlasData.ANTS)
		jobs.append(goANTS)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t\t{0}".format(job)

print '\tdone.'


#############################################################################
#             1.2. Extract rigid part and prepare images                    #
#############################################################################

print '\tExtract rigid part from affine deformation and prepare images'

jobs = []

for patient in btkAtlasData.patients:
	if patient[0] != btkAtlasData.patientReference:
		deformation = '{0}/{1}to{2}'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference)
		affine      = deformation + 'Affine.txt'
		rigid       = deformation + 'Rigid.txt'
		affineOnly  = deformation + 'AffineOnly.txt'
		
		if btkAtlasData.scriptOn:
			# extract rigid part
			affineFile = open(affine, 'r')
			lines      = affineFile.readlines()
			parameters = lines[3].replace('\n','').replace('Parameters: ','').split(' ')
		
			affineMatrix    = tools.parametersToMatrix(parameters[0:9])
			rigidMatrix     = tools.extractRigidPartFromAffineMatrix(affineMatrix)
			rigidParameters = tools.matrixToParameters(rigidMatrix)
		
			rigidFile = open(rigid, 'w')
			lines[3]  = 'Parameters: {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}\n'.format(rigidParameters[0], rigidParameters[1], rigidParameters[2], rigidParameters[3], rigidParameters[4], rigidParameters[5], rigidParameters[6], rigidParameters[7], rigidParameters[8], parameters[9], parameters[10], parameters[11])
			for line in lines:
				rigidFile.write(line)
		
			affineFile.close()
			rigidFile.close()
		
		# compute affine only
		goCompose = '{0} 3 {1} -i {2} {3} > {4}.log 2> {4}.errlog'.format(btkAtlasData.ComposeTransform, affineOnly, rigid, affine, affineOnly+btkAtlasData.ComposeTransform)
		jobs.append(goCompose)
		
		# Warp image with rigid transform
		for modality in btkAtlasData.modalities.keys():
			if btkAtlasData.modalities[modality][btkAtlasData.UseInRegression]:
				movingImage = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], patient[0], modality)
				outputImage = '{0}/{1}to{2}Rigid_{3}.nii.gz'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference, modality)
				reference   = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], btkAtlasData.patientReference, modality)
				
				goWarp = '{0} 3 {1} {2} -R {3} {4} > {5}.log 2> {5}.errlog'.format(btkAtlasData.Warp, movingImage, outputImage, reference, rigid, outputImage+btkAtlasData.Warp)
				jobs.append(goWarp)
				
if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t\t{0}".format(job)

print '\tdone.'


#############################################################################
#                        1.3. Average Deformation                           #
#############################################################################

print '\tComputing deformation from {0} to average...'.format(btkAtlasData.patientReference)

patientReferenceToAverageAffine = '{0}/{1}toAverageAffine.txt'.format(btkAtlasData.templatePath, btkAtlasData.patientReference)
patientReferenceToAverageFields = '{0}/{1}toAverageWarp.nii.gz'.format(btkAtlasData.templatePath, btkAtlasData.patientReference)

goAverageAffine = '{0} 3 {1} '.format(btkAtlasData.AverageAffines, patientReferenceToAverageAffine)
#goAverageFields = '{0} -o {1} '.format(btkAtlasData.AverageFields, patientReferenceToAverageFields)
goAverageFields = 'btkWeightedMean -f -o {0} '.format(patientReferenceToAverageFields)

for patient in btkAtlasData.patients:
	if patient[0] != btkAtlasData.patientReference:
		deformation = '{0}/{1}to{2}'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference)
		affine      = deformation + 'AffineOnly.txt'
		field       = deformation + 'InverseWarp.nii.gz'
	
		goAverageAffine += '{0} '.format(affine)
#		goAverageFields += '{0} '.format(field)
		goAverageFields += '-i {0} '.format(field)
	
goAverageAffine += ' > {0}.log 2> {0}.errlog'.format(patientReferenceToAverageAffine+btkAtlasData.AverageAffines)
#goAverageFields += ' > {0}.log 2> {0}.errlog'.format(patientReferenceToAverageFields+btkAtlasData.AverageFields)
goAverageFields += ' > {0}.log 2> {0}.errlog'.format(patientReferenceToAverageFields)
	
if btkAtlasData.scriptOn:
	os.system(goAverageAffine)
	os.system(goAverageFields)
else:
	print "\t\t{0}".format(goAverageAffine)
	print "\t\t{0}".format(goAverageFields)
	
print '\tdone.'


#############################################################################
#                      1.4. Warp images to average                          #
#############################################################################

print '\tWarping images to average...'

jobs = []

for modality in btkAtlasData.modalities.keys():
	if btkAtlasData.modalities[modality][btkAtlasData.UseInRegression]:
		for patient in btkAtlasData.patients:
			outputImage = '{0}/{1}toAverage_{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], modality)
			reference   = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], btkAtlasData.patientReference, modality)
			fieldToRef  = '{0}/{1}to{2}Warp.nii.gz'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference)
			affineToRef = '{0}/{1}to{2}Affine.txt'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference)
				
			if patient[0] != btkAtlasData.patientReference:
				movingImage = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], patient[0], modality)
			
				goWarp = '{0} 3 {1} {2} -R {3} --use-BSpline {4} {5} {6} {7} > {8}.log 2> {8}.errlog'.format(btkAtlasData.Warp, movingImage, outputImage, reference, patientReferenceToAverageFields, patientReferenceToAverageAffine, fieldToRef, affineToRef, outputImage+btkAtlasData.Warp)
			else:
				movingImage = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], btkAtlasData.patientReference, modality)
				
				goWarp = '{0} 3 {1} {2} -R {3} --use-BSpline {4} {5} > {6}.log 2> {6}.errlog'.format(btkAtlasData.Warp, movingImage, outputImage, reference, patientReferenceToAverageFields, patientReferenceToAverageAffine, outputImage+btkAtlasData.Warp)
				
			jobs.append(goWarp)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t\t{0}".format(job)

print '\tdone.'


#############################################################################
#                       1.5. Average warped images                          #
#############################################################################

print '\tAveraging warped images...'

for modality in btkAtlasData.modalities.keys():
	if btkAtlasData.modalities[modality][btkAtlasData.UseInRegression]:
		outputImage = '{0}/Average_{1}.nii.gz'.format(btkAtlasData.templatePath, modality)
#		goAverage = '{0} 3 {1} 0 '.format(btkAtlasData.AverageImages, outputImage)
		goAverage = 'btkWeightedMean -o {0} '.format(outputImage)
		
		for patient in btkAtlasData.patients:
			image = '{0}/{1}toAverage_{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], modality)
			
#			goAverage += '{0} '.format(image)
			goAverage += '-i {0} '.format(image)
			
#		goAverage += '> {0}.log 2> {0}.errlog'.format(outputImage+btkAtlasData.AverageImages)
		goAverage += '> {0}.log 2> {0}.errlog'.format(outputImage)
	
		if btkAtlasData.scriptOn:
			os.system(goAverage)
		else:
			print "\t\t{0}".format(goAverage)

print '\tdone.'

print 'done.'


print 'Estimation deformation to average and computing template...'

#############################################################################
#                  2.1. Estimate deformation to template                    #
#############################################################################

print '\tEstimating deformations to template...'

jobs = []

for patient in btkAtlasData.patients:
	outputImage = '{0}/{1}to{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], 'Template')
	goANTS      = '{0} 3 '.format(btkAtlasData.ANTS)
	
	for modality in btkAtlasData.modalities.keys():
		if btkAtlasData.modalities[modality][btkAtlasData.UseInRegistration]:
			fixedImage  = '{0}/Average_{1}.nii.gz'.format(btkAtlasData.templatePath, modality)
			if patient[0] == btkAtlasData.patientReference:
				movingImage = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], patient[0], modality)
			else:
				movingImage = '{0}/{1}to{2}Rigid_{3}.nii.gz'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference, modality)
			weight      = btkAtlasData.modalities[modality][btkAtlasData.ModalityWeight]
				
			goANTS += '-m CC[{0},{1},{2},5] '.format(fixedImage, movingImage, weight)
			
	goANTS += '-o {0} -i {1} -t SyN[{2}] --use-Histogram-Matching > {3}.log 2> {3}.errlog'.format(outputImage, btkAtlasData.registrationSteps, btkAtlasData.gradientStep, outputImage+btkAtlasData.ANTS)
	jobs.append(goANTS)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t\t{0}".format(job)

print '\tdone.'


#############################################################################
#                    2.2. Warp images to template                           #
#############################################################################

print '\tWarping images to template...'

jobs = []

for modality in btkAtlasData.modalities.keys():
	if btkAtlasData.modalities[modality][btkAtlasData.UseInRegression]:
		for patient in btkAtlasData.patients:
			outputImage = '{0}/{1}toTemplate_{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], modality)
			reference   = '{0}/Average_{1}.nii.gz'.format(btkAtlasData.templatePath, modality)
			fieldToRef  = '{0}/{1}to{2}Warp.nii.gz'.format(btkAtlasData.templatePath, patient[0], 'Template')
			affineToRef = '{0}/{1}to{2}Affine.txt'.format(btkAtlasData.templatePath, patient[0], 'Template')
			
			if patient[0] == btkAtlasData.patientReference:
				movingImage = '{0}/{1}_{2}.nii.gz'.format(btkAtlasData.modalities[modality][btkAtlasData.ModalityDataPath], patient[0], modality)
			else:
				movingImage = '{0}/{1}to{2}Rigid_{3}.nii.gz'.format(btkAtlasData.templatePath, patient[0], btkAtlasData.patientReference, modality)
			
			goWarp = '{0} 3 {1} {2} -R {3} --use-BSpline {4} {5} > {6}.log 2> {6}.errlog'.format(btkAtlasData.Warp, movingImage, outputImage, reference, fieldToRef, affineToRef, outputImage+btkAtlasData.Warp)
			jobs.append(goWarp)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t\t{0}".format(job)

print '\tdone.'



#############################################################################
#                       2.3. Average warped images                          #
#############################################################################

print '\tAveraging warped images...'

for modality in btkAtlasData.modalities.keys():
	if btkAtlasData.modalities[modality][btkAtlasData.UseInRegression]:
		outputImage = '{0}/Template_{1}.nii.gz'.format(btkAtlasData.templatePath, modality)
#		goAverage = '{0} 3 {1} 0 '.format(btkAtlasData.AverageImages, outputImage)
		goAverage = 'btkWeightedMean -o {0} '.format(outputImage)
		
		for patient in btkAtlasData.patients:
			image = '{0}/{1}toTemplate_{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], modality)
			
#			goAverage += '{0} '.format(image)
			goAverage += '-i {0} '.format(image)
			
#		goAverage += '> {0}.log 2> {0}.errlog'.format(outputImage+btkAtlasData.AverageImages)
		goAverage += '> {0}.log 2> {0}.errlog'.format(outputImage)
	
		if btkAtlasData.scriptOn:
			os.system(goAverage)
		else:
			print "\t\t{0}".format(goAverage)

print '\tdone.'



#############################################################################
#                       2.4. Compute tissue map                             #
#############################################################################

print '\tComputing tissues map...'

outputImage = '{0}/Template_Tissues.nii.gz'.format(btkAtlasData.templatePath)
goTissues   = 'btkBinarizeTissueProbabilityMaps -o {0} --csf {1}_CSF.nii.gz --brainstem {1}_Brainstem.nii.gz --cervelet {1}_Cervelet.nii.gz --white_matter {1}_WM.nii.gz --grey_matter {1}_GM.nii.gz --other {1}_Other.nii.gz {1}{0}.log 2> {1}{0}.errlog'.format(outputImage, btkAtlasData.templatePath+'/Template')

if btkAtlasData.scriptOn:
	os.system(goTissues)
else:
	print "\t\t{0}".format(goTissues)

print '\tdone.'

print 'done.'
