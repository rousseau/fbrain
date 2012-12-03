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


import btkAtlasData
import btkPatientsTools
import btkGaussianKernel
import numpy
import os
import multiprocessing


# Pool used for multiprocessing
pool = multiprocessing.Pool(btkAtlasData.nbOfProcesses)


if btkAtlasData.scriptOn and not(os.path.isdir(btkAtlasData.atlasPath)):
	os.mkdir(btkAtlasData.atlasPath)

# Minimal and maximal age of patients
minAge = btkPatientsTools.minAge(btkAtlasData.patients)
maxAge = btkPatientsTools.maxAge(btkAtlasData.patients)


#############################################################################
#                1. Kernel regression over time on shapes                   #
#############################################################################

print 'Performing kernel regression over time on shapes...'

jobs = []

for time in numpy.arange(minAge, maxAge+btkAtlasData.timeStep/2.0, btkAtlasData.timeStep):
	weightsSum = 0
	weights    = {}

	for patient in btkAtlasData.patients:
		weights[patient[0]] = btkGaussianKernel.Compute(time, float(patient[1]), btkAtlasData.bandwith)
		weightsSum += weights[patient[0]]

	outputAffine = '{0}/AtlasShapeInverse-{1:.4f}.txt'.format(btkAtlasData.atlasPath, time)
	outputField  = '{0}/AtlasShapeInverse-{1:.4f}.nii.gz'.format(btkAtlasData.atlasPath, time)

	goRegressionAffine = '{0}{1} -o {2} '.format(btkAtlasData.BtkBinaryDir, btkAtlasData.WeightedSumAffine, outputAffine)
	goRegressionField  = '{0}{1} -o {2} '.format(btkAtlasData.BtkBinaryDir, btkAtlasData.WeightedSum, outputField)
			
	for patient in btkAtlasData.patients:
		affine = '{0}/{1}toTemplateAffine.txt'.format(btkAtlasData.templatePath, patient[0])
		field  = '{0}/{1}toTemplateInverseWarp.nii.gz'.format(btkAtlasData.templatePath, patient[0])

		goRegressionAffine += '-i {0} -w {1} '.format(affine, weights[patient[0]]/weightsSum)
		goRegressionField  += '-i {0} -w {1} '.format(field, weights[patient[0]]/weightsSum)

	goRegressionAffine += ' > {0}_{1}.log 2> {0}_{1}.errlog'.format(outputAffine, btkAtlasData.WeightedSumAffine)
	goRegressionField  += ' > {0}_{1}.log 2> {0}_{1}.errlog'.format(outputField, btkAtlasData.WeightedSum)

	jobs.append(goRegressionAffine)
	jobs.append(goRegressionField)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print '\t{0}'.format(job)

print 'done.'



#############################################################################
#             2. Kernel regression over time on appearance                  #
#############################################################################

print 'Performing kernel regression over time on appearance...'

for modality in btkAtlasData.modalities.keys():
	if btkAtlasData.modalities[modality][btkAtlasData.UseInRegression]:
		for time in numpy.arange(minAge, maxAge+btkAtlasData.timeStep/2.0, btkAtlasData.timeStep):
			jobs = []
			
			weightsSum = 0
			weights    = {}
			
			# registration to time point and weight computation
			for patient in btkAtlasData.patients:
				outputImage    = '{0}/{1}toAtlas-{2:.4f}_{3}.nii.gz'.format(btkAtlasData.atlasPath, patient[0], time, modality)
				reference      = '{0}/Template_{1}.nii.gz'.format(btkAtlasData.templatePath, modality)
				fieldToSample  = '{0}/AtlasShapeInverse-{1:.4f}.nii.gz'.format(btkAtlasData.atlasPath, time)
				affineToSample = '{0}/AtlasShapeInverse-{1:.4f}.txt'.format(btkAtlasData.atlasPath, time)
				
				movingImage = '{0}/{1}toTemplate_{2}.nii.gz'.format(btkAtlasData.templatePath, patient[0], modality)
			
				goWarp = '{0}{1} 3 {2} {3} -R {4} --use-BSpline -i {5} {6} > {3}_{1}.log 2> {3}_{1}.errlog'.format(btkAtlasData.AntsBinaryDir, btkAtlasData.Warp, movingImage, outputImage, reference, affineToSample, fieldToSample)
				jobs.append(goWarp)
				
				weights[patient[0]] = btkGaussianKernel.Compute(time, float(patient[1]), btkAtlasData.bandwith)
				weightsSum += weights[patient[0]]
				
			if btkAtlasData.scriptOn:
				pool.map(os.system, jobs)
			else:
				for job in jobs:
					print "\t{0}".format(job)

			inputImage  = '{0}/Template_{1}.nii.gz'.format(btkAtlasData.templatePath, modality)
			outputImage = '{0}/Atlas-{1:.4f}_{2}.nii.gz'.format(btkAtlasData.atlasPath, time, modality)


			# weighted sum
			goSum = '{0}{1} -o {2} '.format(btkAtlasData.BtkBinaryDir, btkAtlasData.WeightedSum, outputImage)
			for patient in btkAtlasData.patients:
				inputImage  = '{0}/{1}toAtlas-{2:.4f}_{3}.nii.gz'.format(btkAtlasData.atlasPath, patient[0], time, modality)
				
				goSum += '-i {0} -w {1} '.format(inputImage, weights[patient[0]]/weightsSum)
				
			goSum += '> {0}_{1}.log 2> {0}_{1}.errlog'.format(outputImage, btkAtlasData.WeightedSum)

			if btkAtlasData.scriptOn:
				os.system(goSum)
			else:
				print '\t{0}'.format(goSum)

print 'done.'

########################################################
#             3. Compute tissues maps                  #
########################################################

print 'Computing tissues maps...'

minAge = btkPatientsTools.minAge(btkAtlasData.patients)
maxAge = btkPatientsTools.maxAge(btkAtlasData.patients)

for time in numpy.arange(minAge, maxAge+btkAtlasData.timeStep/2.0, btkAtlasData.timeStep):
	outputImage = '{0}/Atlas-{1:.4f}_Tissues.nii.gz'.format(btkAtlasData.atlasPath, time)
	inputImages = '{0}/Atlas-{1:.4f}'.format(btkAtlasData.atlasPath, time)
	goTissues   = '{0}{1} -o {2} -i {3}_Other.nii.gz -i {3}_GM.nii.gz -i {3}_WM.nii.gz -i {3}_Cervelet.nii.gz -i {3}_Brainstem.nii.gz -i {3}_CSF.nii.gz > {2}_{1}.log 2> {2}_{1}.errlog'.format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeMaps, outputImage, inputImages)

	if btkAtlasData.scriptOn:
		os.system(goTissues)
	else:
		print '\t{0}'.format(goTissues)

print 'done.'

#############################################################################
#                           3. Clean directory                             #
#############################################################################

print 'Cleaning...'

goClean = 'rm -f `ls {0}/* | grep -v -E "^{0}/Atlas-[[:digit:]]{{2}}.[[:digit:]]{{4}}_[[:alnum:]]+.nii.gz$|^{0}/AtlasShapeInverse-[[:digit:]]{{2}}.[[:digit:]]{{4}}.(nii.gz|txt)$"`'.format(btkAtlasData.atlasPath)

if btkAtlasData.scriptOn:
    os.system(goClean)
else:
    print '\t{0}'.format(goClean)

print 'done.'


