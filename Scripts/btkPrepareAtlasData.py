#!/usr/bin/python
########################################################
#        Name: prepareData.py                          #      
# Description: Prepare data before atlas construction  #
#              (requires a data.py configuration file) #
#      Author: Julien Pontabry                         #
#        Date: Juin 2012                               #
########################################################


import data
import tools
import numpy
import os
import multiprocessing


# Pool used for multiprocessing
pool = multiprocessing.Pool(data.nbOfProcesses)


print 'Preparing data for atlas construction...'

########################################################
#            1. Binarize GM, WM and CSF                #
########################################################

print '\tBinarizing GM, WM and CSF...'

jobs = []

for patient in data.patients:
	path = "{0}/GM".format(data.dataPath)

	if data.scriptOn and not(os.path.isdir(path)):
		os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(data.dataPath, patient[0])
	outputImage = "{0}/GM/{1}_GM.nii.gz".format(data.dataPath, patient[0])
	goBinarize  = "btkBinarizeLabels -i {0} -o {1} -l 1 > /dev/null 2> /dev/null".format(inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/WM".format(data.dataPath)

	if data.scriptOn and not(os.path.isdir(path)):
		os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(data.dataPath, patient[0])
	outputImage = "{0}/WM/{1}_WM.nii.gz".format(data.dataPath, patient[0])
	goBinarize  = "btkBinarizeLabels -i {0} -o {1} -l 2 > /dev/null 2> /dev/null".format(inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/Brainstem".format(data.dataPath)

	if data.scriptOn and not(os.path.isdir(path)):
		os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(data.dataPath, patient[0])
	outputImage = "{0}/Brainstem/{1}_Brainstem.nii.gz".format(data.dataPath, patient[0])
	goBinarize  = "btkBinarizeLabels -i {0} -o {1} -l 3 > /dev/null 2> /dev/null".format(inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/Cervelet".format(data.dataPath)

	if data.scriptOn and not(os.path.isdir(path)):
		os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(data.dataPath, patient[0])
	outputImage = "{0}/Cervelet/{1}_Cervelet.nii.gz".format(data.dataPath, patient[0])
	goBinarize  = "btkBinarizeLabels -i {0} -o {1} -l 4 > /dev/null 2> /dev/null".format(inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/CSF".format(data.dataPath)

	if data.scriptOn and not(os.path.isdir(path)):
		os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(data.dataPath, patient[0])
	outputImage = "{0}/CSF/{1}_CSF.nii.gz".format(data.dataPath, patient[0])
	goBinarize  = "btkBinarizeLabels -i {0} -o {1} -l 5 > /dev/null 2> /dev/null".format(inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/Other".format(data.dataPath)

	if data.scriptOn and not(os.path.isdir(path)):
		os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(data.dataPath, patient[0])
	outputImage = "{0}/Other/{1}_Other.nii.gz".format(data.dataPath, patient[0])
	goBinarize  = "btkBinarizeLabels -i {0} -o {1} -l 0 > /dev/null 2> /dev/null".format(inputImage, outputImage)
	jobs.append(goBinarize)

if data.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print job

print '\tdone.'


########################################################
#              2. Blur images                          #
########################################################

print '\tBlurring images...'

jobs = []

for patient in data.patients:
	path = "{0}/GM".format(data.dataPath)
	image  = "{0}/GM/{1}_GM.nii.gz".format(data.dataPath, patient[0])
	goBlur = "btkImageGaussianFilter -i {0} -o {1} > /dev/null 2> /dev/null".format(image, image)
	jobs.append(goBlur)


	path = "{0}/WM".format(data.dataPath)
	image  = "{0}/WM/{1}_WM.nii.gz".format(data.dataPath, patient[0])
	goBlur = "btkImageGaussianFilter -i {0} -o {1} > /dev/null 2> /dev/null".format(image, image)
	jobs.append(goBlur)


	path = "{0}/Brainstem".format(data.dataPath)
	image  = "{0}/Brainstem/{1}_Brainstem.nii.gz".format(data.dataPath, patient[0])
	goBlur = "btkImageGaussianFilter -i {0} -o {1} > /dev/null 2> /dev/null".format(image, image)
	jobs.append(goBlur)


	path = "{0}/Cervelet".format(data.dataPath)
	image  = "{0}/Cervelet/{1}_Cervelet.nii.gz".format(data.dataPath, patient[0])
	goBlur = "btkImageGaussianFilter -i {0} -o {1} > /dev/null 2> /dev/null".format(image, image)
	jobs.append(goBlur)


	path = "{0}/CSF".format(data.dataPath)
	image  = "{0}/CSF/{1}_CSF.nii.gz".format(data.dataPath, patient[0])
	goBlur = "btkImageGaussianFilter -i {0} -o {1} > /dev/null 2> /dev/null".format(image, image)
	jobs.append(goBlur)


	path = "{0}/Other".format(data.dataPath)
	image  = "{0}/Other/{1}_Other.nii.gz".format(data.dataPath, patient[0])
	goBlur = "btkImageGaussianFilter -i {0} -o {1} > /dev/null 2> /dev/null".format(image, image)
	jobs.append(goBlur)

if data.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print job

print '\tdone.'


########################################################
#              3. Normalize images                     #
########################################################

print '\tNormalizing images...'

jobs = []

for patient in data.patients:
	goNorm = "btkProbabilityMapNormalization"

	path = "{0}/GM".format(data.dataPath)
	image   = "{0}/GM/{1}_GM.nii.gz".format(data.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/WM".format(data.dataPath)
	image  = "{0}/WM/{1}_WM.nii.gz".format(data.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/Brainstem".format(data.dataPath)
	image  = "{0}/Brainstem/{1}_Brainstem.nii.gz".format(data.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/Cervelet".format(data.dataPath)
	image  = "{0}/Cervelet/{1}_Cervelet.nii.gz".format(data.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/CSF".format(data.dataPath)
	image  = "{0}/CSF/{1}_CSF.nii.gz".format(data.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)

	path = "{0}/Other".format(data.dataPath)
	image  = "{0}/Other/{1}_Other.nii.gz".format(data.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)

	goNorm += " > /dev/null 2> /dev/null"
	jobs.append(goNorm)

if data.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print job

print '\tdone.'

print 'done.'
