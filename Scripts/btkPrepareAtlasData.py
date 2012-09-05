#!/usr/bin/python
# -*- coding: utf-8 -*-
#############################################################################
#
#  © Université de Strasbourg - Centre National de la Recherche Scientifique
#
#  Date: 01/06/2012
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
import numpy
import os
import multiprocessing


# Pool used for multiprocessing
pool = multiprocessing.Pool(btkAtlasData.nbOfProcesses)


#############################################################################
#            1. Binarize probability maps (GM, WM and CSF)                  #
#############################################################################

print 'Binarizing probability maps (GM, WM and CSF)...'

jobs = []

for patient in btkAtlasData.patients:
	path = "{0}/GM".format(btkAtlasData.dataPath)

	if btkAtlasData.scriptOn:
		if not(os.path.isdir(path)):
			os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(btkAtlasData.dataPath, patient[0])
	outputImage = "{0}/GM/{1}_GM.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBinarize  = "{0}{1} -i {2} -o {3} -l 1 > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeLabels, inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/WM".format(btkAtlasData.dataPath)

	if btkAtlasData.scriptOn:
		if not(os.path.isdir(path)):
			os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(btkAtlasData.dataPath, patient[0])
	outputImage = "{0}/WM/{1}_WM.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBinarize  = "{0}{1} -i {2} -o {3} -l 2 > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeLabels, inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/Brainstem".format(btkAtlasData.dataPath)

	if btkAtlasData.scriptOn:
		if not(os.path.isdir(path)):
			os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(btkAtlasData.dataPath, patient[0])
	outputImage = "{0}/Brainstem/{1}_Brainstem.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBinarize  = "{0}{1} -i {2} -o {3} -l 3 > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeLabels, inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/Cervelet".format(btkAtlasData.dataPath)

	if btkAtlasData.scriptOn:
		if not(os.path.isdir(path)):
			os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(btkAtlasData.dataPath, patient[0])
	outputImage = "{0}/Cervelet/{1}_Cervelet.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBinarize  = "{0}{1} -i {2} -o {3} -l 4 > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeLabels, inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/CSF".format(btkAtlasData.dataPath)

	if btkAtlasData.scriptOn:
		if not(os.path.isdir(path)):
			os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(btkAtlasData.dataPath, patient[0])
	outputImage = "{0}/CSF/{1}_CSF.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBinarize  = "{0}{1} -i {2} -o {3} -l 5 > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeLabels, inputImage, outputImage)
	jobs.append(goBinarize)


	path = "{0}/Other".format(btkAtlasData.dataPath)

	if btkAtlasData.scriptOn:
		if not(os.path.isdir(path)):
			os.mkdir(path)

	inputImage  = "{0}/Tissues/{1}_Tissues.nii.gz".format(btkAtlasData.dataPath, patient[0])
	outputImage = "{0}/Other/{1}_Other.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBinarize  = "{0}{1} -i {2} -o {3} -l 0 > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.BinarizeLabels, inputImage, outputImage)
	jobs.append(goBinarize)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t{0}".format(job)

print 'done.'


#############################################################################
#                      2. Blur probability maps                             #
#############################################################################

print 'Blurring probability maps...'

jobs = []

for patient in btkAtlasData.patients:
	path = "{0}/GM".format(btkAtlasData.dataPath)
	image  = "{0}/GM/{1}_GM.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBlur = "{0}{1} -i {2} -o {3} > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.GaussianFilter, image, image)
	jobs.append(goBlur)


	path = "{0}/WM".format(btkAtlasData.dataPath)
	image  = "{0}/WM/{1}_WM.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBlur = "{0}{1} -i {2} -o {3} > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.GaussianFilter, image, image)
	jobs.append(goBlur)


	path = "{0}/Brainstem".format(btkAtlasData.dataPath)
	image  = "{0}/Brainstem/{1}_Brainstem.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBlur = "{0}{1} -i {2} -o {3} > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.GaussianFilter, image, image)
	jobs.append(goBlur)


	path = "{0}/Cervelet".format(btkAtlasData.dataPath)
	image  = "{0}/Cervelet/{1}_Cervelet.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBlur = "{0}{1} -i {2} -o {3} > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.GaussianFilter, image, image)
	jobs.append(goBlur)


	path = "{0}/CSF".format(btkAtlasData.dataPath)
	image  = "{0}/CSF/{1}_CSF.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBlur = "{0}{1} -i {2} -o {3} > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.GaussianFilter, image, image)
	jobs.append(goBlur)


	path = "{0}/Other".format(btkAtlasData.dataPath)
	image  = "{0}/Other/{1}_Other.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goBlur = "{0}{1} -i {2} -o {3} > {3}_{1}.log 2> {3}_{1}.errlog".format(btkAtlasData.BtkBinaryDir, btkAtlasData.GaussianFilter, image, image)
	jobs.append(goBlur)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t{0}".format(job)

print 'done.'


#############################################################################
#                      3. Normalize probability maps                        #
#############################################################################

print 'Normalizing probability maps...'

jobs = []

for patient in btkAtlasData.patients:
	goNorm = "{0}{1}".format(btkAtlasData.BtkBinaryDir, btkAtlasData.ProbMapNormalization)

	path = "{0}/GM".format(btkAtlasData.dataPath)
	image   = "{0}/GM/{1}_GM.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/WM".format(btkAtlasData.dataPath)
	image  = "{0}/WM/{1}_WM.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/Brainstem".format(btkAtlasData.dataPath)
	image  = "{0}/Brainstem/{1}_Brainstem.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/Cervelet".format(btkAtlasData.dataPath)
	image  = "{0}/Cervelet/{1}_Cervelet.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)


	path = "{0}/CSF".format(btkAtlasData.dataPath)
	image  = "{0}/CSF/{1}_CSF.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)

	path = "{0}/Other".format(btkAtlasData.dataPath)
	image  = "{0}/Other/{1}_Other.nii.gz".format(btkAtlasData.dataPath, patient[0])
	goNorm += " -i {0} -o {0}".format(image)

	goNorm += " > {0}/{1}_{2}.log 2> {0}/{1}_{2}.errlog".format(btkAtlasData.dataPath, patient[0], btkAtlasData.ProbMapNormalization)
	jobs.append(goNorm)

if btkAtlasData.scriptOn:
	pool.map(os.system, jobs)
else:
	for job in jobs:
		print "\t{0}".format(job)

print 'done.'

