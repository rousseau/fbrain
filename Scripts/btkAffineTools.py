#!/usr/bin/python
# -*- coding: utf-8 -*-
#############################################################################
#
#  © Université de Strasbourg - Centre National de la Recherche Scientifique
#
#  Date: 01/12/2012
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


import numpy
import scipy.linalg as linalg


def extractRigidPartFromMatrix(affine):
	"Extract and return the rigid part of the affine matrix."
#	return (affine * affine.getT()) ** (-1/2) * affine
	return linalg.inv(linalg.sqrtm(affine*affine.getT())) * affine

def parametersToMatrix(parameters):
	"Convert affine parameters to a homogeneous matrix."
	s = '{0} {1} {2} ; {3} {4} {5} ; {6} {7} {8}'.format(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5], parameters[6], parameters[7], parameters[8])

	return numpy.matrix(s).getT()

def matrixToParameters(matrix):
	"Convert a homogeneous matrix to affine parameters."
	M = matrix.getT()
	params = [ str(M[0,0]), str(M[0,1]), str(M[0,2]), str(M[1,0]), str(M[1,1]), str(M[1,2]), str(M[2,0]), str(M[2,1]), str(M[2,2]) ]
	i = 0
	
	for elem in params:
		s = elem.replace('(','').replace('+0j','').replace(')','')
		params[i] = s
		i = i+1
	
	return params

