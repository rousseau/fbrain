#!/usr/bin/python
# -*- coding: utf-8 -*-
#############################################################################
#
#  © Université de Strasbourg - Centre National de la Recherche Scientifique
#
#  Date: 26/09/2012
#  Author(s): Youssef Taleb, Marc Schweitzer (marc.schweitzer(at)unistra.fr)
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


import os
import multiprocessing

numberProcess = 4 #To Change

pool = multiprocessing.Pool(numberProcess)

jobs = []
images={}



#############################################################################################
# Images
# images[('PatientName','examNum','Age(week)')] = ['image1',...,'imageN']
#############################################################################################

# Examples :
#images[('DOE_Jo','exam01','32')] = ['19910101_T2_COR_01','19910101_T2_SAG_01','19910101_T2_AXIAL_02','19910101_T2_AXIAL_01']

#############################################################################################
# Processing
#############################################################################################

for k in images.keys():

        go = 'btkMaskExtractionProcess.py '+k[0]+' '+k[1]+' '+k[2]

        for i in images[k]:
                go += ' '+i

        jobs.append( go )        
        print(go)

pool.map(os.system, jobs)
