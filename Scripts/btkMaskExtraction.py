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
#images[('AKB_Tu','exam01','34')] = ['20090805_T2_COR_01','20090805_T2_SAG_01','20090805_T2_AXIAL_02','20090805_T2_AXIAL_01']
#images[('ARS_Hu','exam01','33')] = ['20090121_T2_AXIAL_01','20090121_T2_AXIAL_02','20090121_T2_AXIAL_03','20090121_T2_COR_01','20090121_T2_COR_02','20090121_T2_SAG_01','20090121_T2_SAG_02','20090121_T2_SAG_03']
#images[('AYD_Na','exam01','28')] = ['20090408_T2_AXIAL_01','20090408_T2_COR_01','20090408_T2_COR_02','20090408_T2_SAG_01','20090408_T2_SAG_02','20090408_T2_SAG_03','20090408_T2_SAG_04']
#images[('BAL_In','exam01','33')] = ['20060908_T2_AXIAL_01','20060908_T2_COR_01','20060908_T2_SAG_01']
#images[('BAR_Ma','exam01','37')] = ['20090624_T2_AXIAL_01','20090624_T2_AXIAL_02','20090624_T2_COR_01','20090624_T2_SAG_01','20090624_T2_SAG_02','20090624_T2_SAG_03']
#images[('BER_Sl','exam01','30')] = ['20071205_T2_AXIAL','20071205_T2_COR','20071205_T2_SAG']
#images[('BLA_El','exam01','32')] = ['20100331_T2_AXIAL_01','20100331_T2_AXIAL_02','20100331_T2_AXIAL_03','20100331_T2_AXIAL_04','20100331_T2_COR_01','20100331_T2_COR_02','20100331_T2_SAG_01']
#images[('CHA_Sa','exam01','28')] = ['20070808_T2_Axial_01','20070808_T2_Cor_01','20070808_T2_Cor_02','20070808_T2_Cor_03','20070808_T2_Sag_01']
#images[('DAH_Au','exam01','30')] = ['20060628_T2_Axial_01','20060628_T2_Axial_02','20060628_T2_Cor_01','20060628_T2_Sag_01']
#images[('DER_An','exam01','26')] = ['20090311_T2_AXIAL_01','20090311_T2_COR_01','20090311_T2_COR_02','20090311_T2_SAG_01']
#images[('ELO_Ha','exam01','28')] = ['20070613_T2_Axial_01','20070613_T2_Cor_01','20070613_T2_Sag_01']
#images[('ESC_Mi','exam01','34')] = ['20071003_T2_AXIAL_01','20071003_T2_AXIAL_02','20071003_T2_COR_01','20071003_T2_SAG_01','20071003_T2_SAG_02']
#images[('FRE_St','exam01','28')] = ['20100210_T2_AXIAL_01','20100210_T2_COR_01','20100210_T2_COR_02','20100210_T2_COR_03','20100210_T2_SAG_01','20100210_T2_SAG_02']
#images[('HER_Au','exam01','32')] = ['20091223_T2_AXIAL_01','20091223_T2_COR_01','20091223_T2_COR_02','20091223_T2_SAG_01']
#images[('HIE_Au','exam01','28')] = ['20070801_T2_Axial_01','20070801_T2_Axial_02','20070801_T2_Cor_01','20070801_T2_Sag_01']
#images[('KAR_Nu','exam01','34')] = ['20070718_T2_Axial_01','20070718_T2_Cor_01','20070718_T2_Sag_01']
#images[('KOG_Fa','exam01','29')] = ['20080806_T2_Axial_01','20080806_T2_Axial_02','20080806_T2_Cor_01','20080806_T2_Sag_01']
#images[('KRA_Na','exam01','32')] = ['20070124_T2_Axial_01','20070124_T2_Axial_02','20070124_T2_Axial_03','20070124_T2_Cor_01','20070124_T2_Cor_02','20070124_T2_Cor_03','20070124_T2_Cor_04','20070124_T2_Sag_01','20070124_T2_Sag_02','20070124_T2_Sag_03']
#images[('LEI_El','exam01','30')] = ['20100308_T2_AXIAL_01','20100308_T2_AXIAL_02','20100308_T2_AXIAL_03','20100308_T2_AXIAL_04','20100308_T2_COR_01','20100308_T2_COR_02','20100308_T2_COR_03','20100308_T2_SAG_01','20100308_T2_SAG_02']
#images[('LIP_La','exam01','33')] = ['20090429_T2_AXIAL_01','20090429_T2_COR_01','20090429_T2_SAG_01','20090429_T2_SAG_02','20090429_T2_SAG_03','20090429_T2_SAG_04']
#images[('MAG_Ai','exam01','30')] = ['20061122_T2_Axial_01','20061122_T2_Cor_01','20061122_T2_Sag_01']
#images[('MIH_Ra','exam01','30')] = ['20051123_T2_Axial_01','20051123_T2_Cor_01','20051123_T2_Sag_01']
#images[('NAV_Co','exam01','35')] = ['20080206_T2_Axial_01','20080206_T2_Axial_02','20080206_T2_Cor_01','20080206_T2_Sag_01']
#images[('NEF_Dr','exam01','27')] = ['20070110_T2_Axial_01','20070110_T2_Cor_01','20070110_T2_Sag_01']
#images[('PED_Pr','exam01','34')] = ['20070411_T2_AXIAL_01','20070411_T2_COR_01','20070411_T2_SAG_01']
#images[('RYC_Ca','exam01','27')] = ['20090225_T2_AXIAL_01','20090225_T2_COR_01','20090225_T2_SAG_01']
#images[('STR_Re','exam01','31')] = ['20080730_T2_Axial_01','20080730_T2_Cor_01','20080730_T2_Sag_01']
#images[('TAB_Ar','exam01','27')] = ['20090325_T2_AXIAL_01','20090325_T2_COR_01','20090325_T2_COR_02','20090325_T2_COR_03','20090325_T2_COR_04','20090325_T2_SAG_01']
#images[('TRO_Sa','exam01','32')] = ['20100428_T2_AXIAL_01','20100428_T2_AXIAL_02','20100428_T2_COR_01','20100428_T2_COR_02','20100428_T2_SAG_01','20100428_T2_SAG_02']
#images[('VON_Au','exam01','32')] = ['20060208_T2_Axial_01','20060208_T2_Cor_01','20060208_T2_Cor_02','20060208_T2_Cor_03','20060208_T2_Sag_01']
#images[('WAG_Ma','exam01','31')] = ['20060614_T2_Axial_01','20060614_T2_Axial_02','20060614_T2_Cor_01','20060614_T2_Sag_01','20060614_T2_Sag_02']

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
