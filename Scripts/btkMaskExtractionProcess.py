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


import sys
import os

#usage : script.py subject exam modality image01 image02 ...


subject   = sys.argv[1]
exam	  = sys.argv[2]
age       = sys.argv[3]
images	  = sys.argv[4:len(sys.argv)]
inputimages = []
outputmasks =  []


format = '.nii.gz'

#############################################################################
# Paths
#
#image directory
#images_directory = '/your/image/directory/'
#template_directory ='/your/template/directory'
images_directory = '/home/miv/schweitzer/Donnees_a_traiter/'
template_directory ='/home/miv/schweitzer/Donnees_a_traiter/Templates2/'
# image directory should respect this hierrarchy
#                           |
#                     [PatientName]
#                  /                 \
#              [original]            [processing]
#             /        |                |       \
#         [exam01]   [...]            [exam01]     [...]
#        /   |                     /     |      \
#  [dicom] [nifti]      [diffusion][segmentation][reconstruction]
#
#############################################################################

#############################################################################
# Programs
# Be sure that you have btk, fsl & ants programs installed, and added in your PATH
#############################################################################

btk_intersect      = 'btkExtractMaskUsingBoundingBox'
btk_crop           = 'btkCropImageUsingMask'
btk_apply_mask     = 'btkApplyMaskToImage'
btk_translate      = 'btkTranslateImageOverTemplate'

flirt              = 'time fsl4.1-flirt'
convert_xfm        = 'fsl4.1-convert_xfm'

ants               ='WarpImageMultiTransform'


#############################################################################
# Processing
#############################################################################
#1.Intersect images
print('Intersection of the images \n')
for j in range(len(images)):

 inputimages.append(images_directory+subject+'/original/'+exam+'/nifti/'+images[j]+format)
 outputmasks.append(images_directory+subject+'/processing/'+exam+'/segmentation/'+images[j]+'_crop_mask'+format)

go = btk_intersect

for k in range(len(images)):

 go += ' -i '+inputimages[k]

for k in range(len(images)):

 go += ' -o '+outputmasks[k]

print(go)
print('\n')
os.system(go)
print('\n')

for im in images:

#2. Crop images with outputmask
 print('\n **** Crop of the images **** \n')
 inputimage = images_directory+subject+'/original/'+exam+'/nifti/'+im+format
 inputmask  = images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_crop_mask'+format
 outputimage= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop'+format
 go = btk_apply_mask + ' -i '+inputimage+' -m '+inputmask+' -o '+outputimage
 #print(go)
 os.system(go)


 template=template_directory+'template-'+age+format
 template_mask=template_directory+'mask-'+age+format

#3.Translate images over Template
 print('\n**** Translation of the images over the template **** \n')
 inputimage = images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop'+format
 outputimage= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop_centered-'+age+format
 translation= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_Template-'+age+'_translation.txt'
 go = btk_translate + ' -i '+inputimage+' -r '+template+' -o '+outputimage +' -t '+translation
 print(go)
 os.system(go)

#4. Apply translation to bouding box mask
 print('\n**** Apply the translation to the bounding box mask **** \n')
 translation= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_Template-'+age+'_translation.txt'
 bbox_mask=images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_crop_mask'+format
 bbox_mask_centered=images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_crop_centered_mask-'+age+format
 #go = btk_apply_transform+' -i'
 go = ants+' 3 '+bbox_mask+' '+bbox_mask_centered+' -R '+template+' '+translation+' --use-NN'
 #print(go)
 os.system(go)


#5. Registration with FSL
 print('\n**** Processing registration with FSL **** \n')
 inputimage = images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop_centered-'+age+format
 inputweight= bbox_mask_centered
 referenceimage = template
 referenceweight = template_mask
 outputimage = images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop_registered-'+age+format
 outputmat = images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_template-'+age+'_registration.mat'
 go = flirt+' -in '+inputimage+' -inweight '+inputweight+' -ref '+referenceimage+' -refweight '+referenceweight+' -out '+outputimage+' -omat '+outputmat+' -searchrx -180 180 -searchry -180 180 -searchrz -180 180'
 #print(go)
 os.system(go)


#6.Inverse registration transform
 print('\n**** Inverse the FSL transform **** \n')
 inputmat = images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_template-'+age+'_registration.mat'
 outputmat = images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_inverse_template-'+age+'_registration.mat'
 go = convert_xfm+' -inverse '+inputmat+' -omat '+outputmat
 #print(go)
 os.system(go)


#7. Apply inverse registration transform to Template intracranial mask
 print('\n**** Apply the inverse transform to the template **** \n')
 transform= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_inverse_template-'+age+'_registration.mat'
 referenceimage= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop_centered-'+age+format
 outputmask = images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_centered_intracranial_mask-'+age+format
 go = flirt+' -in '+template_mask+' -applyxfm -init '+transform+' -out '+outputmask+' -paddingsize 0.0 -interp trilinear -ref '+referenceimage
 #print(go)
 os.system(go)


#8. Apply inverse translation to Template intracranial mask
 print('\n**** Apply the inverse translation to the template **** \n')
 inputmask= images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_centered_intracranial_mask-'+age+format
 transform= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_Template-'+age+'_translation.txt'
 referenceimage= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop'+format
 outputmask= images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_mask_estimation-'+age+format
 go = ants+' 3 '+inputmask+' '+outputmask+' -R '+referenceimage+' -i '+translation+' --use-NN'
 #print(go)
 os.system(go)



#9. Crop original image with extended mask then crop extended mask with himself
 print('\n**** Crop the original image with the new template mask **** \n')
 inputimage= images_directory+subject+'/original/'+exam+'/nifti/'+im+format
 inputmask= images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_mask_estimation-'+age+format
 outputimage= images_directory+subject+'/processing/'+exam+'/reconstruction/'+im+'_crop_final-'+age+format
 go = btk_crop+' -i '+inputimage+' -m '+inputmask+' -o '+outputimage
 #print(go)
 os.system(go)


 inputmask = images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_mask_estimation-'+age+format
 outputmask= images_directory+subject+'/processing/'+exam+'/segmentation/'+im+'_mask_estimation_final-'+age+format
 go = btk_crop+' -i '+inputmask+' -m '+inputmask+' -o '+outputmask
 #print(go)
 os.system(go)

print('Done !')

