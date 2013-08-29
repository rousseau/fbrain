#!/usr/bin/python
# -*- coding: utf-8 -*-
#############################################################################
#
#  © Université de Strasbourg - Centre National de la Recherche Scientifique
#
#  Date: 26/09/2012
#  Author(s): Youssef Taleb, Marc Schweitzer (marc.schweitzer(at)unistra.fr)
#             François Rousseau
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
import numpy

#usage : script.py datapath subject exam age modality 

images_directory = sys.argv[1]
subject          = sys.argv[2]
exam	           = sys.argv[3]
age              = sys.argv[4]
modality         = sys.argv[5]
clean_all_temporary_files = sys.argv[6] #Set this option to 1 if you want to remove all temporary files, and to 0 otherwise

temporary_files = []


format = '.nii.gz'

#############################################################################
# Paths
#
######### Image directory :
#
#template_directory ='/your/template/directory'
template_directory = '/Users/rousseau/Data/Atlas_Fetal_London/'
#
##############################################################################
#
# image directory should respect this hierarchy
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

btk_softmask       = 'btkComputeSoftMaskUsingOrthogonalImages'
btk_intersect      = 'btkExtractMaskUsingBoundingBox'
btk_crop           = 'btkCropImageUsingMask'
btk_apply_mask     = 'btkApplyMaskToImage'
btk_translate      = 'btkTranslateImageOverTemplate'

#for MAC OSX
flirt              = 'time flirt'
convert_xfm        = 'convert_xfm'
#for LINUX
#flirt             = 'time fsl4.1-flirt' 
#convert_xfm       = 'fsl4.1-convert_xfm'


ants               ='WarpImageMultiTransform'
copy_header        ='CopyImageHeaderInformation'
image_math         ='ImageMath'
#############################################################################
# Processing
#############################################################################

print('Checking input and output directories')
input_directory = images_directory+'/'+subject+'/original/'+exam+'/nifti'
if os.path.isdir(input_directory)==False:
    print('Error : the input directory is incorrect : '+input_directory)
    exit(1)

print('Looking for input images in input directory:')
images = []
for i in os.listdir(input_directory):
    if i.endswith(format) and i.find(modality)>=0:
        images.append(i.split(".")[0])
        
print('Found input images :')
print images
    
output_directory = images_directory+'/'+subject+'/processing/'+exam+'/mask_estimation'
if os.path.isdir(output_directory)==False:
    os.mkdir(output_directory)

#1.Intersect images
inputimages = []
outputmasks =  []

print('Compute the intersection of the images \n')
for j in range(len(images)):

    inputimages.append(input_directory+'/'+images[j]+format)
    outputmasks.append(output_directory+'/'+images[j]+'_bbox_crop_mask'+format)
    temporary_files.append(output_directory+'/'+images[j]+'_bbox_crop_mask'+format)

go = btk_intersect

for k in range(len(images)):
    go += ' -i '+inputimages[k]
    go += ' -o '+outputmasks[k]

print(go)
print('\n')
os.system(go)
print('\n')

for im in images:

    #2. Apply the intersection mask
    print('\n **** Apply the intersection mask **** \n')
    inputimage = input_directory+'/'+im+format
    inputmask  = output_directory+'/'+im+'_bbox_crop_mask'+format
    outputimage= output_directory+'/'+im+'_bbox_crop'+format
    #go = btk_apply_mask + ' -i '+inputimage+' -m '+inputmask+' -o '+outputimage
    go = image_math+' --in '+inputimage+' --in '+inputmask+' --mul --out '+outputimage
    print(go)
    os.system(go)   
    temporary_files.append(outputimage)

    template=template_directory+'template-'+age+format
    template_mask=template_directory+'mask-'+age+format

    #3. Translate images over Template
    print('\n**** Translation of the images over the template **** \n')
    inputimage = output_directory+'/'+im+'_bbox_crop'+format
    outputimage= output_directory+'/'+im+'_crop_centered-'+age+format
    translation= output_directory+'/'+im+'_Template-'+age+'_translation.txt'
    go = btk_translate + ' -i '+inputimage+' -r '+template+' -o '+outputimage +' -t '+translation
    print(go)
    os.system(go)
    temporary_files.append(outputimage)
    temporary_files.append(translation)
    

    #4. Apply translation to bounding box mask
    print('\n**** Apply the translation to the bounding box mask **** \n')
    translation = output_directory+'/'+im+'_Template-'+age+'_translation.txt'
    bbox_mask   = output_directory+'/'+im+'_bbox_crop_mask'+format
    bbox_mask_centered=output_directory+'/'+im+'_crop_centered_mask-'+age+format
    #go = btk_apply_transform+' -i'
    go = ants+' 3 '+bbox_mask+' '+bbox_mask_centered+' -R '+template+' '+translation+' --use-NN'
    #print(go)
    os.system(go)
    temporary_files.append(bbox_mask_centered)

    #5. Registration with FSL
    print('\n**** Processing registration with FSL **** \n')
    inputimage = output_directory+'/'+im+'_crop_centered-'+age+format
    inputweight= bbox_mask_centered
    referenceimage = template
    referenceweight = template_mask
    outputimage = output_directory+'/'+im+'_crop_registered-'+age+format
    outputmat = output_directory+'/'+im+'_template-'+age+'_registration.mat'
    go = flirt+' -in '+inputimage+' -inweight '+inputweight+' -ref '+referenceimage+' -refweight '+referenceweight+' -out '+outputimage+' -omat '+outputmat+' -searchrx -180 180 -searchry -180 180 -searchrz -180 180'
    #print(go)
    os.system(go)
    temporary_files.append(outputimage)
    temporary_files.append(outputmat)

    #6. Inverse registration transform
    print('\n**** Inverse the FSL transform **** \n')
    inputmat = output_directory+'/'+im+'_template-'+age+'_registration.mat'
    outputmat= output_directory+'/'+im+'_inverse_template-'+age+'_registration.mat'
    go = convert_xfm+' -inverse '+inputmat+' -omat '+outputmat
    #print(go)
    os.system(go)
    temporary_files.append(outputmat)

    #7. Apply inverse registration transform to Template intracranial mask
    print('\n**** Apply the inverse transform to the template **** \n')
    transform= output_directory+'/'+im+'_inverse_template-'+age+'_registration.mat'
    referenceimage= output_directory+'/'+im+'_crop_centered-'+age+format
    outputmask = output_directory+'/'+im+'_centered_intracranial_mask-'+age+format
    go = flirt+' -in '+template_mask+' -applyxfm -init '+transform+' -out '+outputmask+' -paddingsize 0.0 -interp trilinear -ref '+referenceimage
    #print(go)
    os.system(go)
    temporary_files.append(outputmask)


    #8. Apply inverse translation to Template intracranial mask
    print('\n**** Apply the inverse translation to the template **** \n')
    inputmask= output_directory+'/'+im+'_centered_intracranial_mask-'+age+format
    transform= output_directory+'/'+im+'_Template-'+age+'_translation.txt'
    referenceimage= output_directory+'/'+im+'_bbox_crop'+format
    outputmask= output_directory+'/'+im+'_mask_estimation-'+age+format
    go = ants+' 3 '+inputmask+' '+outputmask+' -R '+referenceimage+' -i '+translation+' --use-NN'
    #print(go)
    os.system(go)
    temporary_files.append(outputmask)

    #.. Copy the header of the input image into the estimated mask
    inputimage= input_directory+'/'+im+format
    inputmask = output_directory+'/'+im+'_mask_estimation-'+age+format    
    outputmask= output_directory+'/'+im+'_mask_estimation-'+age+format    
    go = copy_header+' '+inputimage+' '+inputmask+' '+outputmask+' 1 1 1 '
    print(go)
    os.system(go)

    #.. Binarize the mask image
    inputmask = output_directory+'/'+im+'_mask_estimation-'+age+format
    outputmask= output_directory+'/'+im+'_mask_estimation-'+age+format
    go = 'btkBinarizeMask -t 0.1 -i '+inputmask+' -o '+outputmask
    os.system(go)

    #9. Crop original image with extended mask then crop extended mask with himself
    print('\n**** Crop the original image with the new template mask **** \n')
    inputimage= input_directory+'/'+im+format
    inputmask= output_directory+'/'+im+'_mask_estimation-'+age+format
    outputimage= output_directory+'/'+im+'_crop'+format
    go = btk_crop+' -i '+inputimage+' -m '+inputmask+' -o '+outputimage
    print(go)
    os.system(go)

    inputmask = output_directory+'/'+im+'_mask_estimation-'+age+format
    outputmask= output_directory+'/'+im+'_crop_mask'+format
    go = btk_crop+' -i '+inputmask+' -m '+inputmask+' -o '+outputmask
    print(go)
    os.system(go)

# CHECK CONSISTENCY OF THE RESULTS ----------------------
#10. Binarize template mask -------------
inputmask  = template_mask
outputmask = output_directory+'/binarized-template-mask.nii.gz'
go = 'btkBinarizeMask -i '+inputmask+' -o '+outputmask+' -t 0.0000001 '
print(go)
os.system(go)
temporary_files.append(outputmask)

#11. Compute image similarity between images registered on the template
mask = output_directory+'/binarized-template-mask.nii.gz'
for i in images:
  imageA = output_directory+'/'+i+'_crop_registered-'+age+format
  for j in images:
    imageB = output_directory+'/'+j+'_crop_registered-'+age+format
    outputfile = output_directory+'/nc_'+i+'_'+j+'.txt'
    go = 'btkImageSimilarity -a '+imageA+' -m '+mask+' -b '+imageB+' --nc > '+outputfile    
    print go
    os.system(go)
    temporary_files.append(outputfile) 
    
simMatrix = numpy.eye(len(images),len(images))

for i in range(len(images)):       
  for j in range(len(images)):
    filename = output_directory+'/nc_'+images[i]+'_'+images[j]+'.txt'
    L = open(filename, "r").read().splitlines();

    for line in L:
      if "NC" in line:
        simMatrix[i,j] = float(line[5:len(line)])       

print simMatrix

#12.Detection of good/consistent estimates
goodDetections = []
for i in range(len(images)):       
  for j in range(len(images)):
    if i!=j and simMatrix[i,j]>0.9:
      goodDetections.append(images[i])
      break
print 'Good Estimations are:' 
print goodDetections      

#13. Final mask estimation using only good estimates

for im in images:
  inputimages = []
  outputmasks =  []
  inputimages.append(output_directory+'/'+im+'_mask_estimation-'+age+format)
  outputmasks.append(output_directory+'/'+im+'_softfinal_mask'+format)
  for good in goodDetections:
    inputimages.append(output_directory+'/'+good+'_mask_estimation-'+age+format)
    outputmasks.append(output_directory+'/tmp'+format)
    
  go = btk_softmask+' -d 100000 -v 1 -t 2 '

  for k in range(len(inputimages)):
    go += ' -i '+inputimages[k]
    go += ' -o '+outputmasks[k]

  print(go)
  os.system(go)

temporary_files.append(output_directory+'/tmp'+format)

for im in images:
  inputmask = output_directory+'/'+im+'_softfinal_mask'+format
  outputmask= output_directory+'/'+im+'_hardfinal_mask'+format
  go = 'btkBinarizeMask -i '+inputmask+' -o '+outputmask+' -t 0.2 '
  print(go)
  os.system(go)

if clean_all_temporary_files == '1':
    print('Removing temporary files')
    for f in temporary_files:
        os.remove(f)

print('Done !')
