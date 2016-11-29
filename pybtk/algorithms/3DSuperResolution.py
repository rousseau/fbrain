# -*- coding: utf-8 -*-
"""

  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.

"""

import sys
import argparse
from os import path
import os
import nibabel
import numpy as np

#print path.dirname( path.dirname( path.abspath(__file__) ) )
#Add path to pybtk:
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )

from pybtk.io.itk_transforms import read_itk_transform
from pybtk.filters.imagefilters import apply_affine_itk_transform_on_image,gaussian_biais_correction,apply_N4_on_image
from pybtk.reconstruction.psf import compute_psf
from pybtk.reconstruction.utilities import convert_image_to_vector, convert_vector_to_image, convert_vector_to_list_images, convert_list_images_to_vector
from pybtk.reconstruction.observation import compute_H
from pybtk.reconstruction.optim import optimization
from pybtk.reconstruction.ibp import iterativeBackPropagation


if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument('-i', '--input', help='Low-resolution image filename (required)', type=str, action='append', required = True)
  parser.add_argument('-t', '--transform', help='Transform for each input image (optional)', action='append')
  parser.add_argument('-m', '--mask', help='Mask for each input image (optional)', action='append')
  parser.add_argument('-o', '--output', help='Estimated high-resolution image filename (required)', type=str, required = True)
  parser.add_argument('--init', help='Image filename for initialization (required)', type=str, required = True)
  parser.add_argument('-r', '--resolution', help='Resolution of the output HR image (1 for isotropic case, or 3 for anisotropic HR image). If not provided, the minimum size of the input image will be used.', action='append')
  parser.add_argument('-p', '--psf', help='3D PSF type (boxcar (default), gauss)', type=str, default='boxcar')
  parser.add_argument('--maxiter', help='Maximum number of iterations (SR optimization)', type=int, default=10)
  parser.add_argument('--padding', help='Padding value used when no mask is provided', type=float, default=0)
  parser.add_argument('--bias', help='Do bias correction (N4) + local intensity correction', type=bool, default = False)

  args = parser.parse_args()


  ###---- Print Input Information ---------------------------------------------
  print('Number of input images: ',str(len(args.input)))
  for image in args.input:
    print(image)
    
  if args.transform is not None:
    print('Number of input transforms: ',str(len(args.transform)))
    for t in args.transform:
      print(t)
    
    if len(args.transform) != len(args.input):
      print('Please provide the same number of input images and transforms. Exit.\n')
      sys.exit()

  if args.mask is not None :
    print('Number of input masks: ',str(len(args.mask)))
    for m in args.mask:
      print(m)
    
    if len(args.mask) != len(args.input):
      print('Please provide the same number of input images and masks. Exit.\n')
      sys.exit()
  
  print('Output image: ',args.output)
  
  if args.init is not None:
    print('Initialization image: ',args.init)
  else:
    print('No initialization image provided.')
        
  ###---- Load input data -----------------------------------------------------
  inputImages = []
  for i in args.input:
    inputImages.append(nibabel.load(i))
  
  print('Loading Transforms')
  print('Warning: Transforms depend on the HR image.') 
  inputTransforms = []
  if args.transform is not None :
    for t in args.transform:
      inputTransforms.append(read_itk_transform(t))
  else:
    #no transform provided : use identity as transform and zero as center
    m = np.identity(4)
    c = np.array([0, 0, 0, 1])
    for i in args.input:
      inputTransforms.append( (m,c) )

  maskImages = []
  if args.mask is not None :
    for i in args.mask:
      maskImages.append(nibabel.load(i))

  else:
    print('Creating mask images using the following padding value:',str(args.padding))
    for i in range(len(inputImages)):
      data = np.zeros(inputImages[i].get_data().shape)
      data[inputImages[i].get_data() > args.padding] = 1
      maskImages.append(nibabel.Nifti1Image(data, inputImages[i].affine)) 
      #np.nonzero returns index array, so needs to be divide by the dimension of the array (i.e. 3 here)
      print('Percentage of masked values : %.2f '%( np.size(np.nonzero((data))) / (1.0*np.size(data.shape)) * 100.0 / np.size(data) )   )     
  
  HRSpacing = []
  if args.resolution is not None :  
    if len(args.resolution) not in [1,3]:
      print('Please provide 0, 1 or 3 values for image resolution. Exit.\n')
      sys.exit()  
    if len(args.resolution) == 1:
      r = np.array(float(args.resolution[0]))
      HRSpacing = np.array([1, 1, 1]) * r
    else:
      r = [float(i) for i in args.resolution]
      HRSpacing = np.array(r)
      
  else:
    r = float(min(inputImages[0].header['pixdim'][1:4]))
    HRSpacing = np.array([1, 1, 1]) * r
    
  print('Resolution for image reconstruction: ', HRSpacing  )
  
  ###---- Reading HR image as initialization ----------------------------
  initHRImage = nibabel.load( args.init )
  data = np.zeros(initHRImage.get_data().shape)
  data[initHRImage.get_data() > args.padding] = 1
  maskHRImage = nibabel.Nifti1Image(data, initHRImage.affine)
        
  #convert 3D images to stacks of 2D slice images
  #create a transform for each 2D slice
  #create a initialisation using 2D slice stacks
  #compute H for each LR image

  psfList = []
  for i in range(len(inputImages)):
    LRSpacing = np.float32(np.array(inputImages[i].header['pixdim'][1:4]))
    HRpsf = compute_psf(LRSpacing, HRSpacing, args.psf)
    psfList.append(HRpsf)
  
  #Compute H
  HList = []
  for i in range(len(inputImages)):
    HList.append( compute_H(inputImages[i], initHRImage, inputTransforms[i], psfList[i], maskImages[i]) )
   
  #Intensity correction To do
  #N4 on initHR
  #local correction
  #New init HR
  if args.bias == True:
    initHRImage_N4 = apply_N4_on_image(initHRImage, shrink_factor=1)
      
    xN4 = convert_image_to_vector(initHRImage_N4)
    hrN4Data = np.zeros(initHRImage.get_data().shape)
    for i in range(len(inputImages)):
      simu = convert_vector_to_image(HList[i].dot(xN4),inputImages[i])
      im = gaussian_biais_correction(inputImages[i],simu, 5)
        
      warped = apply_affine_itk_transform_on_image(input_image=im,transform=inputTransforms[i][0], center=inputTransforms[i][1], reference_image=initHRImage, order=3)
      hrN4Data += (warped.get_data() / np.float32(len(inputImages)) )
    initHRImage = nibabel.Nifti1Image(hrN4Data, initHRImage.affine)  
             
  #Compute x
  x = convert_image_to_vector(initHRImage)
  maskX = convert_image_to_vector(maskHRImage)
  #Let mask the HR image
  x = x*maskX

  #loop over LR images and stack y and masks
  maskList = []
  yList = []
  for i in range(len(inputImages)):
    y = convert_image_to_vector(inputImages[i])
    m = convert_image_to_vector(maskImages[i])
    maskList.append(m)
    #-------Masked version of y
    yList.append(y*m)
    
#  outputData = optimization(HList,x,yList,args.maxiter,initHRImage)
#  outputData = optimization(HList,x,yList,100,initHRImage)  
#  
#  tmpImage = initHRImage
#  for j in range(10):
#    outputData = optimization(HList,x,yList,100,tmpImage)
#    tmpImage = nibabel.Nifti1Image(outputData, initHRImage.affine)
#    nibabel.save(tmpImage,'toto_'+str(j)+'.nii.gz')
#  
#  outputImage = nibabel.Nifti1Image(outputData, initHRImage.affine)
  
  outputImage = iterativeBackPropagation(initHRImage, inputImages, maskImages, inputTransforms, HList, 10, 3)  
  
  nibabel.save(outputImage,args.output)

  for i in range(len(inputImages)):
  
    nibabel.save(convert_vector_to_image(HList[i].dot(x),inputImages[i]),'simu_'+str(i)+'.nii.gz')
    nibabel.save(convert_vector_to_image(HList[i].dot(x)-yList[i],inputImages[i]),'diff_'+str(i)+'.nii.gz')


  from skimage.restoration import denoise_tv_chambolle
  current = initHRImage
  for i in range(10):
    w = 5
    ibp = iterativeBackPropagation(current, inputImages, maskImages, inputTransforms, HList, 1, 3)
    current = nibabel.Nifti1Image(ibp.get_data(), initHRImage.affine)
    nibabel.save(current,'toto_ibp_w'+str(w)+'_'+str(i)+'.nii.gz')

  current = initHRImage
  for i in range(10):
    w = 5
    ibp = iterativeBackPropagation(current, inputImages, maskImages, inputTransforms, HList, 1, 3)
    cham = denoise_tv_chambolle(ibp.get_data(),weight=w)
    current = nibabel.Nifti1Image(cham, initHRImage.affine)
    nibabel.save(current,'toto_ibp_cham_w'+str(w)+'_'+str(i)+'.nii.gz')

  current = initHRImage
  for i in range(10):
    w = 5
    ibp = iterativeBackPropagation(current, inputImages, maskImages, inputTransforms, HList, 1, 3)
    cham = denoise_tv_chambolle(ibp.get_data(),weight=w)
    current = nibabel.Nifti1Image(cham, initHRImage.affine)
    nibabel.save(current,'toto_ibp_cham_local_w'+str(w)+'_'+str(i)+'.nii.gz')
    for i in range(len(inputImages)):
      xcurrent = convert_image_to_vector(current)
      simu = convert_vector_to_image(HList[i].dot(xcurrent),inputImages[i])
      im = gaussian_biais_correction(inputImages[i],simu, 5)
      inputImages[i] = im
      #yList[i] = convert_image_to_vector(im)

      
      
  for i in range(len(inputImages)):
    x = convert_image_to_vector(current)
    im = convert_vector_to_image(HList[i].dot(x)-yList[i],inputImages[i])
    nibabel.save(im,'diff_'+str(i)+'.nii.gz')
    warped = apply_affine_itk_transform_on_image(input_image=im,transform=inputTransforms[i][0], center=inputTransforms[i][1], reference_image=initHRImage, order=3)
    nibabel.save(im,'err_'+str(i)+'.nii.gz')



      