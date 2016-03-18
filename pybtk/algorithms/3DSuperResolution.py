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
import nibabel
import numpy as np
from time import time
from scipy.sparse import lil_matrix

#print path.dirname( path.dirname( path.abspath(__file__) ) )
#Add path to pybtk:
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )

from pybtk.io.itk_transforms import read_itk_transform
from pybtk.filters.imagefilters import apply_affine_itk_transform_on_image
from pybtk.reconstruction.psf import compute_psf
from pybtk.reconstruction.observation import compute_H, convert_image_to_vector, convert_vector_to_image, convert_vector_to_list_images, convert_list_images_to_vector
from pybtk.reconstruction.optim import optimize, optimizebis, myOptimization

if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument('-i', '--input', help='Low-resolution image filename (required)', type=str, action='append', required = True)
  parser.add_argument('-t', '--transform', help='Transform for each input image (optional)', action='append')
  parser.add_argument('-m', '--mask', help='Mask for each input image (optional)', action='append')
  parser.add_argument('-o', '--output', help='Estimated high-resolution image filename (required)', type=str, required = True)
  parser.add_argument('--init', help='Image filename for initialization (optional)')
  parser.add_argument('-r', '--resolution', help='Resolution of the output HR image (1 for isotropic case, or 3 for anisotropic HR image). If not provided, the minimum size of the input image will be used.', action='append')
  parser.add_argument('-p', '--psf', help='3D PSF type (boxcar (default), gauss)', type=str, default='boxcar')
  parser.add_argument('--maxiter', help='Maximum number of iterations (SR optimization)', type=int, default=10)


  args = parser.parse_args()


  ###---- Print Input Information ---------------------------------------------
  print 'Number of input images: '+str(len(args.input))
  for image in args.input:
    print image
    
  if args.transform is not None:
    print 'Number of input transforms: '+str(len(args.transform)) 
    for t in args.transform:
      print t
    
    if len(args.transform) != len(args.input):
      print 'Please provide the same number of input images and transforms. Exit.\n'
      sys.exit()

  if args.mask is not None :
    print 'Number of input masks: '+str(len(args.mask))
    for m in args.mask:
      print m
    
    if len(args.mask) != len(args.input):
      print 'Please provide the same number of input images and masks. Exit.\n'
      sys.exit()
  
  print 'Output image: '+args.output
  
  if args.init is not None:
    print 'Initialization image: '+args.init
  else:
    print 'No initialization image provided.'
        
  ###---- Load input data -----------------------------------------------------
  inputImages = []
  for i in args.input:
    inputImages.append(nibabel.load(i))
  
  print 'Loading Transforms'
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
    print 'Creating mask images using 0 as a padding value'
    for i in range(len(inputImages)):
      data = np.zeros(inputImages[i].get_data().shape)
      data[np.nonzero(inputImages[i].get_data())] = 1
      maskImages.append(nibabel.Nifti1Image(data, inputImages[i].affine))         
  
  HRSpacing = []
  if args.resolution is not None :  
    if len(args.resolution) not in [1,3]:
      print 'Please provide 0, 1 or 3 values for image resolution. Exit.\n'
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
    
  print 'Resolution for image reconstruction: '
  print HRSpacing  
  
  ###---- Computing initialization if not provided ----------------------------
  if args.init is not None:
    initHRImage = nibabel.load( args.init )

  else:
    print 'Computing initialization for HR image using the stack of input images'
    #Create a reference image using the first image 
    #todo : compute the bounding box containing all input image
    LRSize    = np.float32(np.array(inputImages[0].header['dim'][1:4]))
    LRSpacing = np.float32(np.array(inputImages[0].header['pixdim'][1:4]))
    HRSize    = np.int16(np.ceil(LRSize / HRSpacing * LRSpacing))
    
    refImageData = np.zeros(HRSize, dtype=np.int16)
    s = np.eye(4)
    s[0,0] = inputImages[0].header['pixdim'][1]
    s[1,1] = inputImages[0].header['pixdim'][2]
    s[2,2] = inputImages[0].header['pixdim'][3]
    newAffine = np.dot(inputImages[0].header.get_sform(), np.linalg.inv(s))
    s[0,0] = HRSpacing[0]
    s[1,1] = HRSpacing[1]
    s[2,2] = HRSpacing[2]    
    newAffine = np.dot(newAffine,s)
    referenceImage = nibabel.Nifti1Image(refImageData, newAffine)
    
    initHRImageData = np.zeros(HRSize, dtype=np.float32)
    for i,t in zip(inputImages,inputTransforms):
      tmpImage = apply_affine_itk_transform_on_image(input_image=i,transform=t[0], center=t[1], reference_image=referenceImage, order=3) 
      initHRImageData += ( tmpImage.get_data() / np.float32(len(inputImages)) )
    initHRImage = nibabel.Nifti1Image(initHRImageData, newAffine)
    
    maskHRImageData = np.zeros(HRSize, dtype=np.float32)
    for i,t in zip(maskImages,inputTransforms):
      tmpImage = apply_affine_itk_transform_on_image(input_image=i,transform=t[0], center=t[1], reference_image=referenceImage, order=1) 
      maskHRImageData += ( tmpImage.get_data() / np.float32(len(maskImages)) )
    maskHRImage = nibabel.Nifti1Image(maskHRImageData, newAffine)
        
  #convert 3D images to stacks of 2D slice images
  #create a transform for each 2D slice
  #create a initialisation using 2D slice stacks
  #compute H for each LR image
    
  HRpsf = compute_psf(LRSpacing, HRSpacing, args.psf)
  
  x = convert_image_to_vector(initHRImage)
  maskX = convert_image_to_vector(maskHRImage)
  #Let mask the HR image
  x = x*maskX
  
  #loop over LR images and stack H, y and masks
  HList = [] 
  maskList = []
  yList = []
  for i in range(len(inputImages)):
    y = convert_image_to_vector(inputImages[i])
    m = convert_image_to_vector(maskImages[i])
    maskList.append(m)
    #-------Masked version of y
    yList.append(y*m)
    H = compute_H(inputImages[i], initHRImage, inputTransforms[i], HRpsf, maskImages[i])
    HList.append(H)
  
  #compress x and H
  index = np.nonzero(x)[0]
  xc = x[index]
  HListc = []
  for i in range(len(inputImages)):
    HListc.append(HList[i][:,index])

  
#  res = optimize(H,x,y,args.maxiter)
  #res = optimizebis(HList,x,yList,args.maxiter,args.optim)
  #outputData = res.x.reshape(initHRImage.get_data().shape)
  #decompress res.x
  #res = optimizebis(HListc,xc,yList,args.maxiter)
  #x[index] = res.x
  #print res.message
  
  res,grad = myOptimization(HListc,xc,yList,args.maxiter)
  x[index] = res
  
  outputData = x.reshape(initHRImage.get_data().shape)
  
  outputImage = nibabel.Nifti1Image(outputData, initHRImage.affine)
  nibabel.save(outputImage,args.output)

  for i in range(len(inputImages)):  
    nibabel.save(convert_vector_to_image(HList[i].dot(x),inputImages[i]),'simu_'+str(i)+'.nii.gz')
    nibabel.save(convert_vector_to_image(HList[i].dot(x)-y,inputImages[i]),'diff_'+str(i)+'.nii.gz')




      