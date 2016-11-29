# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 15:11:44 2016

@author: rousseau
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
from pybtk.reconstruction.psf import compute_psf
from pybtk.reconstruction.utilities import convert_image_to_vector, convert_vector_to_image
from pybtk.reconstruction.observation import compute_H


if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument('-r', '--ref',   help='Reference Image filename (i.e. ground truth) (required)', type=str, required = True)
  parser.add_argument('-i', '--input', help='Low-resolution image filename (required), created using btkImageResampling', type=str, required = True)
  parser.add_argument('-t', '--transform', help='Transform for each input image (optional)')
  parser.add_argument('-o', '--output', help='Estimated high-resolution image filename (required)', type=str, required = True)
  parser.add_argument('-p', '--psf', help='3D PSF type (boxcar (default), gauss)', type=str, default='boxcar')
  parser.add_argument('--padding', help='Padding value used when no mask is provided (default is 0)', type=float, default=0)

  args = parser.parse_args()
  
  HRimage = nibabel.load(args.ref)
  LRimage = nibabel.load(args.input)
  inputTransform = None
  if args.transform is not None :
    inputTransform = read_itk_transform(args.transform)
  else:
    #no transform provided : use identity as transform and zero as center
    m = np.identity(4)
    c = np.array([0, 0, 0, 1])
    inputTransform = (m,c) 

  print('Creating mask image using the following padding value:'+str(args.padding))
  data = np.zeros(HRimage.get_data().shape)
  data[HRimage.get_data() > args.padding] = 1
  maskHRImage = nibabel.Nifti1Image(data, HRimage.affine)
  print('Percentage of HR masked values : %.2f '%( np.size(np.nonzero((data))) / (1.0*np.size(data.shape)) * 100.0 / np.size(data) ) )
  data = np.zeros(LRimage.get_data().shape)
  data[LRimage.get_data() > args.padding] = 1
  maskLRImage = nibabel.Nifti1Image(data, LRimage.affine)
  print('Percentage of LR masked values : %.2f '%( np.size(np.nonzero((data))) / (1.0*np.size(data.shape)) * 100.0 / np.size(data) ) )
  
  HRSpacing = np.float32(np.array(HRimage.header['pixdim'][1:4]))  
  LRSpacing = np.float32(np.array(LRimage.header['pixdim'][1:4]))  
  psf = compute_psf(LRSpacing, HRSpacing, args.psf)
  H = compute_H(LRimage, HRimage, inputTransform, psf, maskLRImage)    
  
  x = convert_image_to_vector(HRimage)
  maskX = convert_image_to_vector(maskHRImage)
  #Let mask the HR image
  x = x*maskX

  nibabel.save(convert_vector_to_image(H.dot(x),LRimage),args.output)

  

