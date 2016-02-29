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

print path.dirname( path.dirname( path.abspath(__file__) ) )
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )

from pybtk.io.itk_transforms import read_itk_transform
from pybtk.filters.imagefilters import apply_affine_itk_transform_on_image
from pybtk.reconstruction.psf import compute_psf


if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument('-i', '--input', help='Low-resolution image filename (required)', type=str, action='append', required = True)
  parser.add_argument('-t', '--transform', help='Transform for each input image (optional)', action='append')
  parser.add_argument('-m', '--mask', help='Mask for each input image (optional)', action='append')
  parser.add_argument('-o', '--output', help='Estimated high-resolution image filename (required)', type=str, required = True)
  parser.add_argument('--init', help='Image filename for initialization (optional)')
  parser.add_argument('-r', '--resolution', help='Resolution of the output HR image (1 for isotropic case, or 3 for anisotropic HR image). If not provided, the minimum size of the input image will be used.', action='append')
  parser.add_argument('-p', '--psf', help='3D PSF type (boxcar (default), gauss)', type=str)

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

  #todo : mask image or use padding values
  
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
    
    nibabel.save(initHRImage,args.output)  
    
  #convert 3D images to stacks of 2D slice images
  #create a transform for each 2D slice
  #create a initialisation using 2D slice stacks
  #compute H
    
  if args.psf == None:
    psftype = 'boxcar'
  else:
    psftype = args.psf
  HRpsf = compute_psf(LRSpacing, HRSpacing, psftype)

      