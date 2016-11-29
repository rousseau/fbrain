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
import numpy as np
import nibabel
from os import path
import sys
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )
from pybtk.reconstruction.utilities import convert_image_to_vector, convert_vector_to_image
from pybtk.filters.imagefilters import apply_affine_itk_transform_on_image

def ibpComputeError(x, H, y, hrMaskSum, lrImages, transforms, interpOrder):
  hrError = np.zeros(hrMaskSum.get_data().shape, dtype=np.float32)
  #Loop over LR image to compute the sum of errors
  for i in range(len(lrImages)):
    lrError = convert_vector_to_image(H[i].dot(x)-y[i], lrImages[i])
    tmp2 = apply_affine_itk_transform_on_image(input_image = lrError, transform=transforms[i][0], center=transforms[i][1], reference_image=hrMaskSum, order=interpOrder)
    hrError += tmp2.get_data()
  
  #Normalize the error image  
  hrErrorNorm = np.zeros(hrMaskSum.get_data().shape, dtype=np.float32)
  index = np.nonzero(hrMaskSum.get_data())
  hrErrorNorm[index] = hrError[index] / hrMaskSum.get_data()[index]
  
  return hrErrorNorm
  
def iterativeBackPropagation(hrImage, lrImages, lrMasks, transforms, H, itermax, interpOrder):
  
  #Convert LR images to a list of vectors
  y = []
  for i in range(len(lrImages)):
    y.append(convert_image_to_vector(lrImages[i]) * convert_image_to_vector(lrMasks[i]))

  #Convert HR Image to vector
  x = convert_image_to_vector(hrImage)
  outputImage = nibabel.Nifti1Image(hrImage.get_data(),hrImage.affine) 

  #Compute HR mask
  hrMaskSum=np.zeros(hrImage.get_data().shape, dtype=np.float32)
  for i in range(len(lrImages)):
    tmp1 = apply_affine_itk_transform_on_image(input_image = lrMasks[i], transform=transforms[i][0], center=transforms[i][1], reference_image=hrImage, order=0)
    hrMaskSum += tmp1.get_data()  
  
  
  #index = np.nonzero(hrMaskSum)
  
  for j in range(itermax):
     
    error = ibpComputeError(x, H, y, nibabel.Nifti1Image(hrMaskSum,hrImage.affine), lrImages, transforms, interpOrder)
#    #simulation and error computation
#    hrError = np.zeros(hrImage.get_data().shape, dtype=np.float32)
#
#    for i in range(len(lrImages)):
#      lrError = convert_vector_to_image(H[i].dot(x)-y[i], lrImages[i])
#      tmp2 = apply_affine_itk_transform_on_image(input_image = lrError, transform=transforms[i][0], center=transforms[i][1], reference_image=hrImage, order=interpOrder)
#      hrError += tmp2.get_data()
#    
#    hrError2 = np.zeros(hrImage.get_data().shape, dtype=np.float32)
#    hrError2[index] = hrError[index] / hrMaskSum[index]
#    
    #filter error map
    from skimage.restoration import denoise_tv_chambolle
    hrError2 = denoise_tv_chambolle(error,weight=5)
    
    #update hr image and x
    outputImage = nibabel.Nifti1Image(outputImage.get_data() - hrError2,hrImage.affine)
    nibabel.save(nibabel.Nifti1Image(hrError2,hrImage.affine),'error_iter'+str(j)+'.nii.gz')
    nibabel.save(outputImage,'ibp_iter'+str(j)+'.nii.gz')    
    x = convert_image_to_vector(outputImage)
  
  return outputImage  