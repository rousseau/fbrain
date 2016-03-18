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
from time import time
import numpy as np
from scipy.sparse import lil_matrix
from sklearn.preprocessing import normalize
import nibabel

from numba import jit

from os import path
import sys
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )
from pybtk.registration.affine_transforms import transform_a_set_of_points, convert_itk_transform_to_affine_transform

#Basic operation to convert vector to an image, or a list of images
def convert_image_to_vector(HRImage):
  return np.float32(HRImage.get_data()).flatten()

def convert_vector_to_image(x,HRImage):
  return nibabel.Nifti1Image(x.reshape(HRImage.get_data().shape), HRImage.affine)

def convert_list_images_to_vector(LRImages):
  y = np.array([])
  for image in LRImages:
    y = np.append(y, np.float32(image.get_data()).flatten())
  return y

def convert_vector_to_list_images(y,LRImages):
  list = []
  index=0
  for image in LRImages:
    n = image.get_data().size
    data = y[index:index+n]
    index+=n
    data = data.reshape(image.get_data().shape)
    list.append(nibabel.Nifti1Image(data.reshape(image.get_data().shape), image.affine))
  return list

def compute_H(y, x, w2wtransform, psf, mask):
  """
  Compute the matrix H, y = Hx
  y : low resolution
  x : high resolution
  
  Parameters
  ----------
  y : one nibabel image
    An array containing the lr resolution pixels (list of array) 
  x : one nibabel image
    An array containing the hr resolution pixels (array)
  w2w_transform : ndarray
    Transform from world to world coordinate (y to x), as a 2D array (homogeneous matrix)
  psf : ndarray
    PSF kernel
  mask : one nibabel image
    List of array containing points to process  
  
  Returns
  -------  
  H : sparse matrix
    Sparse matrix H corresponding to the observation model y=Hx
    
  """
  print 'Computing the observation matrix H (for each low-resolution image)'
  t = time()
  
  psfCenter = (np.array(psf.shape) -1.0 )/2.0

  LRSpacing = np.float32(np.array(y.header['pixdim'][1:4]))
  ydata  = np.float32(y.get_data())
  lookup_lr_index = np.arange(ydata.size).reshape(ydata.shape).astype(int)

  HRSize    = np.float32(np.array(x.header['dim'][1:4]))
  HRSpacing = np.float32(np.array(x.header['pixdim'][1:4]))
  xdata     = np.float32(x.get_data())
  lookup_hr_index = np.arange(xdata.size).reshape(xdata.shape)

  #Create a large empty sparse matrix
  H = lil_matrix((np.int(ydata.size),np.int(xdata.size)))
  
  #Compute the inverse of the affine matrix (from X to world coordinate)
  inverseMatrixX = np.linalg.inv(x.affine)
  #tmpH is a temporary variable to avoid many access to the sparse matrix H
  tmpH = np.zeros((xdata.size))  
    
  ratioi = HRSpacing[0]/LRSpacing[0]
  ratioj = HRSpacing[1]/LRSpacing[1]
  ratiok = HRSpacing[2]/LRSpacing[2]
  
  #Precompute PSF shift (expressed in the coordinate system of current LR image
  psfShift = np.zeros((4,psf.size))
  psfValues = np.zeros((psf.size))
  index = 0
  for ipsf in range(psf.shape[0]):
    for jpsf in range(psf.shape[1]):
      for kpsf in range(psf.shape[2]):
        psfShift[0,index] = (ipsf - psfCenter[0])*ratioi
        psfShift[1,index] = (jpsf - psfCenter[1])*ratioj
        psfShift[2,index] = (kpsf - psfCenter[2])*ratiok
        psfValues[index] = psf[ipsf,jpsf,kpsf]
        index+=1

  #convert the transform from itk (LPS) to nibabel (RAS)
  transform, center = convert_itk_transform_to_affine_transform(w2wtransform[0],w2wtransform[1])
  
  #Get index to process in the current LR image
  nonzeroMask = np.nonzero(mask.get_data())
        
  for i,j,k in zip(nonzeroMask[0],nonzeroMask[1],nonzeroMask[2]):
                
    #Compute the index of current voxel in y
    indexlr = lookup_lr_index[i,j,k]

    #Loop over PSF elements (put all psf points into psfPoints array)
    psfPoints = np.array([[i],[j],[k],[0]]) + psfShift

    #Transform index to world coordinate (y)
    #Transform world coordinate (y) to world coordinate (x)
    #Transform world coordinate (x) to index
    #These three transforms are applied using the following function:
    pointset = transform_a_set_of_points(psfPoints,transform,y.affine,inverseMatrixX,center)
    #pointset = transform_a_set_of_points(psfPoints,transform,y[l].affine,np.linalg.inv(x.affine),center)

    #Do psf weights interpolation on the grid of x (HR image)
    #we do here basic linear interpolation of the weights
    floorPoints = np.int32(np.floor(pointset))

    w1 = pointset - floorPoints
    w2 = 1 - w1

    #weights for linear interpolation
    w = np.zeros((8,psf.size))
    w[0,:] = w2[0,:] *  w2[1,:] * w2[2,:]
    w[1,:] = w1[0,:] *  w2[1,:] * w2[2,:]
    w[2,:] = w1[0,:] *  w1[1,:] * w2[2,:]
    w[3,:] = w1[0,:] *  w1[1,:] * w1[2,:]
    w[4,:] = w2[0,:] *  w1[1,:] * w1[2,:]
    w[5,:] = w1[0,:] *  w2[1,:] * w1[2,:]
    w[6,:] = w2[0,:] *  w2[1,:] * w1[2,:]
    w[7,:] = w2[0,:] *  w1[1,:] * w2[2,:]

    #To fill efficiently the sparse matrix H, we need to store nonzero elements index.
    #Indeed, finding non zero element is very long in a vector (size of x)
    tmpHList = []
    
    #Loop over all transformed PSF points to accumulate weights
    for f,weight in zip(floorPoints.T,w.T):
      if (f>=0).all() and (f[0:3]<HRSize-1).all():

        indexhr = lookup_hr_index[f[0],f[1],f[2]]
        tmpH[indexhr] += weight[0]
        tmpHList.append(indexhr)

        indexhr = lookup_hr_index[f[0]+1,f[1],f[2]]
        tmpH[indexhr] += weight[1]
        tmpHList.append(indexhr)
          
        indexhr = lookup_hr_index[f[0]+1,f[1]+1,f[2]]
        tmpH[indexhr] += weight[2]
        tmpHList.append(indexhr)

        indexhr = lookup_hr_index[f[0]+1,f[1]+1,f[2]+1]
        tmpH[indexhr] += weight[3]
        tmpHList.append(indexhr)

        indexhr = lookup_hr_index[f[0],f[1]+1,f[2]+1]
        tmpH[indexhr] += weight[4]
        tmpHList.append(indexhr)

        indexhr = lookup_hr_index[f[0]+1,f[1],f[2]+1]
        tmpH[indexhr] += weight[5]
        tmpHList.append(indexhr)

        indexhr = lookup_hr_index[f[0],f[1],f[2]+1]
        tmpH[indexhr] += weight[6]
        tmpHList.append(indexhr)

        indexhr = lookup_hr_index[f[0],f[1]+1,f[2]]
        tmpH[indexhr] += weight[7]
        tmpHList.append(indexhr)
    
    #Get the list of index to update in H
    tmpHIndex = np.unique(np.asarray(tmpHList))
    if tmpHIndex.size > 0:
      #Fill the corresponding line of the sparse matrix
      H[indexlr,tmpHIndex] = tmpH[tmpHIndex] 
      #Set tmpH to zero
      tmpH[tmpHIndex] = 0
    
  #normalization of rows (using sklearn)
  H_normalized = normalize(H, norm='l1', axis=1)

  print 'Computation done in '+str(time()-t)+' s'
  return H_normalized
