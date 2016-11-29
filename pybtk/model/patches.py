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
from itertools import product
from sklearn.feature_extraction.image import extract_patches

def array_to_patches(arr, patch_shape=(3,3,3), extraction_step=1, normalization=False):
  #Make use of skleanr function extract_patches
  #https://github.com/scikit-learn/scikit-learn/blob/51a765a/sklearn/feature_extraction/image.py
  """Extracts patches of any n-dimensional array in place using strides.
  Given an n-dimensional array it will return a 2n-dimensional array with
  the first n dimensions indexing patch position and the last n indexing
  the patch content. 
  Parameters
  ----------
  arr : 3darray
      3-dimensional array of which patches are to be extracted
  patch_shape : integer or tuple of length arr.ndim
      Indicates the shape of the patches to be extracted. If an
      integer is given, the shape will be a hypercube of
      sidelength given by its value.
  extraction_step : integer or tuple of length arr.ndim
      Indicates step size at which extraction shall be performed.
      If integer is given, then the step is uniform in all dimensions.
  Returns
  -------
  patches : strided ndarray
      2n-dimensional array indexing patches on first n dimensions and
      containing patches on the last n dimensions. These dimensions
      are fake, but this way no data is copied. A simple reshape invokes
      a copying operation to obtain a list of patches:
      result.reshape([-1] + list(patch_shape))
  """
  
  patches = extract_patches(arr, patch_shape, extraction_step)
  patches = patches.reshape(-1, patch_shape[0],patch_shape[1],patch_shape[2])
  patches = patches.reshape(patches.shape[0], -1) 
  if normalization==True:
    patches = patches.astype(np.float32)
    patches -= np.mean(patches, axis=0)
    patches /= np.std(patches, axis=0)
  #print('%.2d patches have been extracted' % patches.shape[0])  
  return patches

def patches_to_array(patches, patch_shape, array_shape):
  #Adapted from 2D reconstruction from sklearn
  #https://github.com/scikit-learn/scikit-learn/blob/51a765a/sklearn/feature_extraction/image.py
  """
  Patches are assumed to overlap and the image is constructed by filling in
  the patches from left to right, top to bottom, averaging the overlapping
  regions.
  """
  patches = patches.reshape(len(patches),*patch_shape)
  i_x, i_y, i_z = array_shape
  p_x, p_y, p_z = patch_shape
  array = np.zeros(array_shape)
  # compute the dimensions of the patches array
  n_x = i_x - p_x + 1
  n_y = i_y - p_y + 1
  n_z = i_z - p_z + 1
  for p, (i, j, k) in zip(patches, product(range(n_x), range(n_y), range(n_z))):
      array[i:i + p_x, j:j + p_y, k:k + p_z] += p
  
  for (i, j, k) in product(range(i_x), range(i_y), range(i_z)):
      array[i, j, k] /= float(min(i + 1, p_x, i_x - i) * min(j + 1, p_y, i_y - j) * min(k + 1, p_z, i_z - k))
  return array   
  
  
