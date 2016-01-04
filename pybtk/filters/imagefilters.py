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

from scipy.ndimage.interpolation import map_coordinates

from transforms import transform_a_set_of_points, convert_itk_transform_to_affine_transform

def apply_affine_itk_transform_on_image(input_image, transform, center, reference_image=None, order=None):
  """
  
  """
  input_data = np.float32(input_image.get_data())
      
  if reference_image is None:
    reference_image = input_image

  #set the interpolation order to 1 if not specified
  if order is None:
    order = 1
    
  ref_data = np.float32(reference_image.get_data())    
  
  #create index for the reference space
  i = np.arange(0,ref_data.shape[0])
  j = np.arange(0,ref_data.shape[1])
  k = np.arange(0,ref_data.shape[2])
  iv,jv,kv = np.meshgrid(i,j,k,indexing='ij')
  
  iv = np.reshape(iv,(-1))
  jv = np.reshape(jv,(-1))
  kv = np.reshape(kv,(-1))
  
  #convert the transform from itk (LPS) to nibabel (RAS) 
  nb_transform, nb_center = convert_itk_transform_to_affine_transform(transform,center)


  #compute the coordinates in the input image
  pointset = np.zeros((4,iv.shape[0]))
  pointset[0,:] = iv
  pointset[1,:] = jv
  pointset[2,:] = kv
  pointset[3,:] = np.ones((iv.shape[0]))

  #no need to specify the center here because it is included in itk-based nb_transform  
  #pointset = transform_a_set_of_points(pointset,nb_transform,reference_image.affine,np.linalg.inv(input_image.affine))
  pointset = transform_a_set_of_points(pointset,nb_transform,reference_image.affine,np.linalg.inv(input_image.affine),nb_center)

  #compute the interpolation
  val = np.zeros(iv.shape)            
  map_coordinates(input_data,[pointset[0,:],pointset[1,:],pointset[2,:]],output=val,order=order)
  
  output_data = np.reshape(val,ref_data.shape)
 
  return nibabel.Nifti1Image(output_data, reference_image.affine)
