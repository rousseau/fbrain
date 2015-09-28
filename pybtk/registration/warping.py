# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 17:12:24 2015

@author: rousseau
"""
import numpy as np
import nibabel

from scipy.ndimage.interpolation import map_coordinates

from affine_transforms import transform_a_set_of_points, convert_itk_transform_to_affine_transform

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
