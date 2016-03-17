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
import math

def euler_to_matrix(ax,ay,az):
  """
  compute the homogeneous rotation matrix from Euler angles using the xyz axis
  """
  mat = np.identity(4)
  
  cosx = math.cos(math.pi*ax/180)
  cosy = math.cos(math.pi*ay/180)
  cosz = math.cos(math.pi*az/180)
  sinx = math.sin(math.pi*ax/180)
  siny = math.sin(math.pi*ay/180)
  sinz = math.sin(math.pi*az/180)
  
  mat[0,0] = cosy * cosz
  mat[0,1] = cosy * sinz 
  mat[0,2] = - siny
  mat[1,0] = sinx * siny * cosz - cosx * sinz
  mat[1,1] = sinx * siny * sinz + cosx * cosz 
  mat[1,2] = sinx * cosy 
  mat[2,0] = cosx * siny * cosz  + sinx * sinz 
  mat[2,1] = cosx * siny * sinz  - sinx * cosz
  mat[2,2] = cosx * cosy 

  return mat
  

def compute_affine_matrix_by_composition(translation=None, angles=None, scale=None, shear=None):
  """
  compute the homogeneous affine matrix by composing each matrix
  slower computation than direct computation
  """
  mat = np.identity(4)

  if translation is not None:
    mat[:3,3] = translation[:3]
  if angles is not None:
    rot = euler_to_matrix(angles[0], angles[1], angles[2])
    mat = np.dot(mat,rot)
  if shear is not None:
    sh = np.identity(4)
    sh[0,1] = shear[0]
    sh[0,2] = shear[1]
    sh[1,2] = shear[2]
    mat = np.dot(mat,sh)
  if scale is not None:
    sc = np.identity(4)
    sc[0,0] = scale[0]
    sc[1,1] = scale[1]
    sc[2,2] = scale[2]
    mat = np.dot(mat,sc)
    
  return mat  
  
  
def compute_affine_matrix(translation=None, angles=None, scale=None, shear=None):
  """
  compute the affine matrix using a direct computation
  faster computation than numpy matrix multiplication
  """
  
  mat = np.identity(4)
  gx,gy,gz = 0.0,0.0,0.0
  sx,sy,sz = 1.0,1.0,1.0
  if translation is not None:
    mat[:3,3] = translation[:3]
  if angles is not None:
    ax = math.pi*angles[0]/180.0
    ay = math.pi*angles[1]/180.0
    az = math.pi*angles[2]/180.0
    cosx = math.cos(ax)
    cosy = math.cos(ay)
    cosz = math.cos(az)
    sinx = math.sin(ax)
    siny = math.sin(ay)
    sinz = math.sin(az)
  if shear is not None:
    gx = shear[0]
    gy = shear[1]
    gz = shear[2]
  if scale is not None:
    sx = scale[0]
    sy = scale[1]
    sz = scale[2]
    
  mat[0,0] = sx * cosy * (cosz + (gy*sinz) )
  mat[0,1] = sy * (cosy * (sinz + (gx * gy * cosz)) - (gz * siny) )
  mat[0,2] = sz * ( (gx * cosy * cosz) - siny)
  mat[1,0] = sx * (sinx * siny * (cosz + gy * sinz) - cosx * (sinz + (gy * cosz) ))
  mat[1,1] = sy * (sinx * siny * (sinz + (gx * gz * cosz) ) + cosx * (cosz - (gx * gy * sinz)) + (gz * sinx * cosy))
  mat[1,2] = sz * (sinx * cosy + (gx * (sinx * siny * cosz - cosx * sinz)))
  mat[2,0] = sx * (cosx * siny * (cosz + (gy * sinz)) + sinx * (sinz - (gy * cosz) ))
  mat[2,1] = sy * (cosx * siny * (sinz + (gx * gz * cosz)) - sinx * (cosz - (gx * gz * sinz)) + (gz * cosx * cosy) )
  mat[2,2] = sz * (cosx * cosy + (gx * ( (cosx * siny * cosz) + (sinx * sinz) )) )
  
  return mat


def transform_a_point(point, transform, matrix1=None, matrix2=None, center=None):
  """
  transform a point from image 1 to image 2 using the input transform matrix
  
  Parameters
  ----------
  point : array
    An array containing the point coordinates in image 1
  transform : ndarray
    Transform expressed in world coordinate, as a 2D array (homogeneous matrix)
  matrix1 : ndarray
    Transform from image 1 to world coordinate, as a 2D array (homogeneous matrix)
  matrix2 : ndarray
    Transform from world coordinate to image 2, as a 2D array (homogeneous matrix)
  center : array
    An array containing the coordinate of the center of rotation, expressed in world coordinate
  
  Returns
  -------  
  output : array
    An array containing the output point coordinates in image 2
    
  """
  output = point
  
   #from image 1 coordinate to world coordinate
  if matrix1 is not None:
    output = np.dot(matrix1,output)
  
  #apply a translation with respect to the center of the rotation if specified    
  if center is not None:
    output[:3] = output[:3] - center[:3]
    
  #apply the transform
  output = np.dot(transform,output)

  #apply the inverse translation wrt the center of the rotation if specified
  if center is not None:
    output[:3] = output[:3] + center[:3]
    
  #from world coordinate to image 2 coordinate
  if matrix2 is not None:
    output = np.dot(matrix2,output)
  
  return output  
    
def transform_a_set_of_points(points, transform, matrix1=None, matrix2=None, center=None):
  """
  same as transform_a_point but used for a set of points
  
  Parameters
  ----------
  points : ndarray
    An 2D array containing coordinates (each column is a point)
  transform : ndarray
    Transform expressed in world coordinate, as a 2D array (homogeneous matrix)
  matrix1 : ndarray
    Transform from image 1 to world coordinate, as a 2D array (homogeneous matrix)
  matrix2 : ndarray
    Transform from world coordinate to image 2, as a 2D array (homogeneous matrix)
  center : array
    An array containing the coordinate of the center of rotation, expressed in world coordinate
  
  Returns
  -------  
  output : ndarray
    An 2D array containing the output coordinates of the set of points in image 2
  
  """
  output = points
  
   #from image 1 coordinate to world coordinate
  if matrix1 is not None:
    output = np.dot(matrix1,output)
  
  #apply a translation with respect to the center of the rotation if specified    
  if center is not None:
    output[:3] = output[:3] - center[:3,np.newaxis]    

  #apply the transform
  output = np.dot(transform,output)

  #apply the inverse translation wrt the center of the rotation if specified
  if center is not None:
    output[:3] = output[:3] + center[:3,np.newaxis]
    
  #from world coordinate to image 2 coordinate
  if matrix2 is not None:
    output = np.dot(matrix2,output)
  
  return output  
  
  
def read_itk_transform( transform_file ):
  '''
  modified from https://gist.github.com/haehn/5614966
  '''

  print 'reading the following ITK transform:'
  print transform_file
  
  # read the transform
  transform = None
  with open( transform_file, 'r' ) as f:
    for line in f:

      # check for Parameters:
      if line.startswith( 'Transform:' ):
        print 'type for transform:'
        print line.split( ': ' )[1]

      if line.startswith( 'Parameters:' ):
        values = line.split( ': ' )[1].split( ' ' )

        # filter empty spaces and line breaks
        values = [float( e ) for e in values if ( e != '' and e != '\n' )]
        # create the upper left of the matrix
        transform_upper_left = np.reshape( values[0:9], ( 3, 3 ) )
        # grab the translation as well
        translation = values[9:]

      # check for FixedParameters:
      if line.startswith( 'FixedParameters:' ):
        values = line.split( ': ' )[1].split( ' ' )

        # filter empty spaces and line breaks
        values = [float( e ) for e in values if ( e != '' and e != '\n' )]
        # setup the center
        center = np.vstack( (np.reshape(values,(3,1)), [1]) )

  center = np.reshape(center,(4,))
  # add the [0, 0, 0] line
  transform = np.vstack( ( transform_upper_left, [0, 0, 0] ) )
  # and the [offset, 1] column
  
  transform = np.hstack( ( transform, np.vstack( (np.reshape(translation,(3,1)), [1]) ) ) )

  return (transform,center)

def convert_itk_transform_to_affine_transform_without_output_center(transform,center):
  """
  Convert a LPS ITK transform into a RAS nibabel transform
  
  Parameters
  ----------
  transform : ndarray
    Transform expressed in ITK world coordinate (LPS), as a 2D array (homogeneous matrix)
  center : array
    An array containing the coordinate of the center of rotation, expressed in world coordinate
  
  Returns
  -------  
  matrix : ndarray
    Transform expressed in nibabel world coordinate (RAS), as a 2D array (homogeneous matrix)
    Note : the center of the output transform is then (0,0,0)
 
  """
  matrix = np.identity(4)

  #from itk transform, compute an affine transform for RAS nibabel data
  RAS2LPS = np.diag([-1, -1, 1, 1])
  
  #copy rotation matrix
  matrix[0:3,0:3] = transform[0:3,0:3]
  #compute offset taking into account center of rotation
  offset = center+transform[:,3] 
  offset[3] = 1
  offset[0:3] = offset[0:3] - np.dot(transform[0:3,0:3],np.reshape(center[0:3],(3,)))
  matrix[:,3] = offset
  
  matrix = np.dot(RAS2LPS, np.dot(matrix,RAS2LPS) )

  return matrix  
  
def convert_itk_transform_to_affine_transform(transform,center):  
  """
  Convert a LPS ITK transform into a RAS nibabel transform
  
  Parameters
  ----------
  transform : ndarray
    Transform expressed in ITK world coordinate (LPS), as a 2D array (homogeneous matrix)
  center : array
    An array containing the coordinate of the center of rotation, expressed in world coordinate
  
  Returns
  -------  
  RAStransform : ndarray
    Transform expressed in nibabel world coordinate (RAS), as a 2D array (homogeneous matrix)
  RAScenter : array
    An array containing the coordinate of the center of rotation, expressed in world coordinate (RAS)   
  """
  #from itk transform, compute an affine transform for RAS nibabel data
  RAS2LPS = np.diag([-1, -1, 1, 1])
  
  #compute new matrix in RAS space
  RAStransform = np.dot(RAS2LPS, np.dot(transform, RAS2LPS))
  #compute new center in RAS space
  RAScenter = np.dot(RAS2LPS,center)

  return RAStransform, RAScenter  

def convert_affine_transform_to_itk_transform(transform,center):
  """
  Convert a RAS nibabel transform into a LPS ITK transform
  """  
  RAS2LPS = np.diag([-1, -1, 1, 1])
  ITKtransform = np.dot(RAS2LPS, np.dot(transform, RAS2LPS))
  ITKcenter = np.dot(RAS2LPS,center)

  return ITKtransform, ITKcenter  
  
def write_itk_transform( transform_file, transform, center ):
  f = open(transform_file, 'w')
  f.write('#Insight Transform File V1.0\n')
  f.write('Transform: MatrixOffsetTransformBase_double_3_3\n')
  line = 'Parameters: '+np.str(transform[0,0])+' '+np.str(transform[0,1])+' '+np.str(transform[0,2])+' '
  line+= np.str(transform[1,0])+' '+np.str(transform[1,1])+' '+np.str(transform[1,2])+' '
  line+= np.str(transform[2,0])+' '+np.str(transform[2,1])+' '+np.str(transform[2,2])+' '
  line+= np.str(transform[0,3])+' '+np.str(transform[1,3])+' '+np.str(transform[2,3])+'\n'
  f.write(line)
  line = 'FixedParameters: '+np.str(center[0])+' '+np.str(center[1])+' '+np.str(center[2])+'\n'
  f.write(line)
  f.close()