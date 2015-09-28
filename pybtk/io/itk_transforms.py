# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 17:22:10 2015

@author: rousseau
"""

import numpy as np
  
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