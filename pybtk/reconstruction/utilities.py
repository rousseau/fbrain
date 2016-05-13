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

