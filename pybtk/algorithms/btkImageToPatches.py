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
from sklearn.externals import joblib

#print path.dirname( path.dirname( path.abspath(__file__) ) )
#Add path to pybtk:
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )

import pybtk.filters.imagefilters as filters

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--ps', help='Patch size (for sparse coding)', type=int, default=3)
  parser.add_argument('--sub', help='Subsampling rate for patch-based sparse coding', type=float, default = 1.0)
  parser.add_argument('--padding', help='Padding value used when no mask is provided', type=float, default=0)

  parser.add_argument('-i', '--input', help='Input image filename', type=str, required=True)
  parser.add_argument('-m', '--mask',  help='Input mask filename ', type=str)
  parser.add_argument('-o', '--output',help='Output patches to save into pickle file', type=str, required=True)

  args = parser.parse_args()

  print("Loading input image")
  image = nibabel.load(args.input)

  mask = None  
  if args.mask is not None:
    mask = nibabel.load(args.mask)
  else:
    mask = filters.create_mask_image_using_padding_value(image,args.padding)
  
  #compressed means : remove 0 patches (or out of mask) and do possibly subsampling  
  p = filters.image_to_compressed_patches(image,mask,args.ps,args.sub)  
  
  #Save into pickle file
  #Issue related to files larger than 4Go on Mac -> workaround : use subsampling
  joblib.dump(p, args.output, compress=1)
    