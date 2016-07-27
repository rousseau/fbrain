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

import argparse
import nibabel
import numpy as np
from scipy.stats import rice

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--input', help='Input image', type=str, required = True)
  parser.add_argument('-o','--output', help='Output image', type=str, required = True)
  parser.add_argument('-m','--mean', help='Mean of added gaussian noise', type=float, default = 0.0)
  parser.add_argument('-s','--std', help='Standard deviation of added gaussian noise', type=float, default = 1.0)
  parser.add_argument('-r','--rician', help='Rician noise is added', type=bool, default=False)
  parser.add_argument('-b',help='b parameter for Rician distribution', type=float, default = 1.0)
  parser.add_argument('-n','--normalize', help='Normalized data between [0,1]', type=bool, default=False)
  

  args = parser.parse_args()
  
  # Create random dictionary
  np.random.seed(12345)

  input_image = nibabel.load(args.input)
  data = input_image.get_data().astype(float)
  
  maxval = 0
  if args.normalize==True:
    maxval = np.max(data)
    data /= maxval
    
  if args.rician==False:
    data += np.random.normal(args.mean, args.std, data.shape)
  else:
    data += rice.rvs(args.b, loc=args.mean, scale=args.std, size=data.shape)
    

  if args.normalize==True:
    data *= maxval
    
  nibabel.save(nibabel.Nifti1Image(data,input_image.affine),args.output)  
    
  """A Rice continuous random variable. (from scipy code)
    Notes
    -----
    The probability density function for `rice` is::
        rice.pdf(x, b) = x * exp(-(x**2+b**2)/2) * I[0](x*b)
    for ``x > 0``, ``b > 0``.
  """
 
  