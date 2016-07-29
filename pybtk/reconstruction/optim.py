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
from scipy.optimize import minimize, approx_fprime,line_search
from time import time
import nibabel
from os import path
import sys
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )
from pybtk.reconstruction.utilities import convert_image_to_vector


def lossL2(x,H,y):
  """
  Loss L2 : Hx-y
  
  Parameters
  ----------
  H : list of sparse matrix
  y : list of low resolution images (array)
  x : high resolution (array)
  
  Returns
  -------  
  res : loss    
  """      
  
  res = 0
  for i in range(len(y)):
    tmp = H[i].dot(x) - y[i]
  
    res += np.linalg.norm(tmp) 
  return res

def lossL2prime(x,H,y):
  """
  Analytical computation of the gradient of loss L2 : Hx-y
  
  Parameters
  ----------
  H : list of sparse matrix
  y : list of low resolution images (array)
  x : high resolution (array)
  
  Returns
  -------  
  res : gradient    
  """      
  #It could be done also using approx_fprime
  #epsilon = np.ones(x.shape) * 0.01 * (np.max(x)-np.min(x));
  #grad = approx_fprime(x, f, epsilon,H,y)

  res = np.zeros(x.shape)
  #analytic mean gradient over each low resolution image y
  for i in range(len(y)):
    tmp = H[i].dot(x) - y[i]
    res += 2.0 * H[i].transpose().dot(tmp)
  res /= len(y)  
  return res    

def LossSpectrumprime(x,magRef,index):
  im = np.zeros(magRef.shape).reshape(-1)
  im[index] = x
  im = im.reshape(magRef.shape)
  # ax = np.hanning(magRef.shape[0]).reshape((-1,1,1))
  # ay = np.hanning(magRef.shape[1]).reshape((1,-1,1))
  # az = np.hanning(magRef.shape[2]).reshape((1,1,-1))
  # window = ax*ay*az
  # window = np.ones(I.shape)
  # im = im*window

  epsilon = 10**-10
  fftIm = np.fft.fftn(im)
  magIm = np.abs(fftIm + epsilon)
  newIm = np.real( np.fft.ifftn( magRef / magIm * fftIm) )
  return (im-newIm).reshape((-1))[index]
  
def optimize_L_BFGS_B(H,x,y, maxiter):
  """
  Optimization for super resolution using L-BFGS-B method
  
  Parameters
  ----------
  H : list of sparse matrix
  y : list of low resolution images (array)
  x : high resolution (array)
  
  Returns
  -------  
  result : results of optimization
    result.x contains the estimation of x
    
  """      
  maxfun  = 1000000 #10000000
  maxcor = 10
  minimize_options={'maxiter': maxiter,'ftol': 0,'gtol': 0, 'maxfun': maxfun, 'disp':True, 'maxcor': maxcor}

  result = minimize(lossL2,x, args=(H,y), method='L-BFGS-B', bounds=[(0, np.max(y))] * x.size, options=minimize_options)  
  print(result.message)
  return result.x

def computeAlpha(x,grad):
  #Simple rule to define alpha 
  if np.max(np.abs(grad)) > 0:
    alpha = 0.05 * (np.max(x)-np.min(x))/np.max(np.abs(grad))
  else:
    alpha = 0
  return alpha

def optimization(HList,x,yList,maxiter,initHRImage):
  #prepare the data for optimization 

  #compress x and H : to perform the optimizaton only on nonzero points (i.e reduce the dimension of the parameter vector)
  index = np.nonzero(x)[0]
  xc = x[index]
  HListc = []
  for i in range(len(HList)):
    HListc.append(HList[i][:,index])

  from skimage.restoration import denoise_tv_chambolle
  HRDenoised = denoise_tv_chambolle(initHRImage.get_data(),weight=10)  
  
  nibabel.save(nibabel.Nifti1Image(HRDenoised,initHRImage.affine),'cham.nii.gz')
  
  xRef = convert_image_to_vector(nibabel.Nifti1Image(HRDenoised, initHRImage.affine))[index]
  
  res = optimize(HListc,xc,yList,maxiter,index,xRef,lambdaL2=0.6)
  #decompress res.x
  x[index] = res  
  return x.reshape(initHRImage.get_data().shape)

    
def optimize(H,x,y,maxiter,index,xRef,lambdaL2=0.5):
  print('Doing super-resolution optimization')
  t = time()
  miny = np.min(y[0])
  maxy = np.max(y[0])
  print('bounds of y : '+str(miny)+', '+str(maxy))
  
  iteration = 0
  maxdiff = np.ones(len(y)) * (maxy-miny)
  threshold = 0.01 * (maxy-miny)
  
  while iteration<maxiter and np.max(maxdiff) > threshold:

    gradL2 = lossL2prime(x,H,y)
      
    #Find alpha that satisfies strong Wolfe conditions.
    #http://scipy.github.io/devdocs/generated/scipy.optimize.line_search.html#scipy.optimize.line_search
    res = line_search(lossL2, lossL2prime, x, -gradL2, args=(H,y))
    alphaL2 = res[0]
    if alphaL2 is None:
      alphaL2 = computeAlpha(x,gradL2) 
    
    update = alphaL2*gradL2    
          
    if xRef is not None:
      gradDenoising =2.0* (x-xRef)
      alphaDenoising = computeAlpha(x,gradDenoising)
      update = (1-lambdaL2)*alphaDenoising*gradDenoising + lambdaL2*alphaL2*gradL2   
       
    #Update high resolution image  
    x = x - update

    #Threshold on Maxdiff or update magnitude ?    
    for i in range(len(y)):
      maxdiff[i] = np.max(H[i].dot(x) - y[i])
    
    #Use bounds to limit intensity range of x
    x[x<miny] = miny
    x[x>maxy] = maxy
    
    iteration+=1
    if iteration==maxiter:
      print('Maximum number of iterations is reached')
  
  print('Optimization done in '+str(time()-t)+' s, in '+str(iteration)+' iterations')
  
  return x