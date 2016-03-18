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
  print result.message
  return result.x
    
def optimize(H,x,y,maxiter):
  print 'Doing super-resolution optimization'
  t = time()
  miny = np.min(y[0])
  maxy = np.max(y[0])
  print 'bounds of y : '+str(miny)+', '+str(maxy)
  
  iteration = 0
  grad = np.zeros(x.shape)
  maxdiff = np.ones(len(y)) * (maxy-miny)
  threshold = 0.01 * (maxy-miny)
  
  while iteration<maxiter and np.max(maxdiff) > threshold:
    grad = lossL2prime(x,H,y)
      
    #Find alpha that satisfies strong Wolfe conditions.
    #http://scipy.github.io/devdocs/generated/scipy.optimize.line_search.html#scipy.optimize.line_search
    res = line_search(lossL2, lossL2prime, x, -grad, args=(H,y))
    alpha = res[0]
    if alpha is None:
      #Simple rule to define alpha  
      alpha = 0.05 * (np.max(x)-np.min(x))/np.max(np.abs(grad))
      
    #Update high resolution image  
    x = x - alpha * grad
    
    for i in range(len(y)):
      maxdiff[i] = np.max(H[i].dot(x) - y[i])
    
    #Use bounds to limit intensity range of x
    x[x<miny] = miny
    x[x>maxy] = maxy
    
    iteration+=1
    if iteration==maxiter:
      print 'Maximum number of iterations is reached'
  
  print 'Optimization done in '+str(time()-t)+' s, in '+str(iteration)+' iterations'
  
  return x, grad