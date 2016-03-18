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

def optimize(H,x,y, maxiter):
  
  def lossL2(x):
    tmp = H.dot(x) - y
    
    f_val = np.linalg.norm(tmp) 
    #g_val =  2.0 * H.transpose().dot(tmp)
    #Gradient vector makes the algorithm stop to early (instability?)
    
    #return [f_val,g_val]
    return f_val
  
  maxfun  = 100 #10000000
  maxcor = 100
  minimize_options={'maxiter': maxiter,'ftol': 0,'gtol': 0, 'maxfun': maxfun, 'disp':True, 'maxcor': maxcor}

  result = minimize(lossL2,x, method='L-BFGS-B', bounds=[(0, np.max(y))] * x.size, options=minimize_options)  
  return result
  
def optimizebis(H,x,y,maxiter):
  """
  Optimization for super resolution
  
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

  def lossL2(x):
    f_val = 0
    for i in range(len(y)):
      tmp = H[i].dot(x) - y[i]
    
      f_val += np.linalg.norm(tmp) 
    #g_val =  2.0 * H.transpose().dot(tmp)
    #Gradient vector makes the algorithm stop to early (instability?)
    
    #return [f_val,g_val]
    return f_val
  minimize_options={'maxiter': maxiter}  
  
  maxfun  = 1000000 #10000000
  #maxcor = 10
  #minimize_options={'maxiter': maxiter,'ftol': 0,'gtol': 0, 'maxfun': maxfun, 'disp':True, 'maxcor': maxcor}
  minimize_options={'maxiter': maxiter,'ftol': 0,'gtol': 0, 'maxfun': maxfun, 'disp':True}
    

  #TODO : modify the boundary of max intensity using all input LR images (instead of just y[0])
  result = minimize(lossL2,x, method='L-BFGS-B', bounds=[(0, np.max(y[0]))] * x.size, options=minimize_options)  
  
  return result 
  
def myOptimization(H,x,y,maxiter):
  print 'Doing super-resolution optimization'
  t = time()
  miny = np.min(y[0])
  maxy = np.max(y[0])
  print 'bounds of y : '+str(miny)+', '+str(maxy)
  
  def f(x,H,y):
    f_val = 0
    for i in range(len(y)):
      tmp = H[i].dot(x) - y[i]
    
      f_val += np.linalg.norm(tmp) 
    return f_val
    
  def fprime(x,H,y):
    #t = time()
    epsilon = np.ones(x.shape) * 0.01 * (np.max(x)-np.min(x));
    grad = approx_fprime(x, f, epsilon,H,y)
    #print 'time to compute grad : '+str(time()-t)
    return grad
  
  def fprime2(x,H,y):
    #t = time()
    grad = np.zeros(x.shape)
    #analytic mean gradient over each low resolution image y
    for i in range(len(y)):
      tmp = H[i].dot(x) - y[i]
      grad += 2.0 * H[i].transpose().dot(tmp)
    grad /= len(y)  
    #print 'time to compute grad2 : '+str(time()-t)
    return grad    
  

  iteration = 0
  grad = np.zeros(x.shape)
  maxdiff = np.ones(len(y)) * (maxy-miny)
  threshold = 0.01 * (maxy-miny)
  
  while iteration<maxiter and np.max(maxdiff) > threshold:
    grad = fprime2(x,H,y)
      
    #Find alpha that satisfies strong Wolfe conditions.
    #http://scipy.github.io/devdocs/generated/scipy.optimize.line_search.html#scipy.optimize.line_search
    res = line_search(f, fprime2, x, -grad, args=(H,y))
    alpha = res[0]
    if alpha is None:
      #Simple rule to define alpha  
      alpha = 0.05 * (np.max(x)-np.min(x))/np.max(np.abs(grad))
      
    #Update high resolution image  
    x = x - alpha * grad
    
    for i in range(len(y)):
      maxdiff[i] = np.max(H[i].dot(x) - y[i])
    #print 'max diff : '+str(np.max(maxdiff))
    #print 'iteration : '+str(iteration+1)
    
    #Use bounds to limit intensity range of x
    x[x<miny] = miny
    x[x>maxy] = maxy
    
    iteration+=1
    if iteration==maxiter:
      print 'Maximum number of iterations is reached'
  
  print 'Optimization done in '+str(time()-t)+' s'
  
  return x, grad