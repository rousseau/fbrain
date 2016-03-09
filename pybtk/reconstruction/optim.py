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
from scipy.optimize import minimize
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
  
def optimizebis(H,x,y,maxiter,method):
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
  print 'Optimization method : '+method

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
  
  if method=='L-BFGS-B':
    maxfun  = 1000000 #10000000
    #maxcor = 10
    #minimize_options={'maxiter': maxiter,'ftol': 0,'gtol': 0, 'maxfun': maxfun, 'disp':True, 'maxcor': maxcor}
    minimize_options={'maxiter': maxiter,'ftol': 0,'gtol': 0, 'maxfun': maxfun, 'disp':True}
    

    #TODO : modify the boundary of max intensity using all input LR images (instead of just y[0])
    result = minimize(lossL2,x, method='L-BFGS-B', bounds=[(0, np.max(y[0]))] * x.size, options=minimize_options)  
  if method=='CG':
    result = minimize(lossL2,x, method='CG', options=minimize_options)  
  if method=='Nelder-Mead':
    result = minimize(lossL2,x, method='Nelder-Mead', options=minimize_options)  
  if method=='Powell':
    result = minimize(lossL2,x, method='Powell', options=minimize_options)  
    
    
  return result  