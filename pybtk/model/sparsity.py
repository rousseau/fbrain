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
from sklearn.decomposition import MiniBatchDictionaryLearning,SparseCoder
from sklearn.externals import joblib
import matplotlib.pyplot as plt

from os import path
import sys
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )
import pybtk.model.patches as mp

def scskl_dico_learning(list_pickled_array,n_atoms,maxepoch=5,maxiter=100):
  D = None
  for e in range(maxepoch):
    for a in list_pickled_array:
      data = joblib.load(a)
      dico = MiniBatchDictionaryLearning(n_components=n_atoms, n_iter=maxiter, dict_init=D)
      D = dico.fit(data).components_.astype(np.float32)
  return D      
      
      
def scskl_reconstruction(data,mask,D):
  output = np.zeros(data.shape)
  fmap = np.zeros((D.shape[0]))
  #fdata = np.zeros((data.shape[0],data.shape[1],data.shape[2],D.shape[0]))

  px = np.int(np.around(np.power(D.shape[1],1/3))) #patch size (is assumed to be isotropic)
  hpx = np.floor(px/2).astype(int)
  nblock = 2 # number of block per dimension
  subsize = np.ceil(np.array(data.shape) / nblock).astype(int)
  
  med = np.median(data)
  currentblock = 1
  
  for x in range(np.ceil(data.shape[0]/subsize[0]).astype(int)):
    xmin = x*subsize[0]
    xmax = np.min((data.shape[0],(x+1)*subsize[0]))
    for y in range(np.ceil(data.shape[1]/subsize[1]).astype(int)):
      ymin = y*subsize[1]
      ymax = np.min((data.shape[1],(y+1)*subsize[1]))
      for z in range(np.ceil(data.shape[2]/subsize[2]).astype(int)):  
        zmin = z*subsize[2]
        zmax = np.min((data.shape[2],(z+1)*subsize[2]))
        
        print('Processing block : ',currentblock)
        currentblock+=1
        
        #Enlarge subimage to take into account block effect due to non-overlapping patches
        xmin2 = np.max((0,xmin-hpx))
        xmax2 = np.min((data.shape[0],xmax+hpx))
        ymin2 = np.max((0,ymin-hpx))
        ymax2 = np.min((data.shape[1],ymax+hpx))
        zmin2 = np.max((0,zmin-hpx))
        zmax2 = np.min((data.shape[2],zmax+hpx))
      
        subdata = data[xmin2:xmax2,ymin2:ymax2,zmin2:zmax2]
        submask = mask[xmin2:xmax2,ymin2:ymax2,zmin2:zmax2]
        p = mp.array_to_patches(subdata,patch_shape=(px,px,px),normalization=False)
        pm = mp.array_to_patches(submask,patch_shape=(px,px,px),normalization=False)
        #remove patch we dont want to process
        index = ~np.all(pm==0,axis=1)
        subp = p[index]
        subp -= med
        
        if subp.shape[0] > 0:
          print('Number of patches to process: ',subp.shape[0])
          #Currently, there is a bug when using n_jobs>1 (https://github.com/scikit-learn/scikit-learn/issues/5956)
          coder = SparseCoder(dictionary=D, transform_algorithm='omp')
          code = coder.transform(subp).astype(np.float32)
          fmap += np.sum((np.fabs(code)>0),axis=0)
          subp = np.dot(code, D)          
          subp += med 
          p[index] = subp 
          suboutput = mp.patches_to_array(patches=p, patch_shape=(px,px,px), array_shape=subdata.shape)
          
          tmpoutput = np.empty(data.shape)
          tmpoutput[xmin2:xmax2,ymin2:ymax2,zmin2:zmax2]= suboutput      
          output[xmin:xmax,ymin:ymax,zmin:zmax] = tmpoutput[xmin:xmax,ymin:ymax,zmin:zmax]
          
#          for a in range(D.shape[0]):
#            for s in range(subp.shape[0]):
#              subp[s,:] = code[s,a]
#            p.fill(0)
#            p[index] = subp
#            fa = mp.patches_to_array(patches=p, patch_shape=(px,px,px), array_shape=subdata.shape)
#            to = np.empty(data.shape)
#            to[xmin2:xmax2,ymin2:ymax2,zmin2:zmax2]= fa      
#            fdata[xmin:xmax,ymin:ymax,zmin:zmax,a] = to[xmin:xmax,ymin:ymax,zmin:zmax]
              
            
          
  #plt.bar(range(0,D.shape[0]), fmap)
  #print('points in mask: ',np.sum(mask!=0))
  #print('Number of non zero elements: ',np.sum(fmap)/np.sum(mask!=0))
  #plt.show()
      
#  return (output,fdata)
  return output
    
