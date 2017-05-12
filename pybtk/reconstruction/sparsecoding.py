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

from sklearn.decomposition import MiniBatchDictionaryLearning
from sklearn.decomposition import SparseCoder
import multiprocessing

def parallel_sc(args):
    
  (dico,p)=args
  coder = SparseCoder(dictionary=dico, transform_algorithm='omp')
  #by default, number of non zero coefficients is 0.1 * n_features
  #Zeyde et al.: 3 non zero coefficients !
  code = coder.transform(p).astype(np.float32)    
  return code 

class SparseCoding():
  def __init__(self, n_atoms=100, alpha=0.1, batch_size_learning=32, batch_size_coding=524288, n_iter=1000, n_jobs=32):
    self.n_atoms = n_atoms
    self.alpha = alpha
    self.batch_size_learning = batch_size_learning
    self.batch_size_coding = batch_size_coding
    self.n_iter = n_iter
    self.Dhr = None
    self.Dlr = None
    self.n_jobs = n_jobs
    
  def fit(self, Xhr, Xlr):
    X = np.concatenate((Xhr,Xlr),axis=1)
    dico = MiniBatchDictionaryLearning(n_components=self.n_atoms, n_iter=self.n_iter, alpha = self.alpha, batch_size=self.batch_size_learning, verbose=False)
    D = dico.fit(X).components_.astype(np.float32)
    n_features = Xhr.shape[1]
    self.Dhr = np.copy(D[:,0:n_features])
    self.Dlr = np.copy(D[:,n_features:2*n_features])
  
  def transform(self, X):
    Y = np.zeros(X.shape)
    #Split the data wrt to the max batch size for coding
    n_samples = X.shape[0]
    n_split = np.int(n_samples/self.batch_size_coding)
    Xsplit = np.array_split(X,n_split)

    p_index = 0
    #We split again the data to limit the use of RAM (however, it's slower than direct computation)
    for samples in Xsplit:
      print('{0:.0f}'.format(p_index*100.0/X.shape[0])),

      tmpsplit = np.array_split(samples,np.int(samples.shape[0]/self.n_jobs))
      params = [None]*len(tmpsplit)
      for i in range(len(tmpsplit)):
        params[i] = [self.Dlr,tmpsplit[i]]

      print(' .'),
      pool_sc = multiprocessing.Pool(processes=self.n_jobs)
      code = pool_sc.map(parallel_sc,params)
      pool_sc.close()
      pool_sc.join()
      print('.'),
      for c in code:
        Y[p_index:p_index+c.shape[0],:] = np.dot(c, self.Dhr)
        p_index += c.shape[0]
      print('.'),
    return Y

