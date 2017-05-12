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
from sklearn.cluster import MiniBatchKMeans, Birch
from sklearn import linear_model
import numpy as np

class ClusterAndRegress():
  def __init__(self, n_clusters=1000, clustering='MiniBatchKMeans'):

    self.n_clusters=n_clusters
    if clustering == 'MiniBatchKMeans':
      self.batch_size = 10 * self.n_clusters
      self.clustering=MiniBatchKMeans(n_clusters=self.n_clusters, batch_size = self.batch_size)
    if clustering == 'Birch':
      self.clustering = Birch(n_clusters=self.n_clusters)      
    self.regressors=[None]*n_clusters
    
  def fit(self,Xlr,Xhr):
        
    y_pred = self.clustering.fit_predict(Xlr)

    for i in range(self.n_clusters):
      x = Xlr[y_pred==i,:]
      y = Xhr[y_pred==i,:]
      self.regressors[i] = linear_model.LinearRegression()
      self.regressors[i].fit(x,y)

  def transform(self,X):
    y_lr = self.clustering.predict(X)   
    p = np.zeros(X.shape)
    for c in range(self.n_clusters):
      index = y_lr==c
      p[index,:] = self.regressors[c].predict(X[index,:])
    return p