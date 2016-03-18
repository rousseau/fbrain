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

def boxcar(x,threshold):
  if np.abs(x) < threshold:
    return 1.0
  else:
    return 0.0
      
def gaussian(x,sigma):
  return np.exp(-x*x/(2*sigma*sigma))  
    
def compute_psf(lowResolution, highResolution, psftype='boxcar'):
  print 'Computing PSF'
  print 'Type of PSF : '+psftype
  #default case : boxcar
  #Appropriate size when PSF is a boxcar
  psfShape = (np.int16(np.ceil(lowResolution / highResolution)))
  
  sigma = np.zeros(3)
  if psftype == 'gauss':
    #size for gaussian kernel
    truncate = 2 #standard deviation cut : 1 = 68%, 2 = 95%, 3 : 99%
    #Sigma size from Kainz et al. TMI 2015   
    sigma[0] = 1.2 * lowResolution[0] / 2.3548
    sigma[1] = 1.2 * lowResolution[1] / 2.3548
    sigma[2] = lowResolution[2] / 2.3548
    psfShape = 2 * (truncate * sigma / highResolution + 0.5).astype(int) + 1
  if psftype == 'boxcar':
    sigma[0:3] = lowResolution[0:3]/2

    
  HRpsf = np.zeros(psfShape)
  
  print 'Shape of PSF : '
  print HRpsf.shape
  #center of HR psf
  cx = (HRpsf.shape[0]-1.0) / 2.0
  cy = (HRpsf.shape[1]-1.0) / 2.0
  cz = (HRpsf.shape[2]-1.0) / 2.0
  
  #the oversampling rate is used to provide a more accurate estimate of the PSF value
  oversampling = 25
  
  for x in range(HRpsf.shape[0]):
    for y in range(HRpsf.shape[1]):
      for z in range(HRpsf.shape[2]):
        
        val = 0
        
        #perform oversampling
        for xs in range(oversampling):
          for ys in range(oversampling):
            for zs in range(oversampling):
              
              xpsf = x - cx - 0.5 + 1.0 * xs / oversampling + 1.0 / (2*oversampling)
              ypsf = y - cy - 0.5 + 1.0 * ys / oversampling + 1.0 / (2*oversampling)
              zpsf = z - cz - 0.5 + 1.0 * zs / oversampling + 1.0 / (2*oversampling)
               
              #convert (xpsf,ypsf,zpsf) into world space (in mm)
              xlr = xpsf*highResolution[0]
              ylr = ypsf*highResolution[1]
              zlr = zpsf*highResolution[2]

              if psftype == 'gauss':
                val += gaussian(xlr,sigma[0]) * gaussian(ylr, sigma[1]) * gaussian(zlr, sigma[2])
              else:
                val += boxcar(xlr,sigma[0])*boxcar(ylr,sigma[1])*boxcar(zlr,sigma[2])
              
              
        HRpsf[x,y,z] = val/(oversampling*oversampling*oversampling)  
  HRpsf /= np.sum(HRpsf)      
  print HRpsf
  return HRpsf  