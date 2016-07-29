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

import pybtk.model.sparsity as sparsity
import pybtk.filters.imagefilters as filters


if __name__ == '__main__':

  parser = argparse.ArgumentParser()

  parser.add_argument('-i', '--input',  help='Input image pickle filename ', type=str, action='append')
  parser.add_argument('-o', '--output', help='Output dictionary filename', type=str)
  parser.add_argument('-a', '--natoms', help='Number of atoms', type=int, default=100)
  parser.add_argument('-e', '--nepoch', help='Number of epochs', type=int, default=5)
  parser.add_argument('-m', '--miter' , help='Number of max iterations', type=int, default=100)

  args = parser.parse_args()
  
  D = sparsity.scskl_dico_learning(args.input,args.natoms,maxepoch=args.nepoch,maxiter=args.miter)
  joblib.dump(D, args.output, compress=1)
  