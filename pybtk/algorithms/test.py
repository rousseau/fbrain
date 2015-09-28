# -*- coding: utf-8 -*-

import sys
from os import path

print path.dirname( path.dirname( path.abspath(__file__) ) )
sys.path.append( path.dirname( path.dirname( path.dirname( path.abspath(__file__) ) ) ) )

from pybtk.io.itk_transforms import read_itk_transform
