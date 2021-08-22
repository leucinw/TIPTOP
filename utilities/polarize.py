''' 
MIT License

Copyright (c) 2021 Chengwen Liu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import numpy as np
from pybel import *

def getpolar(xyz, databasedir):
  polarDict = []
  fname, _ = os.path.splitext(xyz)
  for mol in readfile('xyz',xyz):
    polarDict = {}
    commentsDict = {} 
    natoms = len(mol.atoms)
    atoms = list(range(1, natoms+1)) 
    polarDict = dict.fromkeys(atoms, 0)
    lines = open(os.path.join(databasedir, "amoebaplus_polar.prm")).readlines()
    for line in lines:
      d = line.split()
      smt = d[0]
      pol = d[2]
      smarts = Smarts(smt)
      match = smarts.findall(mol)
      if match:
        for i in range(len(match)):	
          polarDict[match[i][0]] = pol
  return polarDict 
