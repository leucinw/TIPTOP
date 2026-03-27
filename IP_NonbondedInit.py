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
import sys
import numpy as np
try:
  from openbabel import pybel 
except:
  import pybel
import networkx as nx

# Print general set of nonbonded parameters for AMOEBA/AMOEBA+, excluding multipole/polarize
# Given the tinker xyz file
# Assumption: atom type == atom class

def print_initial_parameters(txyz):
  alltypes = [] 
  allelements = []
  vdws = {'H':   'vdw TYPE   3.100   0.020   1.00', \
          'C':   'vdw TYPE   3.800   0.100', \
          'N':   'vdw TYPE   3.800   0.100', \
          'O':   'vdw TYPE   3.800   0.100', \
          'F':   'vdw TYPE   3.800   0.100', \
          'CL':  'vdw TYPE   3.800   0.100', \
          'BR':  'vdw TYPE   3.800   0.100', \
          'S':   'vdw TYPE   3.800   0.100', \
          'P':   'vdw TYPE   4.000   0.100', \
          'I':   'vdw TYPE   4.400   0.100', \
         }
  chgtrns= {'H':   'chgtrn TYPE   9.900   4.800', \
            'C':   'chgtrn TYPE  11.000   3.500', \
            'N':   'chgtrn TYPE  11.500   3.400', \
            'O':   'chgtrn TYPE  11.500   3.300', \
            'F':   'chgtrn TYPE  11.500   3.300', \
            'CL':  'chgtrn TYPE  11.500   3.200', \
            'BR':  'chgtrn TYPE  11.500   3.200', \
            'S':   'chgtrn TYPE  11.500   3.200', \
            'P':   'chgtrn TYPE  11.500   3.200', \
            'I':   'chgtrn TYPE  12.500   3.200', \
         }
  chgpens= {'H':   'chgpen TYPE   1.000   3.100', \
            'C':   'chgpen TYPE   6.000   3.600', \
            'N':   'chgpen TYPE   7.000   3.600', \
            'O':   'chgpen TYPE   8.000   3.600', \
            'F':   'chgpen TYPE   9.000   3.500', \
            'CL':  'chgpen TYPE  17.000   3.500', \
            'BR':  'chgpen TYPE  35.000   3.500', \
            'S':   'chgpen TYPE  16.000   3.500', \
            'P':   'chgpen TYPE  15.000   3.500', \
            'I':   'chgpen TYPE  53.000   3.500', \
         }
  lines = open(txyz).readlines()
  for line in lines[1:]:
    d = line.split()
    if d[5] not in alltypes:
      alltypes.append(d[5]) 
      allelements.append(str(d[1].upper())) 
  
  for atype, element in zip(alltypes, allelements): 
    print(chgtrns[element].replace('TYPE', str("%4d"%int(atype))))
  for atype, element in zip(alltypes, allelements): 
    print(chgpens[element].replace('TYPE', str("%4d"%int(atype))))
  for atype, element in zip(alltypes, allelements): 
    print(vdws[element].replace('TYPE', str("%4d"%int(atype))))
  return


if __name__ == "__main__":
  usage = "python lvalence_init.py txyz"
  if len(sys.argv) != 2: 
    sys.exit(f"Usage of this program: {usage}")
  txyzfile = sys.argv[1]
  if not os.path.isfile(txyzfile):
    sys.exit(f"{txyzfile} does not exist")
  print_initial_parameters(txyzfile)