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

def main():
  template = sys.argv[1]
  fname = sys.argv[2]
  rootdir = os.path.join(os.path.split(__file__)[0])
  cmd = f"babel -ipdb {fname}.pdb -oxyz  {fname}_H.xyz -h"
  os.system(cmd)
  cmd = f"babel -ixyz {fname}.xyz -otxyz {fname}.txyz"
  os.system(cmd)
  cmd = f"python {rootdir}/IP_MatchTXYZ.py {template} {fname}_H.xyz"
  os.system(cmd)

  # here the idea is to add H back and run IP_MatchTXYZ 
  # added Hs are always appended in the end
  # assign types for the needed fragment
  lines = open(f"{fname}.txyz").readlines()
  types = np.loadtxt(f"{fname}_H.txyz_2", usecols=(5,), dtype='str', unpack=True, skiprows=1)[:len(lines)]
  
  with open(fname + ".txyz_2", 'w') as f:
    f.write(lines[0])
    for i in range(1, len(lines)):
      d = lines[i].split()
      line = ' '.join(d[:5] + [types[i-1]] + d[6:]) + "\n"
      f.write(line)
  return

if __name__ == "__main__":
  main()
