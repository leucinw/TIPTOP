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

''' assign the BEST torsion parameters from database '''

import sys
import argparse
import numpy as np
from pybel import *

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'txyz') 
  parser.add_argument('-k', dest = 'key')
  parser.add_argument('-m', dest = 'mode', type=int)
  args = vars(parser.parse_args())
  txyz = args["txyz"] 
  key = args["key"] 
  mode = args["mode"] 
  
  numbers, types = np.loadtxt(txyz, usecols=(0,5), dtype="str", unpack=True, skiprows=1)
  number_type = dict(zip(numbers, types))
  
  type_class_dict = {}
  for line in open(key).readlines():
    d = line.split()
    if "atom " in line:
      if d[1] in types:
        type_class_dict[d[1]] = d[2]
  lines = open(txyz).readlines()
  if len(lines[0].split()) == 1:
    with open(txyz, "w") as f:
      f.write(lines[0].split("\n")[0] + " comments\n")
      for line in lines[1:]: 
        f.write(line)
  rootdir = os.path.join(os.path.split(__file__)[0])
  torsionfile = os.path.join(rootdir, "database.ParmGen", "torsion.dat")
  torsionsmarts = []
  for line in open(torsionfile).readlines():
    if ("#" not in line) and (len(line) > 10):
      torsionsmarts.append(line.split()[0])
  
  tor_prm_key =  {}
  lines = open(key).readlines()
  idx = 0
  for line in lines:
    if ("ASSIGN" in line):
      idx = lines.index(line) 
  for i in range(idx, len(lines)):
    if ("torsion" in lines[i]) and ("#" not in lines[i]):
      d = lines[i].split()
      tor = '  '.join(d[1:5])
      prm = '  '.join(d[5:])
      tor_prm_key[tor] = prm
  
  if mode == 1:
    tor_prm_database =  {}
    lines = open(torsionfile).readlines()
    for line in lines:
      if ("#" not in line) and (len(line) > 10):
        d = line.split()
        smt = d[0] 
        prm = '  '.join(d[1:])
        tor_prm_database[smt] = prm
  
  tmp = []
  for mol in readfile('txyz', txyz):
    for smt in torsionsmarts:
      smarts = Smarts(smt)
      results = smarts.findall(mol)
      if results:
        for r in results:
          types  = [type_class_dict[number_type[str(r[i])]] for i in range(4)]
          typesr = [type_class_dict[number_type[str(r[i])]] for i in range(3,-1,-1)]
          tor  = '  '.join(types)
          torr = '  '.join(typesr)
  
          # learning mode
          if mode == 0:
            if (tor in tor_prm_key) and (smt not in tmp):
              print(f"{smt:<40s}", tor_prm_key[tor])
              tmp.append(smt)
            if (torr in tor_prm_key) and (smt not in tmp):
              print(f"{smt:<40s}", tor_prm_key[torr])
              tmp.append(smt)
  
          # working mode
          elif mode == 1:
            if (smt in tor_prm_database) and (tor in tor_prm_key) and (tor not in tmp):
              print(f"torsion {tor} {tor_prm_database[smt]}")
              tmp.append(tor)
            if (smt in tor_prm_database) and (torr in tor_prm_key) and (torr not in tmp): 
              print(f"torsion {torr} {tor_prm_database[smt]}")
              tmp.append(torr)
           
          else:
            sys.exit(f"Error: mode {mode} is not supported")
  return

if __name__ == "__main__":
  main()
