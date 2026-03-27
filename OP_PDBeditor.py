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
import argparse
import numpy as np

# color
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\033[93m'

def trmPDB(inp,otp):
  dna_libs = ['DC', 'DG', 'DT', 'DA', 'DC3', 'DG3', 'DT3', 'DA3', 'DC5', 'DG5', 'DT5', 'DA5']
  rna_libs = ['C', 'G', 'U', 'A', 'C3', 'G3', 'U3', 'A3', 'C5', 'G5', 'U5', 'A5']
  protein_libs = ['ALA', 'GLY',]
  water_libs = ['TIP3', 'WAT']
  ion_libs = ['CLA', 'POT', 'SOD', 'NA', 'CL', 'K']
  residues = []
  res2lib = {'dna':dna_libs, 'rna':rna_libs, 'protein':protein_libs, 'water':water_libs, 'ion':ion_libs}
  for s in sel:
    if s in list(res2lib.keys()):
      residues += res2lib[s] 
    else:
      residues += [s]
  print(YELLOW + "I will check the existence of " + ' '.join(residues) + ENDC) 
  with open(otp, 'w') as fo:
    iframe = 0 
    f = open(inp, 'r')
    while True:
      line = f.readline()
      if not line:
        break
      if "MODEL " in line:
        iframe += 1
      if iframe%interval == 1:
        if line[17:21].strip() in residues:
          fo.write(line)
        if ("MODEL " in line):
          print(GREEN + f"Writing frame {iframe} ... " + ENDC)
          fo.write(f"MODEL{1 + int(iframe/interval):>6d}\n")
        if ("ENDMDL" in line):
          fo.write(line)
    f.close()
  return 

def ctlPDB(inp, opt):
  massdict = {"H": 1.007941, "C":12.01074, "N":14.006703, "O":15.999405, "S":32.0648, "P":30.973761998, \
              "F": 18.99840, "CL":35.4529, "BR":79.9035, "NA":22.989769, "K":39.0983, "MG":24.3051, "ZN":65.38}
    
  return 

def avgPDB(inp, opt):
  # readin the first frame pdb
  pdblines = []
  f = open(inp, 'r')
  iline = 0
  atom1 = 'somenonblankstrings'
  while True:
    line = f.readline()
    if (atom1 in line) and (iline > 1):
      break
    if ("ATOM " in line) or ("HETATM " in line):
      iline += 1
      pdblines.append(line)
    if iline == 1:
      atom1 = line[:27]
  f.close()

  # average the coordinate
  idxs = []
  xyzs = []
  for i in range(len(pdblines)):
    idxs.append(i)
    xyzs.append([0.0, 0.0, 0.0])
  coords = dict(zip(idxs, xyzs))
  f = open(inp, 'r')
  while True:
    line = f.readline()
    if not line:
      break
    if atom1 in line:
      iline = 0
    if ("ATOM " in line) or ("HETATM " in line):
      x1 = float(line[30:38])/nframe
      y1 = float(line[38:46])/nframe
      z1 = float(line[46:54])/nframe
      x0, y0, z0 = coords[iline]
      coords[iline] = [x0+x1, y0+y1, z0+z1]
      iline += 1
  f.close()
  
  # write the averaged PDB
  with open(opt, 'w') as f:
    for i in range(len(pdblines)):
      line = pdblines[i]
      coord = coords[i]
      xyz = "%8.3f%8.3f%8.3f"%(coord[0], coord[1], coord[2])
      f.write(line[:30] + xyz + line[54:])
  return

def cmbPDB(inps, opt):
  with open(opt, 'w') as f:
    iframe = 0
    for inp in inps:
      print(GREEN + f"Dealing with {inp} ..." + ENDC)
      fi = open(inp, 'r')
      while True:
        line = fi.readline()
        if not line:
          break
        if "MODEL " in line:
          iframe += 1
          f.write(f"MODEL{iframe:>6d}\n")
        else:
          f.write(line)
      fi.close()
  return

if __name__ == "__main__":
  actions = ["trim", "average"]
  global sel, interval, nframe 
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', dest = 'ipdbs', nargs='+', help = "input pdb file(s)", required=True)  
  parser.add_argument('-o', dest = 'opdbs', nargs='+', help = "output pdb file(s)", required=True)  
  parser.add_argument('-a', dest = 'act', help = "actions to take. Choose from [trim, average, centralize, combine]", type=str.lower, required=True)  
  parser.add_argument('-s', dest = 'sel', nargs='+', help = "selections. Can be [protein, dna, rna, water, ion] or given residue name", type=str.lower, default=[])  
  parser.add_argument('-f', dest = 'interval', help = "interval between frame. Default is 1", default=1, type=int)  
  parser.add_argument('-n', dest = 'nframe', help = "number of frames. Default is 1", default=1, type=int)  
  args = vars(parser.parse_args())
  ipdbs = args["ipdbs"]
  opdbs = args["opdbs"]
  act = args["act"]
  sel = args["sel"]
  interval = args["interval"]
  nframe = args["nframe"] 

  nipdb = len(ipdbs)
  nopdb = len(opdbs)
  if act in ['trim', 'average']:
    if nipdb != nopdb:
      sys.exit(RED + f"Error: Numbers of Input ({nipdb}) and Output ({nopdb}) Differ" + ENDC)
    else:
      if act in ['trim',]:
        for ipdb, opdb in zip(ipdbs, opdbs):
          trmPDB(ipdb, opdb)
      elif act in ['centralize',]:
        for ipdb, opdb in zip(ipdbs, opdbs):
          ctlPDB(ipdb, opdb)
      elif act in ['average']:
        for ipdb, opdb in zip(ipdbs, opdbs):
          if nframe == 1:
            sys.exit(RED + "Error: please provide the number of frames for averaging" + ENDC)
          else:
            print(YELLOW + f"Caution: Averaging Over {nframe} Frames" + ENDC)
            avgPDB(ipdb, opdb)
  elif act in ['combine']:
    if nopdb != 1:
      sys.exit(RED + f"Error: Number of Output ({nopdb}) Must be ONE for combine Action" + ENDC)
    else:
      cmbPDB(ipdbs, opdbs[0])
  else:
    sys.exit(RED + "Currently supported action: trim, average, combine" + ENDC)
