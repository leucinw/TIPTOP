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

def readpdb(inp):
  atomlines = []
  for line in open(inp).readlines():
    if ("ATOM " in line) or ("HETATM" in line):
      atomlines.append(line)
  return atomlines

def readarc(inp):
  farc = open(inp, 'r')
  iline = 0
  atomlines = []
  boxlines = []
  while True:
    line = farc.readline()
    if not line:break
    if iline == 0:
      natom = int(line.split()[0])
    if pbc:
      if iline%(natom+2) == 1:
        boxlines.append(line)
        print(GREEN + f"Reading {arc} frame: {int(iline/(natom+2)) + 1}" + ENDC)
      if iline%(natom+2) > 1:
        atomlines.append(line)
    else:
      if iline%(natom+1) == 0:
        print(GREEN + f"Reading {arc} frame: {int(iline/(natom+1)) + 1}" + ENDC)
      if iline%(natom+1) > 0:
        atomlines.append(line)
    iline += 1
  farc.close()
  return atomlines, boxlines, natom

def txyzpdb(pdblines, arclines):
  # default: same atom order in pdb and arc
  index_list = [i for i in range(natom)]
  # prepare mapping.lst if atom order differs
  ## hint: use IP_MatchTXYZ.py with -s option
  if os.path.isfile("mapping.lst"):
    for line in open("mapping.lst").readlines():
      d = line.split()
      index_list[int(d[0])] = int(d[1])
  else:
    print(YELLOW + "Warning: assuming your pdb and arc have same atom order" + ENDC)
  
  # write a new pdb
  with open("new.pdb", 'w') as f:
    nframe = 1
    for i in range(len(arclines)):
      j = i%natom
      k = int(i/natom)*natom
      if j == 0:
        f.write(f"MODEL{nframe:>6d}\n")
        print(GREEN + f"Writing frame {nframe}" + ENDC)
        nframe += 1
      d = arclines[k+index_list[j]].split()
      x,y,z = float(d[2]), float(d[3]), float(d[4])
      xyzstr = "%8.3f%8.3f%8.3f"%(x,y,z)
      newline = pdblines[j][:30] + xyzstr + pdblines[j][54:]
      f.write(newline)
      if j == (natom-1):
        f.write("ENDMDL\n")
  return True 

if __name__ == "__main__":
  global pbc,pdb,arc,natom
  parser = argparse.ArgumentParser()
  parser.add_argument('-pdb', dest = 'pdb', help = "Template pdb file", required=True)  
  parser.add_argument('-arc', dest = 'arc', help = "Tinker trajectory file", required=True)  
  parser.add_argument('-pbc', dest = 'pbc', help = "PBC in arc. True/False", default=True, type=bool)  
  args = vars(parser.parse_args())
  pdb = args["pdb"]
  arc = args["arc"]
  pbc = args["pbc"]
  
  if not os.path.isfile(arc):
    sys.exit(RED + f"Error: {arc} does not exist" + ENDC)
  if not os.path.isfile(pdb):
    sys.exit(RED + f"Error: {pdb} does not exist" + ENDC)
 
  pdblines = readpdb(pdb)
  arclines, boxlines, natom = readarc(arc)
  # no. of natoms 
  natom1 = len(pdblines)
  if natom != natom1:
    sys.exit(RED + "Number of atoms in PDB and ARC not the same" + ENDC)
  signal = txyzpdb(pdblines, arclines) 
  if signal:
    print(GREEN + f"Your {arc} has been converted to new.pdb" + ENDC)
  else:
    print(RED + f"Sorry, I could not convert your {arc} to pdb" + ENDC)
