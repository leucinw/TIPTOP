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

''' convert PDB structure to TINKER xyz 
    1. Amino acids and nucleotides are from amoebabio18.prm
    2. Glycans are generated using Poltype
    3. Lipids are taken from literature and generated using Poltype
    4. Can add more support by adding `residue` in database directory

    Usage: python IP_PDB2txyz.py your.pdb mode
           ~ mode can be chosen from [A,B,C]
           ~ suggest run a,b,c sequencially, to make sure every step is good
           ~ can also run ABC at the same time if you are confident
           !!! your.pdb MUST have hydrogens and protonation state correctly!!!
'''

import os
import sys
import numpy as np
import concurrent.futures 

''' write a pdb for obabel '''
def prepare(pdb):
  lines = open(pdb).readlines()
  if lines[0][77] == ' ':
    with open(pdb, "w") as f:
      for line in lines:
        line = line[:76] + ' '+ line[12:16].split()[0][0] + "\n"
        f.write(line)
    os.system(f"babel {pdb} {pdb.replace('pdb', 'txyz')}")
  return

''' split the pdb into pdb/xyz'''
def splitpdb(pdb):
  resnames = []
  resids = []
  xyzs = []
  lines = open(pdb).readlines()
  curresid = lines[0][22:26]
  curresnm = lines[0][17:21]
  xyzstr = []
  pdbstr = []
  number = 0
  for line in lines:
    if (line[22:26] == curresid) and (curresnm == line[17:21]):
      atom = line.split()[-1] 
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      xyzstr.append(f"{atom:>6s}{x:12.6f}{y:12.6f}{z:12.6f}\n")
      pdbstr.append(line)
    else:
      xyzname = "%04d"%number + f"_{curresnm.split()[0]}.xyz"
      pdbname = "%04d"%number + f"_{curresnm.split()[0]}.pdb"
      xyzs.append(xyzname)
      with open(xyzname, 'w') as f:
        f.write("%6d"%len(xyzstr) + "\n\n")
        for s in xyzstr:
          f.write(s)
      with open(pdbname, 'w') as f:
        for s in pdbstr:
          f.write(s)
      number += 1
      curresid = line[22:26]
      curresnm = line[17:21]
      atom = line.split()[-1] 
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      xyzstr = [(f"{atom:>6s}{x:12.6f}{y:12.6f}{z:12.6f}\n")]
      pdbstr = [line,]
  
  xyzname = "%04d"%number + f"_{curresnm.split()[0]}.xyz"
  pdbname = "%04d"%number + f"_{curresnm.split()[0]}.pdb"
  xyzs.append(xyzname)
  with open(xyzname, 'w') as f:
    f.write("%6d"%len(xyzstr) + "\n\n")
    for s in xyzstr:
      f.write(s)
  with open(pdbname, 'w') as f:
    for s in pdbstr:
      f.write(s)
  with open('xyzlist', 'w') as f:
    for s in xyzs:
      f.write(s + '\n')
  return 

''' read the template database '''
def readdatabase():
  database = {}
  dbdir = os.path.join(rootdir, 'database') 
  ffs = os.listdir(dbdir)
  for ff in ffs:
    dirname = os.path.join(rootdir, 'database', ff)
    if os.path.isdir(dirname):
      fs = os.listdir(dirname)
      for f in fs:
        if f.endswith(".txyz"):
          resname = f.split(".txyz")[0]
          if ff not in database:
            database[ff] = [resname]
          else:
            database[ff] += [resname]
  return database

''' convert xyz to tinker xyz '''
def xyztxyz(xyz):
  found = False
  resname = xyz.split(".xyz")[0].split("_")[-1]
  
  # glycan is special
  isGlycan = False
  if resname in ['AMAN', 'BMAN', 'AFUC', 'ANE5', 'BGAL', 'BGLN', 'AGAN']:
    isGlycan = True 

  for key, value in database.items():
    if (resname in value):
      found = True
      if resname in value:
        template = os.path.join(rootdir, 'database.PDB2txyz', key, resname+".txyz")
      t = xyz.replace("xyz", "txyz_2")
      if not os.path.isfile(t):
        if isGlycan:
          os.system(f"python {rootdir}/IP_MatchTXYZ_Glycan.py {template} {xyz.split('.')[0]}")
        else:
          os.system(f"python {rootdir}/IP_MatchTXYZ.py {template} {xyz}")
  if not found:
    sys.exit(f"Error: {resname} not found in database")
  
  return 

''' generate the final txyz '''
def connect(txyz, txyzs):
  fname = txyz + "_2"
  alltypes = []
  for t in txyzs:
    types = list(np.loadtxt(t, usecols=(5), dtype="str", unpack=True, skiprows=1))
    alltypes += types
  lines = open(txyz).readlines()
  if (len(lines)-1) != len(alltypes):
    sys.exit("Error: number of atoms not the same!!!")

  with open(fname, "w") as f:
    f.write(f"{lines[0].split()[0]:>6s}" + "  Generated by TIPTOP/IP_PDB2txyz.py\n")
    for i in range(1,len(lines)):
      d = lines[i].split()
      d[5] = alltypes[i-1]
      constr = ''.join(["%10s"%s for s in d[6:]]) 
      line = f'{d[0]:>10s}{d[1]:>3s}{float(d[2]):12.4f}{float(d[3]):12.4f}{float(d[4]):12.4f}{d[5]:>5s}{constr}\n'
      f.write(line)
  return

def main():
  global database
  global rootdir
  pdb = sys.argv[1]
  mode = sys.argv[2].upper()
  rootdir = os.path.join(os.path.split(__file__)[0])
  database = readdatabase()
  txyz = pdb.replace("pdb", "txyz")
  
  if not os.path.isfile(txyz):
    os.system(f"babel {pdb} {txyz}")
  if 'A' in mode: 
    prepare(pdb)
    splitpdb(pdb)
  if 'B' in mode:
    xyzs = np.loadtxt("xyzlist", dtype='str')
    jobs = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
      results = [executor.submit(xyztxyz, xyz) for xyz in xyzs]
      for f in concurrent.futures.as_completed(results):
        jobs.append(f.result())
    txyzs = [xyz.replace('xyz', 'txyz_2') for xyz in xyzs]
    with open('txyzlist', 'w') as f:
      for s in txyzs:
        f.write(s + '\n')
  if 'C' in mode:
    txyzs = np.loadtxt("txyzlist",dtype='str')
    connect(txyz,txyzs)
  return

if __name__ == "__main__":
  main()