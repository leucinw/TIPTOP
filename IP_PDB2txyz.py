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
  if lines[1][77] == ' ':
    with open(pdb, "w") as f:
      for line in lines:
        line = line[:76] + ' '+ line[12:16].split()[0][0] + "\n"
        f.write(line)
    os.system(f"babel {pdb} {pdb.replace('pdb', 'txyz')}")
  return

''' split input pdb file into fragment pdbs 
    solvent and ions are ignored
'''
def splitpdb(pdb):
  resnames = []
  resids = []
  pdbs = []
  lines = open(pdb).readlines()
  curresid = lines[0][22:26]
  curresnm = lines[0][17:21]
  pdbstr = []
  pdbstrs = []
  number_res = 0
  number_atm = 0
  for line in lines:
    if ("TER" not in line) and ("END" not in line) and (curresnm not in database['solvent']):
      if (line[22:26] == curresid) and (curresnm == line[17:21]):
        atom = line.split()[-1] 
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        pdbstr.append(line)
        pdbstrs.append(line)
        number_atm += 1
      else:
        pdbname = "%04d"%number_res + f"_{curresnm.split()[0]}.pdb"
        print(f"writing {pdbname}")
        pdbs.append(pdbname)
        with open(pdbname, 'w') as f:
          for s in pdbstr:
            f.write(s)
        number_res += 1
        curresid = line[22:26]
        curresnm = line[17:21]
        pdbstr = [line,]
        pdbstrs.append(line)
        number_atm += 1
  with open('pdblist', 'w') as f:
    for s in pdbs:
      f.write(s + '\n')

  # write solvent and ions if there are
  txyzstr = []
  for line in lines:
    if ("TER" not in line):
      curresnm = line[17:21].strip()
      if curresnm in database['solvent']: 
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom = line[13:17].strip()
        if atom == 'OH2':
          line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 349 {number_atm+1} {number_atm+2}\n"
          txyzstr.append(line_s)
          number_atm += 1
        elif atom == 'H1':
          line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 350 {number_atm-1} \n"
          txyzstr.append(line_s)
          number_atm += 1
        elif atom == 'H2':
          line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 350 {number_atm-2} \n"
          number_atm += 1 
          txyzstr.append(line_s)
        elif atom in ['POT', 'K', 'K+']: 
          line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 353\n"
          txyzstr.append(line_s)
          number_atm += 1
        elif atom in ['CLA', 'Cl', 'Cl-']: 
          line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 361\n"
          txyzstr.append(line_s)
          number_atm += 1
        else:
          sys.exit(f'Could not recognize {atom}')
  txyzname = "solvent.txyz"
  natom = len(txyzstr)
  if natom > 0: 
    with open(txyzname, 'w') as f:
      for s in txyzstr:
        f.write(s)

  with open('bio.pdb', 'w') as f:
    for s in pdbstrs[:-1]:
      f.write(s)
  return 

''' read the template database '''
def readdatabase():
  database = {}
  dbdir = os.path.join(rootdir, 'database.PDB2txyz') 
  ffs = os.listdir(dbdir)
  for ff in ffs:
    dirname = os.path.join(rootdir, 'database.PDB2txyz', ff)
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
def pdbtxyz(pdb):
  resname = pdb.split(".pdb")[0].split("_")[-1]
  
  # glycan is special
  isGlycan = False
  if resname in ['AMAN', 'BMAN', 'AFUC', 'ANE5', 'BGAL', 'BGLN', 'AGAN']:
    isGlycan = True 
  
  for key, value in database.items():
    template = os.path.join(rootdir, 'database.PDB2txyz', key, resname + "*.txyz")
    if (resname in value) or (resname + "_1" in value):
      t = pdb.replace("pdb", "txyz_2")
      if not os.path.isfile(t):
        if isGlycan:
          os.system(f"python {rootdir}/IP_MatchTXYZ_Glycan.py {template} {pdb.split('.')[0]}")
        if (not isGlycan):
          os.system(f"python {rootdir}/IP_MatchTXYZ.py -t {template} -d {pdb}")
  return 

''' generate the final txyz '''
def connect(txyz, txyzs):
  fname = txyz + "_2"
  alltypes = []
  for t in txyzs:
    types = np.loadtxt(t, usecols=(5,), dtype="str", unpack=True, skiprows=1)
    if types.ndim == 0: 
      alltypes += [str(types), ]
    else:
      alltypes += list(types)
  lines = open(txyz).readlines()
  if (len(lines)-1) != len(alltypes):
    sys.exit("Error: number of atoms not the same!!!")
  
  if os.path.isfile('solvent.txyz'):
    nlines = len(open('solvent.txyz').readlines())
  else:
    nlines = 0
  with open(fname, "w") as f:
    f.write(f"{int(lines[0].split()[0]) + nlines:>10d}" + "  Generated by TIPTOP/IP_PDB2txyz.py\n")
    for i in range(1,len(lines)):
      d = lines[i].split()
      d[5] = alltypes[i-1]
      constr = ''.join(["%10s"%s for s in d[6:]]) 
      line = f'{d[0]:>10s}{d[1]:>3s}{float(d[2]):12.4f}{float(d[3]):12.4f}{float(d[4]):12.4f}{d[5]:>5s}{constr}\n'
      f.write(line)
  return

if __name__ == "__main__":
  global database
  global rootdir
  pdb = sys.argv[1]
  
  # special names in CHARMM-GUI file
  ress_in_pdbs = {"CYT":"DC ", "GUA":"DG ", "ADE":"DA ", "THY":"DT "} 
  for res, res_ in ress_in_pdbs.items():
    cmd = f"sed 's/{res}/{res_}/g' -i {pdb}"
    os.system(cmd)
  
  mode = sys.argv[2].upper()
  rootdir = os.path.join(os.path.split(__file__)[0])
  database = readdatabase()

  
  if 'A' in mode: 
    prepare(pdb)
    splitpdb(pdb)
  if 'B' in mode:
    pdbs = np.loadtxt("pdblist", dtype='str')
    jobs = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
      results = [executor.submit(pdbtxyz, pdb) for pdb in pdbs]
      for f in concurrent.futures.as_completed(results):
        jobs.append(f.result())
    txyzs = [pdb.replace('pdb', 'txyz_2') for pdb in pdbs]
    with open('txyzlist', 'w') as f:
      for s in txyzs:
        f.write(s + '\n')
  if 'C' in mode:
    if not os.path.isfile('bio.pdb'):
      txyz = pdb.replace("pdb", "txyz")
    else:
      txyz = 'bio.txyz'
      pdb = 'bio.pdb'
    os.system(f"babel {pdb} {txyz}")
    txyzs = np.loadtxt("txyzlist",dtype='str')
    connect(txyz,txyzs)
    if os.path.isfile('solvent.txyz'):
      os.system(f"cat {txyz}_2 solvent.txyz > final.txyz")
    else:
      os.system(f"cp {txyz}_2 final.txyz")