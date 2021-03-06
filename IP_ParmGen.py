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
import subprocess
import numpy as np
from pybel import *
import concurrent.futures 
from utilities.extrabasis import *
from utilities.polarize import *

# color
RED = '\033[91m'
ENDC = '\033[0m'
GREEN = '\033[92m'
YELLOW = '\033[93m'

def f_GAUSSIAN(mode):
  print(GREEN + f"\n  Doing {mode} ...\n" + ENDC)
  qm = "MP2"
  if mode == "OPT":
    jt = "OPT(MaxCycle=200)"
    bf = "6-31G*"
    bf_iodine = "6-311gd.gbs"
    bf_other = "6-311gd.gbs"
    comfile = optcomfile
    logfile = optlogfile
    chkfile = optchkfile
    extra = " IOP(5/13=1) "
    pseudo = ' '
    xyzref = "initxyz.xyz" 
  if mode == "DMA":
    jt = "SP Density=MP2 "
    bf = "6-311G(d,p)"
    bf_iodine = "def2-svp.gbs"
    bf_other = "6-311gdp.gbs"
    comfile = dmacomfile
    logfile = dmalogfile
    chkfile = dmachkfile
    pseudo = "Pseudo=Read "
    extra = ' '
    xyzref = optxyzfile 
  if mode == "ESP":
    jt = "SP Density=MP2 "
    bf = "aug-cc-pvtz"
    bf_iodine = "def2-tzvppd.gbs"
    bf_other = bf_iodine 
    comfile = espcomfile
    logfile = esplogfile
    chkfile = espchkfile
    pseudo = "Pseudo=Read "
    extra = ' '
    xyzref = optxyzfile
  atoms = np.loadtxt(xyzref, usecols=(0,), dtype='str', skiprows=2, unpack=True) 
  xs,ys,zs = np.loadtxt(xyzref, usecols=(1,2,3), dtype='float', skiprows=2, unpack=True) 
  mem  = f"%Mem={memory}GB \n"
  dk = f"{disk}GB"
  if "I" not in atoms:
    key = ' '.join(["#P", qm + "/"+ bf, extra, jt, f"MaxDisk={dk}\n\n"])
  else:
    key = ' '.join(["#P", qm + "/Gen ", pseudo, extra, jt, f"MaxDisk={dk}\n\n"])
  with open(comfile, 'w') as fout:
    fout.write("%chk=" + chkfile+"\n" + mem + nproc + key + comment + "\n" + chgspin)
    for atom, x,y,z in zip(atoms, xs, ys, zs):
      formatted = "%5s%12.6f%12.6f%12.6f\n"%(atom, x, y, z)
      fout.write(formatted)
    fout.write("\n")
    # special case for Iodine
    if "I" in atoms:
      uniqueatoms = list(set(atoms))
      for ua in uniqueatoms:
        if ua == "I":
          basislines = getBasis(ua, os.path.join(rootdir, "database.ParmGen", bf_iodine))
          for bl in basislines:
            fout.write(bl)
        else:
          basislines = getBasis(ua, os.path.join(rootdir, "database.ParmGen", bf_other))
          for bl in basislines:
            fout.write(bl)
      fout.write("\n")
      if mode in 'DMA ESP'.split(): 
        for ua in uniqueatoms:
          ecplines = getECP(ua, os.path.join(rootdir, "database.ParmGen", bf_iodine))
          for el in ecplines:
            fout.write(el)
        fout.write("\n")
  subprocess.run(f"g09 {comfile} {logfile}; wait", shell=True)
  if mode in 'DMA ESP'.split():
    subprocess.run(f"formchk {chkfile}", shell=True)
  return

def f_TINKER():
  print(GREEN + "\n  Doing TINKER ...\n" + ENDC)
  if not os.path.isfile(dmafchkfile):
    subprocess.run(f"formchk {dmachkfile}",shell=True)
  if not os.path.isfile(espfchkfile):
    subprocess.run(f"formchk {espchkfile}",shell=True)
  
  with open(gdmain, "w") as f:
    f.write("Title %s gdmain\n\n"%prefix)
    f.write("File %s density MP2\n"%dmafchkfile)
    f.write("Angstrom\nAU\nMultipoles\nSwitch 0\nLimit 2\n")
    f.write("Punch %s\n"%dmapunchfile)
    f.write("Radius H 0.65\nRadius S 0.85\nRadius P 0.95\nRadius Cl 1.10\n\n")
    f.write("Start\n\nFinish\n")
  gdmacmd = "gdma < %s > %s"%(gdmain, gdmaout)
  subprocess.run(gdmacmd,shell=True)

  polars = getpolar(optxyzfile, os.path.join(rootdir, "database.ParmGen"), polarset)
  with open(poleditin, "w") as f:
    f.write("\nA\n")
    #f.write("\nP\n")
    for p in polars:
      f.write(f"{p} {polars[p]}\n")
    f.write("\n2\nY\n\nY\n")
  if os.path.isfile(prefix + ".key"):
    subprocess.run("rm -f %s.key"%prefix, shell=True)
  poledit1 = "poledit.x 1 %s %s < %s" %(gdmaout, prmheader, poleditin)
  subprocess.run(poledit1, shell=True)

  subprocess.run("echo 'parameters %s' | cat - %s.key > %s.key_2"%(prmheader, prefix, prefix), shell=True) 
  valencecmd = "valence.x 1 %s.xyz %s %s.key_2 > %s"%(prefix, optlogfile, prefix, valenceout)
  subprocess.run(valencecmd, shell=True)

  if os.path.isfile(gridfile):
    subprocess.run("rm -rf %s"%gridfile, shell=True)
  if os.path.isfile(potfile):
    subprocess.run("rm -rf %s"%potfile, shell=True)
  potential1 = "potential.x 1 %s.xyz -k %s.key_2"%(prefix, prefix)
  subprocess.run(potential1, shell=True)
  cubegen0 = "cubegen 0 potential=MP2 %s %s -5 h < %s"%(espfchkfile, cubefile, gridfile) 
  subprocess.run(cubegen0, shell=True)
  potential2 = "potential.x 2 %s"%cubefile
  subprocess.run(potential2, shell=True)
  with open(prefix + ".key_3", "w") as f:
    for line in open(prefix+".key_2").readlines():
      if "PARAMETERS" in line.upper():
        f.write(line)
        f.write("potential-offset 1.0\nfix-monopole\n")
      else:
        f.write(line)
    for line in open(valenceout).readlines():
      if ("bond " in line) or ("angle " in line) or ("anglep " in line) or ("opbend " in line) or ("strbnd" in line) or ("vdw " in line) or ("torsion " in line):
        f.write(line[1:])
  if os.path.isfile(prefix + ".key_4"):
    subprocess.run("rm -f %s.key_4"%prefix, shell=True)
  potential6 = "potential.x 6 %s.xyz -k %s.key_3 %s N 0.1 > %s.potfitlog" %(prefix, prefix, potfile, prefix)
  subprocess.run(potential6, shell=True)
  
  lines = open(prefix+".key_4").readlines()
  natom = 0
  for line in lines:
    if "atom " in line:
      natom += 1
    if "multipole " in line:
      idxmpole = lines.index(line)
      break
  for n in range(idxmpole, natom*5+idxmpole, 1):
    lines[n] = "\n" 
  with open(prefix+".key_5", "w") as f:
    for line in lines:
      if line != "\n":
        f.write(line)
  atomtyper = os.path.join(rootdir, "IP_AtomTyper.py")
  subprocess.run("python %s -xyz %s.xyz -prm %s.key_5 -idx %s"%(atomtyper, prefix, prefix, atmidx), shell=True)
  os.rename(f"{prefix}.key_5_1", "final.prm")
  os.rename(f"{prefix}.xyz_1", "final.xyz")
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-molecule', dest = 'molecule', required=True, help='input file in mol2 format')  
  parser.add_argument('-charge', dest = 'charge', default = "0", help='total charge')
  parser.add_argument('-spin',  dest = 'spin', default = "1", help='total spin')
  parser.add_argument('-index', dest = 'index', default = "401", help='amoeba atom-type starting index')
  parser.add_argument('-qmonly', dest = 'qmonly', default = 0, help='do QM only; 0 or 1', type=int)
  parser.add_argument('-polarset', dest = 'polarset', default = "amoeba", help='element based `amoeba` polarizability or updated `amoeba+` set')
  parser.add_argument('-optonly', dest = 'optonly', default = 0, help='do QM opt only; 0 or 1', type=int)
  parser.add_argument('-disk',  dest = 'disk', default = "200", help='QM disk [unit GB]')
  parser.add_argument('-memory', dest = 'memory', default = "100", help='QM memory [unit GB]')
  args = vars(parser.parse_args())
  
  infile = args["molecule"]
  if os.path.isfile(f"{infile}"):
    subprocess.run(f"babel -imol2 {infile} -oxyz initxyz.xyz", shell=True) 
  else:
    sys.exit(RED + f"Error: could not find {infile} "+ ENDC)

  global prefix
  prefix, ext = os.path.splitext(infile)
  
  global polarset
  polarset = args["polarset"]

  global memory, disk
  memory = args["memory"]
  disk = args["disk"]

  global chgspin, atmidx, comment, nproc 
  xyz = ".xyz"
  charge = args["charge"]
  spin = args["spin"]
  atmidx = args["index"]
  optonly = args["optonly"]
  qmonly = args["qmonly"]
  if optonly == 1:
    qmonly = 1
  chgspin = " ".join([charge, spin, "\n"])
  comment = " Generated by TIPTOP/IP_ParmGen.py\n"
  nproc = "%Nproc=10\n"

  global cubefile, gridfile, potfile, xyzfile
  cubefile = f"{prefix}.cube"
  gridfile = f"{prefix}.grid"
  potfile =  f"{prefix}.pot"
  xyzfile =  f"{prefix}.xyz"

  global optcomfile, optchkfile, optlogfile, optxyzfile
  optcomfile =  f"{prefix}-opt.com"
  optchkfile =  f"{prefix}-opt.chk"
  optlogfile =  f"{prefix}-opt.log"
  optxyzfile =  f"{prefix}-opt.xyz"

  global dmacomfile, dmachkfile, dmafchkfile, dmapunchfile, dmalogfile
  dmacomfile  =  f"{prefix}-dma.com"
  dmachkfile  =  f"{prefix}-dma.chk"
  dmalogfile  =  f"{prefix}-dma.log"
  dmafchkfile  = f"{prefix}-dma.fchk"
  dmapunchfile = f"{prefix}-dma.punch"

  global espcomfile, espchkfile, espfchkfile, esplogfile
  espcomfile =  f"{prefix}-esp.com"
  espchkfile =  f"{prefix}-esp.chk"
  espfchkfile = f"{prefix}-esp.fchk"
  esplogfile =  f"{prefix}-esp.log"

  global gdmain, gdmaout, poleditin, valenceout
  gdmain     =  f"{prefix}.gdmain"
  gdmaout    =  f"{prefix}.gdmaout"
  poleditin  =  f"{prefix}.poleditin"
  valenceout =  f"{prefix}.valenceout"

  global rootdir, prmheader
  rootdir = os.path.join(os.path.split(__file__)[0])
  prmheader  = os.path.join(rootdir, "database.ParmGen", "amoeba_header.prm")
  #prmheader  = os.path.join(rootdir, "database.ParmGen", "amoebaplus21_header.prm")

  if (not os.path.isfile(optlogfile)):
    f_GAUSSIAN("OPT")
  
  modes = [] 
  opt_finished = "Normal termination" in open(optlogfile).readlines()[-1]
  if (optonly == 0):
    if (opt_finished) and not os.path.isfile(optxyzfile): 
      subprocess.run("babel -ig09 %s -oxyz %s"%(optlogfile, optxyzfile), shell=True)
    if (opt_finished and not os.path.isfile(dmalogfile)):
      modes.append("DMA") 
    if (opt_finished and not os.path.isfile(esplogfile)):
      modes.append("ESP") 
  jobs = []
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(f_GAUSSIAN, mode) for mode in modes]
    for f in concurrent.futures.as_completed(results):
      jobs.append(f.result())
  if (qmonly == 0):
    dma_finished = "Normal termination" in open(dmalogfile).readlines()[-1]
    esp_finished = "Normal termination" in open(esplogfile).readlines()[-1]
    if dma_finished and esp_finished:
      f_TINKER()
    else:
      sys.exit(RED + "Skip f_TINKER(), dma/esp not finished!" + ENDC)

  print(GREEN + "\n  ParmGen jobs finished successfully!\n" + ENDC)

if __name__ == "__main__":
  main()
