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
from utilities.extrabasis import *

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
    xyzref = "input.xyz" 
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
  if mode in "DMA ESP".split():
    if not os.path.isfile(optxyzfile):
      subprocess.run("babel -ig09 %s -oxyz %s"%(optlogfile, optxyzfile), shell=True)
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
    f.write("Title %s gdmain\n\n"%fname)
    f.write("File %s density MP2\n"%dmafchkfile)
    f.write("Angstrom\nAU\nMultipoles\nSwitch 0\nLimit 2\n")
    f.write("Punch %s\n"%dmapunchfile)
    f.write("Radius H 0.60\nRadius S 0.85\nRadius P 0.95\nRadius Cl 1.10\n\n")
    f.write("Start\n\nFinish\n")
  gdmacmd = "gdma < %s > %s"%(gdmain, gdmaout)
  subprocess.run(gdmacmd,shell=True)

  with open(poleditin, "w") as f:
    f.write("\nA\n\n2\nY\n\nY\n")
  if os.path.isfile(fname + ".key"):
    subprocess.run("rm -f %s.key"%fname, shell=True)
  poledit1 = "poledit.x 1 %s %s < %s" %(gdmaout, prmheader, poleditin)
  subprocess.run(poledit1, shell=True)

  subprocess.run("echo 'parameters %s' | cat - %s.key > %s.key_2"%(prmheader, fname, fname), shell=True) 
  valencecmd = "valence.x 1 %s.xyz %s %s.key_2 > %s"%(fname, optlogfile, fname, valenceout)
  subprocess.run(valencecmd, shell=True)

  if os.path.isfile(gridfile):
    subprocess.run("rm -rf %s"%gridfile, shell=True)
  if os.path.isfile(potfile):
    subprocess.run("rm -rf %s"%potfile, shell=True)
  potential1 = "potential.x 1 %s.xyz -k %s.key_2"%(fname, fname)
  subprocess.run(potential1, shell=True)
  cubegen0 = "cubegen 0 potential=MP2 %s %s -5 h < %s"%(espfchkfile, cubefile, gridfile) 
  subprocess.run(cubegen0, shell=True)
  potential2 = "potential.x 2 %s"%cubefile
  subprocess.run(potential2, shell=True)
  with open(fname + ".key_3", "w") as f:
    for line in open(fname+".key_2").readlines():
      if "PARAMETERS" in line.upper():
        f.write(line)
        f.write("potential-offset 1.0\nfix-monopole\n")
      else:
        f.write(line)
    for line in open(valenceout).readlines():
      if ("bond " in line) or ("angle " in line) or ("anglep " in line) or ("opbend " in line) or ("strbnd" in line) or ("vdw " in line) or ("torsion " in line):
        f.write(line[1:])
  if os.path.isfile(fname + ".key_4"):
    subprocess.run("rm -f %s.key_4"%fname, shell=True)
  potential6 = "potential.x 6 %s.xyz -k %s.key_3 %s N 0.1 > %s.potfitlog" %(fname, fname, potfile, fname)
  subprocess.run(potential6, shell=True)
  
  lines = open(fname+".key_4").readlines()
  natom = 0
  for line in lines:
    if "atom " in line:
      natom += 1
    if "multipole " in line:
      idxmpole = lines.index(line)
      break
  for n in range(idxmpole, natom*5+idxmpole, 1):
    lines[n] = "\n" 
  with open(fname+".key_5", "w") as f:
    for line in lines:
      if line != "\n":
        f.write(line)
  atomtyper = os.path.join(rootdir, "IP_AtomTyper.py")
  subprocess.run("python %s -xyz %s.xyz -prm %s.key_5 -idx %s"%(atomtyper, fname, fname, atmidx), shell=True)
  os.rename(f"{fname}.key_5_1", "final.prm")
  os.rename(f"{fname}.xyz_1", "final.xyz")
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-input', dest = 'input', required=True, help='input file, can be xyz')  
  parser.add_argument('-charge', dest = 'charge', default = "0", help='total charge')
  parser.add_argument('-spin',  dest = 'spin', default = "1", help='total spin')
  parser.add_argument('-index', dest = 'index', default = "401", help='amoeba atom-type starting index')
  parser.add_argument('-qmonly', dest = 'qmonly', default = 0, help='do QM only; 0 or 1', type=int)
  parser.add_argument('-disk',  dest = 'disk', default = "200", help='QM disk [unit GB]')
  parser.add_argument('-memory', dest = 'memory', default = "100", help='QM memory [unit GB]')
  args = vars(parser.parse_args())
  
  infile = args["input"]
  if os.path.isfile(f"{infile}"):
    subprocess.run(f"cp {infile} {infile}.b",shell=True)
    subprocess.run(f"mv {infile} input.xyz",shell=True)
  else:
    if os.path.isfile(f"{infile}.b"):
      print(YELLOW + f"Warning: using {infile}.b as input file" + ENDC)
      subprocess.run(f"cp {infile}.b input.xyz",shell=True)
    else:
      sys.exit(RED + f"Error: could not find {infile} or {infile}.b"+ ENDC)

  global fname, inputxyz
  fname, ext = os.path.splitext(infile)
  inputxyz = "input.xyz"
 
  global memory, disk
  memory = args["memory"]
  disk = args["disk"]

  global chgspin, atmidx, comment, nproc 
  xyz = ".xyz"
  charge = args["charge"]
  spin = args["spin"]
  atmidx = args["index"]
  qmonly = args["qmonly"]
  chgspin = " ".join([charge, spin, "\n"])
  comment = " Generated by TIPTOP/IP_ParmGen.py\n"
  nproc = "%Nproc=10\n"

  global cubefile, gridfile, potfile, xyzfile
  cubefile = f"{fname}.cube"
  gridfile = f"{fname}.grid"
  potfile =  f"{fname}.pot"
  xyzfile =  f"{fname}.xyz"

  global optcomfile, optchkfile, optlogfile, optxyzfile
  optcomfile =  f"{fname}-opt.com"
  optchkfile =  f"{fname}-opt.chk"
  optlogfile =  f"{fname}-opt.log"
  optxyzfile =  f"{fname}-opt.xyz"

  global dmacomfile, dmachkfile, dmafchkfile, dmapunchfile, dmalogfile
  dmacomfile  =  f"{fname}-dma.com"
  dmachkfile  =  f"{fname}-dma.chk"
  dmalogfile  =  f"{fname}-dma.log"
  dmafchkfile  = f"{fname}-dma.fchk"
  dmapunchfile = f"{fname}-dma.punch"

  global espcomfile, espchkfile, espfchkfile, esplogfile
  espcomfile =  f"{fname}-esp.com"
  espchkfile =  f"{fname}-esp.chk"
  espfchkfile = f"{fname}-esp.fchk"
  esplogfile =  f"{fname}-esp.log"

  global gdmain, gdmaout, poleditin, valenceout
  gdmain     =  f"{fname}.gdmain"
  gdmaout    =  f"{fname}.gdmaout"
  poleditin  =  f"{fname}.poleditin"
  valenceout =  f"{fname}.valenceout"

  global rootdir, prmheader
  rootdir = os.path.join(os.path.split(__file__)[0])
  prmheader  = os.path.join(rootdir, "database.ParmGen", "amoeba_header.prm")

  if (not os.path.isfile(optlogfile)):
    f_GAUSSIAN("OPT")
  opt_finished = "Normal termination" in open(optlogfile).readlines()[-1]
  if (opt_finished and not os.path.isfile(dmalogfile)):
    f_GAUSSIAN("DMA")
  if (opt_finished and not os.path.isfile(esplogfile)):
    f_GAUSSIAN("ESP")

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
