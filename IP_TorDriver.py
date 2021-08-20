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

# 1. Rotate torsional bond and generate QM input
# 2. Fitting to QM obtaining torsional parameters 

# Usage: python ltorsion.py xyz key qmMethod torInterval fitnonezero mode 

helpinfo =\
''' 
   1.0 run ltorsion.py with Mode 1: generate torsional files, including QM and MM\n
   2.0 run ltorsion.py with Mode 2: fitting torsional parameters to the QM energy generated in 1.5\n
'''
 
import sys
import os
import subprocess
import time
import argparse
import numpy as np
from math import e
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

RED = '\33[91m'
GREEN = '\33[92m'
ENDC = '\033[0m'

class Torsion():
  # init
  def __init__(self, torstring, angle):
    if angle >= 0:
      self.fname = "TOR-%s-ANG_+%03d"%(torstring,angle)  
    else:
      self.fname = "TOR-%s-ANG_-%03d"%(torstring,abs(angle))
    self.com = self.fname + ".com"
    self.xyz = self.fname + ".xyz"
    self.key = self.fname + ".key"
    self.log = self.fname + ".log"
    self.out = self.fname + ".out"
    self.ana = self.fname + ".ana"
    self.opt = self.fname + ".opt"
    self.angle = angle
    self.torstr = " ".join(torstring.split("-"))

  # write gaussian input, .com file
  def writeQM(self):
    Atoms = np.loadtxt(inputxyz, usecols=(1,) , dtype="str", unpack=True, skiprows=1)
    Xs, Ys,Zs = np.loadtxt(inputxyz, usecols=(2,3,4) , dtype="float", unpack=True, skiprows=1)
    if not os.path.isfile(self.com):
      with open(self.com, 'w') as f:
        f.write("%chk="+"%s.chk\n"%self.fname)
        f.write("%Nproc=8\n")
        f.write("%Mem=100GB\n")
        f.write("#p %s/%s MaxDisk=100GB IOP(5/13=1) Opt=(ModRedundant, MaxCyc=200)\n\n"%(optmethod, optbasis))
        f.write("Restrained Optimization for torsion %s\n\n"%self.fname)
        f.write("%s %s\n"%(charge, spin))
        for atom, x, y, z in zip(Atoms, Xs, Ys, Zs):
          f.write("%3s%12.6f%12.6f%12.6f\n"%(atom, x,y,z))
        f.write("\n")
        f.write("%s =%6.1f B\n"%(self.torstr, self.angle))
        f.write("%s F\n\n"%(self.torstr))
        f.write("--Link1--\n")
        f.write("%chk="+"%s.chk\n"%self.fname)
        f.write("%Nproc=8\n")
        f.write("%Mem=100GB\n")
        f.write("#p %s/%s Geom=Checkpoint MaxDisk=100GB SP\n\n"%(spmethod,spbasis))
        f.write("Single-point energy for torsion %s\n\n"%self.fname)
        f.write("%s %s\n\n"%(charge, spin))

  # getQM energy from Normal terminated .log file 
  def getQM(self):
    if not os.path.isfile(self.log):
      print(RED + self.log + "does not exist!" + ENDC)
    else:
      lines = open(self.log).readlines()
      if not ("Normal termination" in lines[-1]):
        e = 0.0
      else:
        combinedString = ''.join([line[1:-1] for line in lines[-200:]])
        if spmethod.upper() == "MP2": 
          e = -float(combinedString.split("MP2=-")[1].split("\\")[0])*hartree2kcal
        else:
          e = -float(combinedString.split("HF=-")[1].split("\\")[0])*hartree2kcal
    return e

  # write key for minimize 
  def writeKeymin(self, allTorlistAtom):
    atmnos, atmtyps = np.loadtxt(inputxyz, usecols=(0,5) , dtype="int", unpack=True, skiprows=1)
    AtomTypeDict = {}
    for atmno, atmtyp in zip(atmnos, atmtyps):
      if atmno not in AtomTypeDict:
        AtomTypeDict[atmno] = atmtyp
    subprocess.run("python /home/liuchw/bin/lconvert.py -it g09 -ot txyz -i %s -o %s -tp %s 2>err"%(self.log, self.xyz, inputxyz), shell=True)
    print(GREEN + f"converted {self.xyz}" + ENDC)
    lines = open(inputkey).readlines()
    with open(self.key, "w") as f:
      for line in lines:
        if ("torsion " not in line) or ("#" in line[0]):
          f.write(line)
        else:
          ss = line.split()
          tortypelist = [int(ss[1]), int(ss[2]), int(ss[3]), int(ss[4])]
          flag = 0
          if(all(x in atmtyps for x in tortypelist)): 
            flag = 1
          if flag == 1:
            t1,t2,t3 = float(ss[5]), float(ss[8]), float(ss[11])
            if (t1 != 0):
              ss[5] = "0.000 "
            if (t2 != 0):
              ss[8] = "0.000 "
            if (t3 != 0):
              ss[11] = "0.000 "
            f.write("   ".join(ss) + "\n")
          else:
            f.write(line)
      # restrain the dihedral angle to be fit
      for tor in allTorlistAtom:
        tor0 = Torsion.dihedral(tor[0], tor[1], tor[2], tor[3], self.xyz)
        tor = [str(t) for t in tor]
        f.write("restrain-torsion %s   5.0  %6.1f \n"%("   ".join(tor), float(tor0)))
    return 

  # getMM energy from tinker .ana file 
  def getMM(self):
    for line in open(self.ana).readlines():
      if "Total Potential Energy :" in line:
        e = float(line.split()[-2])
    return e

  # write a template key file  
  @staticmethod
  def writeKeytemplate(tempkey):
    atmnos, atmtyps = np.loadtxt(inputxyz, usecols=(0,5) , dtype="int", unpack=True, skiprows=1)
    atmnames = np.loadtxt(inputxyz, usecols=(1,) , dtype="str", unpack=True, skiprows=1)
    hydrogentypes = []
    for an, at in zip(atmnames, atmtyps):
      if an == "H":
        hydrogentypes.append(at)
    AtomTypeDict = {}
    for atmno, atmtyp in zip(atmnos, atmtyps):
      if atmno not in AtomTypeDict:
        AtomTypeDict[atmno] = atmtyp
    lines = open(inputkey).readlines()
    nprm = 0
    with open(tempkey, "w") as f:
      for line in lines:
        if ("torsion " not in line) or ("#" in line[0]):
          f.write(line)
        else:
          ss = line.split()
          tortypelist = [int(ss[1]), int(ss[2]), int(ss[3]), int(ss[4])]
          flag = 0
          if(all(x in atmtyps for x in tortypelist)): 
            flag = 1
          if flag == 1:
            t1,t2,t3 = float(ss[5]), float(ss[8]), float(ss[11])
            a1,a4 = int(ss[1]), int(ss[4])
            if (a1 in hydrogentypes) and (a4 in hydrogentypes): # H-X-X-H 
              ss[5]  = "0.000 "
              ss[8]  = "0.000 "
              ss[11] = "0.400 "
            else:
              ss[5]  = "PRM_%02d"%(nprm+0)
              ss[8]  = "PRM_%02d"%(nprm+1)
              ss[11] = "PRM_%02d"%(nprm+2)
              nprm += 3
            f.write("   ".join(ss) + "\n")
          else:
            f.write(line)
    return 

  # get torsion list
  @staticmethod
  def getTorlist():
    atmnos, atmtyps = np.loadtxt(inputxyz, usecols=(0,5) , dtype="int", unpack=True, skiprows=1)
    AtomTypeDict = {}
    for atmno, atmtyp in zip(atmnos, atmtyps):
      if atmno not in AtomTypeDict:
        AtomTypeDict[atmno] = atmtyp
    lines = open(inputkey).readlines()
    hasPN = False
    for line in lines:
      if "parameters none" in line:
        hasPN = True
    if not hasPN:
      with open(inputkey, "w") as f:
        f.write("parameters none\n")
        for line in lines:
          f.write(line)
    if not os.path.isfile("final.out"):
      subprocess.run("analyze.x %s -k %s EP > final.out && wait"%(inputxyz,inputkey), shell=True)
    idx = 0
    lines = open("final.out").readlines()
    for line in lines: 
      if "Torsional Angle Parameters :" in line:
        idx = lines.index(line) + 4
    uniqueTorlistType = []
    uniqueTorlistAtom = []
    allTorlistAtom = []
    torPrmListAtom = []
    for i in range(idx, len(lines), 1):
      ss = lines[i].split()
      if len(ss) == 0:
        break
      else:
        torprm = [0.0, 0.0, 0.0]
        if "0/1" in ss:
          torprmidx = ss.index("0/1") - 1
          torprm[0] = float(ss[torprmidx])
        if "180/2" in ss:
          torprmidx = ss.index("180/2") - 1
          torprm[1] = float(ss[torprmidx])
        if "0/3" in ss:
          torprmidx = ss.index("0/3") - 1
          torprm[2] = float(ss[torprmidx])
        torlist0 = [int(ss[1]), int(ss[2]), int(ss[3]), int(ss[4])]
        allTorlistAtom.append(torlist0)
        torlist1 = "-".join(list(map(str,[AtomTypeDict[int(ss[1])], AtomTypeDict[int(ss[2])], AtomTypeDict[int(ss[3])], AtomTypeDict[int(ss[4])]])))
        torlist2 = "-".join(list(map(str,[AtomTypeDict[int(ss[4])], AtomTypeDict[int(ss[3])], AtomTypeDict[int(ss[2])], AtomTypeDict[int(ss[1])]])))
        if (torlist1 not in uniqueTorlistType) and (torlist2 not in uniqueTorlistType):
          uniqueTorlistType.append(torlist1)
          uniqueTorlistAtom.append(torlist0)
          torPrmListAtom.append(torprm)
    return uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom, allTorlistAtom
  
  # get the initial torsion value
  @staticmethod
  def dihedral(n0, n1, n2, n3, xyzfile):
    Xs, Ys,Zs = np.loadtxt(xyzfile, usecols=(2,3,4) , dtype="float", unpack=True, skiprows=1)
    p0 = np.array([Xs[n0-1], Ys[n0-1], Zs[n0-1]]) 
    p1 = np.array([Xs[n1-1], Ys[n1-1], Zs[n1-1]]) 
    p2 = np.array([Xs[n2-1], Ys[n2-1], Zs[n2-1]]) 
    p3 = np.array([Xs[n3-1], Ys[n3-1], Zs[n3-1]]) 
    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x)) 

  # write .sh files
  @staticmethod
  def writeRunMin():
    files = os.listdir()
    with open("runMin.sh", "w") as fo:
      for f in files:
        if ("TOR" in f) and (".xyz" in f):
          fo.write("minimize.x %s -k %s 0.01 > %s &\n"%(f, f.replace(".xyz",".key"), f.replace(".xyz", ".out")))
      fo.write("wait\n")
      for f in files:
        if ("TOR" in f) and (".xyz" in f):
          fo.write("mv %s %s &\n"%(f.replace(".xyz", ".xyz_2"), f.replace(".xyz", ".opt"))) 
      fo.write("wait\n")
    with open("runAna.sh", "w") as fo:
      for f in files:
        if ("TOR" in f) and (".xyz" in f):
          fo.write("analyze.x %s -k tinker.key E > %s &\n"%(f.replace(".xyz", ".opt"), f.replace(".xyz", ".ana")))
      fo.write("wait\n")
    return

# get QM input and runMin.sh
def DataPrep():
  uniqueTorlistAtom, _, _, allTorlistAtom = Torsion.getTorlist()
  subprocess.run("rm -f runMin.sh runAna.sh", shell=True)
  for tor in uniqueTorlistAtom:
    torstring = "-".join(list(map(str,tor)))
    tor0 = int(Torsion.dihedral(tor[0], tor[1], tor[2], tor[3], inputxyz))
    start = int(npoint*interval/2) 
    for ang in range(tor0-start, tor0+start+1, interval):
      if ang <= -180:
        ang += 360
      if ang > 180:
        ang -= 360
      tor = Torsion(torstring, ang)
      tor.writeQM()
      tor.writeKeymin(allTorlistAtom)
  Torsion.writeRunMin()
  print(GREEN + "doing contrain optimization" + ENDC)
  subprocess.run("sh runMin.sh", shell=True)
  Torsion.writeKeytemplate("template.key")
  return 

# cost Function 
def costFuncTor(params):
  subprocess.run("rm -f TOR*.ana && wait",shell=True)
  uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom, _ = Torsion.getTorlist()
  files = os.listdir(os.getcwd())
  with open("tinker.key", "w") as key2:
    for line in open("template.key").readlines():
      if ("torsion " not in line):
        key2.write(line)
      else:
        sss = line.split()
        if "PRM_" in sss[5]:
          idx = int(sss[5].split("_")[1])
          prmstr = "PRM_%02d"%idx
          line = line.replace(prmstr, str("%15.8f"%params[idx]))
        if "PRM_" in sss[8]:
          idx = int(sss[8].split("_")[1])
          prmstr = "PRM_%02d"%idx
          line = line.replace(prmstr, str("%15.8f"%params[idx]))
        if "PRM_" in sss[11]:
          idx = int(sss[11].split("_")[1])
          prmstr = "PRM_%02d"%idx
          line = line.replace(prmstr, str("%15.8f"%params[idx]))
        key2.write(line)
  subprocess.run("sh runAna.sh",shell=True)
  grepstr = "grep 'Total Potential Energy :' TOR*.ana > grep.dat" 
  nline1 = len(open("runAna.sh").readlines()) - 1
  while True:
    subprocess.run(grepstr, shell=True)
    nline2 = len(open("grep.dat").readlines())
    if nline1 == nline2:
      break
    else:
      time.sleep(1.0)
  QM = []; MM = []
  for tor,torT,prm in zip(uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom):
    torstring = "-".join(list(map(str,tor)))
    tor0 = int(Torsion.dihedral(tor[0], tor[1], tor[2], tor[3], inputxyz))
    start = int(npoint*interval/2) 
    for ang in range(tor0-start, tor0+start+1, interval):
      if ang <=-180:
        ang += 360
      if ang > 180:
        ang -= 360
      tor = Torsion(torstring, ang)
      qm = tor.getQM()
      if (qm != 0.0):
        QM.append(qm)
        MM.append(tor.getMM())
  QM = np.array(QM) - np.array(QM).min() 
  MM = np.array(MM) - np.array(MM).min() 
  for i in range(len(QM)):
    if (QM[i] > 20.0) or (MM[i]>20.0):
      QM[i] = 0.0 
      MM[i] = 0.0
  print(GREEN + " current RMSE is %10.5f"%np.sqrt(np.square(QM-MM).mean()) + ENDC)
  weight = []
  for q in QM:
    if 5.0 <= q <= 10.0:
      weight.append(1.0)
    elif q < 5.0:
      # [2.718, 1]
      weight.append(e*np.exp(-q/5.0))
    else:
      # [1.0, 0.367]
      weight.append(np.exp((1.0-q/10.0)))
  return (QM-MM)*weight

# fitting all torsion parameters together
def fittingTorsion():
  lines = open("template.key").readlines()
  for line in lines:
    if "PRM_" in line:
      nparams = int(line.split("PRM_")[-1].split()[0]) + 1
  x0 = np.zeros(nparams)
  lowerBound = []
  upperBound = []
  for x in x0:
    lowerBound.append(-8.0)
    upperBound.append(+8.0)
  if restrain == 1:
    print(RED + "Doing torsion fitting with restraint" + ENDC)
    ret = least_squares(costFuncTor, x0, bounds = (lowerBound, upperBound), verbose=2, diff_step=0.0001, ftol=0.0001, gtol=0.0001, xtol=0.0001)
  else:
    print(RED + "Doing torsion fitting without restraint" + ENDC)
    ret = least_squares(costFuncTor, x0, verbose=2, diff_step=0.0001, ftol=0.0001, gtol=0.0001, xtol=0.0001)

# plot torsion energy profile and write tinker.key
def plotData():
  uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom, _ = Torsion.getTorlist()
  for tor,torT,prm in zip(uniqueTorlistAtom, uniqueTorlistType, torPrmListAtom):
    QM = []; MM = []; Ang = []
    torstring = "-".join(list(map(str,tor)))
    tor0 = int(Torsion.dihedral(tor[0], tor[1], tor[2], tor[3], inputxyz))
    start = int(npoint*interval/2) 
    for ang in range(tor0-start, tor0+start+1, interval):
      if ang <= -180:
        ang += 360
      if ang > 180:
        ang -= 360
      tor = Torsion(torstring, ang)
      qm = tor.getQM()
      if (qm != 0.0):
        QM.append(qm)
        MM.append(tor.getMM())
        Ang.append(tor.angle)
    QM = np.array(QM) - np.array(QM).min()
    MM = np.array(MM) - np.array(MM).min()
    for i in range(len(QM)):
      if (QM[i] > 20.0) or (MM[i] >20.0):
        QM[i] = 0.0
        MM[i] = 0.0
    fig, ax = plt.subplots()

    QM = [x for _,x in sorted(zip(Ang,QM))]
    MM = [x for _,x in sorted(zip(Ang,MM))]
    Ang = sorted(Ang)
    ax.plot(Ang, QM, marker='o', markersize=10, fillstyle='none', label="QM")
    ax.plot(Ang, MM, marker='x', markersize=10, label="MM")
    ax.set(xlabel='Torsional Angle (deg.)', ylabel='Energy (kcal/mol)', title='Torsion-%s'%torstring)
    ax.grid()
    ax.legend(loc="best")
    fig.savefig("torsion_fitting-%s.png"%torstring)
    plt.close() 
  return

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-xyz',  dest = 'inputxyz', required=True,          type=str, help="tinker xyz file")  
  parser.add_argument('-key',  dest = 'inputkey', required=True,          type=str, help="tinker key/prm file")  
  parser.add_argument('-optm', dest = 'optmethod',default="MP2",          type=str, help="OPT method for QM torsion scan")
  parser.add_argument('-optb', dest = 'optbasis', default="6-311++G(d,p)",type=str, help="OPT basis for QM torsion scan")
  parser.add_argument('-spm',  dest = 'spmethod', default="MP2",          type=str, help="SP method for QM torsion scan")
  parser.add_argument('-spb',  dest = 'spbasis',  default="aug-cc-pvtz",  type=str, help="SP basis for QM torsion scan")
  parser.add_argument('-chrg', dest = 'charge',   default="0",            type=str, help="Total charge of the molecule")
  parser.add_argument('-spin', dest = 'spin',     default="1",            type=str, help="Total spin of the molecule")
  parser.add_argument('-mode', dest = 'mode',     required=True,          type=int, help="Mode 1: prepare Gaussian .com, tinker key and submit files;\
                                                                                          Mode 2: tune torsion parameters to fit QM and plot results")
  parser.add_argument('-npt',  dest = 'npoint',   default=12,             type=int, help="Number of datapoints per torsional angle. Default: 12")
  parser.add_argument('-intv', dest = 'interval', default=30,             type=int, help="Torsional angle interval in degree. Default: 30")
  parser.add_argument('-rstr', dest = 'restrain', default=0,              type=int, help="Use restrain in torsion fitting? 0 for NO; 1 for YES")
  args = vars(parser.parse_args())

  global inputxyz,inputkey
  global optmethod,spmethod
  global optbasis,spbasis
  global charge,spin
  global npoint,interval,restrain
  global hartree2kcal,degree2radian
  inputxyz = args["inputxyz"] 
  inputkey = args["inputkey"] 
  optmethod = args["optmethod"] 
  spmethod = args["spmethod"] 
  optbasis = args["optbasis"] 
  spbasis = args["spbasis"] 
  charge = args["charge"] 
  spin = args["spin"]
  mode = args["mode"] 
  npoint = args["npoint"] 
  interval = args["interval"] 
  restrain = args["restrain"]
  hartree2kcal = 627.5094740631 
  degree2radian = np.pi/180.0
  # Execute different MODE
  if mode == 1:
    DataPrep()
  # After mode 1, QM jobs should be run
  # python ~/bin/lsub.py -i TOR*.com
  elif mode == 2:
    fittingTorsion()
    plotData()
  else:
    print(GREEN + helpinfo + ENDC)

if __name__ == '__main__':
  main()
