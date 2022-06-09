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


''' Usage: python IP_MatchTXYZ.py -t template.txyz -d dealwith.(t)xyz 
    ~ Assign the atom types of the template.txyz file to the dealwith txyz/xyz/pdb file.
    ~ Especially suitable for the case that their atoms are in different order.
    ~ It uses connectivity information only to match the structures
    ~ make sure you can run obabel command before you use this program
'''

import os
import sys
import argparse
import numpy as np

def readTXYZ(TXYZ):
  atoms=[]
  coord=[]
  order=[]
  types=[]
  connections=[]
  for line in open(TXYZ).readlines()[1:]: 
    data=line.split()
    order.append(data[0])
    atoms.append(data[1])
    coord.append([float(data[2]), float(data[3]), float(data[4])])
    types.append(data[5])
    connections.append(data[6:])
  return atoms,coord,order,types,connections

def fingerprint(TXYZ):
  fprints = []
  atoms, elements = np.loadtxt(TXYZ, usecols=(0,1), dtype='str', skiprows=1, unpack=True)
  connections = []
  for line in open(TXYZ).readlines()[1:]:
    d = line.split()
    connections.append(d[6:])
  
  atom_ele_dict = dict(zip(atoms, elements))
  atom_con_dict = {}
  for atom, con in zip(atoms,connections):
    con_ele = [atom_ele_dict[c] for c in con] 
    constr = ''.join(sorted(con_ele)) 
    atom_con_dict[atom] = constr

  level = 5 
  if level > 1:
    atom_con_dict2 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict[c] for c in cons])
      atom_con_dict2[atom] = ''.join(sorted(newstr))

  # level 3 is good for chain molecules 
  if level > 2:
    atom_con_dict3 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict2[c] for c in cons])
      atom_con_dict3[atom] = ''.join(sorted(newstr))

  # level 4 is needed for ring molecules 
  if level > 3:
    atom_con_dict4 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict3[c] for c in cons])
      atom_con_dict4[atom] = ''.join(sorted(newstr))
  
  if level > 4:
    atom_con_dict5 = {}
    for atom, con in zip(atoms,connections):
      eles = []
      cons = []
      for c in con:
        eles.append(atom_ele_dict[c])
        cons.append(c)
      cons = [x for _,x in sorted(zip(eles,cons))]
      newstr = ''.join([atom_con_dict4[c] for c in cons])
      atom_con_dict5[atom] = ''.join(sorted(newstr))
  
  for atom in atoms:
    fprints.append(atom_ele_dict[atom] + '-' + str(''.join(sorted(atom_con_dict[atom] + \
    atom_con_dict2[atom] + atom_con_dict3[atom] + atom_con_dict4[atom] + atom_con_dict5[atom]))))
  
  return fprints

def comparefingerprints(fp1, fp2):
  match = True
  newidx = []
  if sortatom:
    for i in fp1:
      if i in fp2:
        idx = fp2.index(i)
        newidx.append(idx)
        fp2[idx] = ' '
      else:
        match = False
        break
  else:
    for i in fp2:
      if i in fp1:
        idx = fp1.index(i)
        newidx.append(idx)
        fp1[idx] = ' '
      else:
        match = False
        break
  return match, newidx

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', dest = 'template', nargs='+', help = "Template txyz file(s)", required=True)  
  parser.add_argument('-d', dest = 'dealwith', help = "File to deal-with, can be .xyz/.txyz/.pdb", required=True)  
  parser.add_argument('-s', dest = 'sortatom', help = "sort atoms according to template file", default=False, type=bool)  
  args = vars(parser.parse_args())
  templates = args["template"]
  dealwith  = args["dealwith"]
  global sortatom
  sortatom  = args["sortatom"]
  
  if os.path.splitext(dealwith)[1] == ".xyz":
    lines = open(dealwith).readlines()
    if not (int(lines[0].split()[0]) == len(lines)-1):
      xyz = dealwith
      dealwith = dealwith.replace("xyz", "txyz")
      obstr = "obabel -ixyz %s -otxyz -O %s"%(xyz, dealwith)
      os.system(obstr)
  if os.path.splitext(dealwith)[1] == ".pdb":
    xyz = dealwith
    dealwith = dealwith.replace("pdb", "txyz")
    obstr = "obabel -ipdb  %s -otxyz -O %s"%(xyz, dealwith)
    os.system(obstr)
  
  match = False
  newidx = []
  fname = dealwith + "_2"
  for template in templates:
    atoms1, coord1, _, types1, connections1 =  readTXYZ(dealwith)
    atoms2, coord2, _, types2, connections2 =  readTXYZ(template)
    if len(atoms1) == len(atoms2):
      fp1 = fingerprint(template)
      fp2 = fingerprint(dealwith)
      match, newidx = comparefingerprints(fp1, fp2)
    if not match: 
      print(f"Could not match {template} and {dealwith}")
    else:
      print(f"Matched {template}. New file generated {dealwith}_2")
      with open(fname, 'w') as f:
        f.write("%3s\n"%len(atoms1))
        for i in range(len(newidx)):
          idx = int(newidx[i])
          if sortatom:
            f.write("%3s%3s%12.6f%12.6f%12.6f  %s   %s\n"%(i+1,atoms1[idx], coord1[idx][0], coord1[idx][1], coord1[idx][2], types2[i], '  '.join(connections2[i])))
          else:
            f.write("%3s%3s%12.6f%12.6f%12.6f  %s   %s\n"%(i+1,atoms1[i], coord1[i][0], coord1[i][1], coord1[i][2], types2[idx], '  '.join(connections1[i])))
      with open("mapping.lst", 'w') as f:
        for i in range(len(newidx)):
          f.write("%6d%6d\n"%(i, newidx[i]))
      print("The atom order mapping has been saved in mapping.lst")
