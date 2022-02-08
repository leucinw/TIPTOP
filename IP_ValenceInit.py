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
from pybel import *
import networkx as nx

# Print the zero parameters for bond/angle/strbnd/opbend/torsion
# Given the tinker xyz file
# Assumption: atom type == atom class

def calc_bond_angle(coords):
  
  # bond
  # coords = [[x1,y1,z1], [x2,y2,z2]]
  result = 0
  ndata = len(coords)
  if ndata == 2:
    coord1 = np.array(coords[0])
    coord2 = np.array(coords[1])
    result = np.sqrt(np.square(coord1-coord2).sum())
  
  # angle
  # coords = [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
  if ndata == 3:
    coord1 = np.array(coords[0])
    coord2 = np.array(coords[1])
    coord3 = np.array(coords[2])
    vec21 = coord1 - coord2 
    vec23 = coord3 - coord2
    dot = np.dot(vec21, vec23)
    vec21norm = np.linalg.norm(vec21) 
    vec23norm = np.linalg.norm(vec23) 
    result = 180.0/np.pi * (np.arccos(dot/(vec21norm*vec23norm)))
  return result

def print_initial_parameters(txyz):
  atom2type = {}
  g = nx.Graph()
  nodes = []
  edges = []
  lines = open(txyz).readlines()
  coords = []
  if len(lines[0].split()) == 1:
    os.system(f"sed '1 s/$/ xxx/g' -i {txyz}")
  for line in lines[1:]:
    d = line.split()
    coords.append([float(d[2]), float(d[3]), float(d[4])])
    if int(d[0]) not in atom2type:
      atom2type[int(d[0])] = int(d[5])
    if d[0] not in nodes: 
      nodes.append(int(d[0]))
    for c in d[6:]:
      s = sorted([int(d[0]), int(c)])
      if s not in edges:
        edges.append(s)
  g.add_nodes_from(nodes)
  g.add_edges_from(edges)
 
  # find bonds and angles
  bonds = {} 
  angles = {} 
  for edge in g.edges:
    (n1, n2) = edge
    t1, t2 = atom2type[n1], atom2type[n2]
    comb = '-'.join([str(s) for s in sorted([t1,t2])])
    b = calc_bond_angle([coords[n1-1], coords[n2-1]])
    if comb not in bonds:
      bonds[comb] = [b]
    else:
      bonds[comb] += [b]
  
  for node in g.nodes:
    adjs = list(g.adj[node])
    if len(adjs) > 1:
      for i in range(len(adjs)-1):
        for j in range(i+1, len(adjs)):
          n1, n2, n3 = adjs[i], node, adjs[j]
          t1, t2, t3 = atom2type[n1], atom2type[n2], atom2type[n3]
          if t1 > t3:
            t1, t3 = t3, t1
            n1, n3 = n3, n1
          a = calc_bond_angle([coords[n1-1], coords[n2-1], coords[n3-1]])
          comb = '-'.join([str(s) for s in [t1,t2,t3]])
          if comb not in angles:
            angles[comb] = [a]
          else:
            angles[comb] += [a]
  
  # find torsions
  torsions = []
  for edge in g.edges:
    [a1, a2] = edge
    if (g.degree[a1] >=2) and (g.degree[a2] >= 2):
      a1_adjs = list(g.adj[a1])
      a2_adjs = list(g.adj[a2])
      for x in a1_adjs:
        for y in a2_adjs:
          if (x != a2) and (y != a1):
            if ([atom2type[x], atom2type[a1], atom2type[a2], atom2type[y]] not in torsions) \
            and ([atom2type[y], atom2type[a2], atom2type[a1], atom2type[x]] not in torsions) :
              torsions.append([atom2type[x], atom2type[a1], atom2type[a2], atom2type[y]])
  
  # find the tri- center
  tricentertypes = []
  for mol in readfile("txyz", txyz):
    natoms = len(mol.atoms)
    for i in range(natoms):
      hyb = mol.atoms[i].hyb
      if (hyb == 2) and (g.degree[i+1] == 3):
        if atom2type[i+1] not in tricentertypes:
          tricentertypes.append(str(atom2type[i+1]))
  # print  
  for bond in bonds:
    b = bond.split('-') 
    comb = '-'.join([str(s) for s in b])
    val = np.array(bonds[comb]).mean()
    print(f"bond   {b[0]}  {b[1]}  0.00 {val:10.4f}") 
  for angle in angles:
    a = angle.split('-')
    comb = '-'.join([str(s) for s in a])
    val = np.array(angles[comb]).mean()
    if a[1] in tricentertypes:
      print(f"anglep {a[0]}  {a[1]}  {a[2]}  0.00 {val:10.4f}")
    else:
      print(f"angle  {a[0]}  {a[1]}  {a[2]}  0.00 {val:10.4f}")
  
  for angle in angles:
    a = angle.split('-')
    print(f"strbnd {a[0]}  {a[1]}  {a[2]}  0.00 0.00") 
  tmp = []
  for bond in bonds:
    b = bond.split('-')
    if (b[0] in tricentertypes) and ([b[1], b[0]] not in tmp):
      print(f"opbend {b[1]}  {b[0]}  0  0  0.00  0.00")
      tmp.append([b[1], b[0]])
    if b[1] in tricentertypes and ([b[0], b[1]] not in tmp):
      print(f"opbend {b[0]}  {b[1]}  0  0  0.00  0.00")
      tmp.append([b[0], b[1]])
  for tor in torsions:
    [a,b,c,d] = tor
    print(f"torsion {a}  {b}  {c}  {d}  0.00        0.0   1   0.00        180.0   2    0.00   0.0   3")
  return


if __name__ == "__main__":
  usage = "python lvalence_init.py txyz"
  if len(sys.argv) != 2: 
    sys.exit(f"Usage of this program: {usage}")
  txyzfile = sys.argv[1]
  if not os.path.isfile(txyzfile):
    sys.exit(f"{txyzfile} does not exist")
  print_initial_parameters(txyzfile)