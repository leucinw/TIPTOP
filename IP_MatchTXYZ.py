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
    ~ Assign the atom types of the template.txyz file to the dealwith txyz/xyz/pdb/sdf file.
    ~ Especially suitable for the case that their atoms are in different order.
    ~ It uses connectivity information only to match the structures
    ~ make sure you can run obabel command before you use this program
'''

import os
import sys
import shutil
import subprocess
import argparse
import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso


def _run(cmd, check=True):
  ''' Run a command using subprocess; cmd is a list of arguments. '''
  result = subprocess.run(cmd, capture_output=True, text=True)
  if check and result.returncode != 0:
    sys.exit(f"Command failed: {' '.join(cmd)}\n{result.stderr}")
  return result


def _replace_ext(filepath, new_ext):
  ''' Replace the file extension, e.g. .pdb -> .txyz '''
  base, _ = os.path.splitext(filepath)
  return base + new_ext


def readTXYZ(TXYZ):
  if not os.path.isfile(TXYZ):
    sys.exit(f"Error: TXYZ file not found: {TXYZ}")
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

def txyz2graph(txyz):
  G = nx.Graph()
  nodes = []
  edges = []
  lines = open(txyz).readlines()
  natom = int(lines[0].split()[0])
  for line in lines[-natom:]:
    ss = line.split()
    if ss[0] not in nodes:
      nodes.append(ss[0])
    for s in ss[6:]:
      ds = [ss[0], s]
      if ds not in edges:
        edges.append(ds)

  G.add_nodes_from(nodes)
  G.add_edges_from(edges)
  return G

def matchgraphs(G1, G2):
  GM = iso.GraphMatcher(G1, G2)
  if not GM.is_isomorphic():
    return [False, {}]
  else:
    return [True, GM.mapping]

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-t', dest = 'template', nargs='+', help = "Template txyz file(s)", required=True)
  parser.add_argument('-d', dest = 'dealwith', help = "File to deal-with, can be .xyz/.txyz/.pdb/.sdf", required=True)
  parser.add_argument('-s', dest = 'sortatom', help = "sort atoms according to template file", default=False, type=bool)
  args = vars(parser.parse_args())
  templates = args["template"]
  dealwith  = args["dealwith"]
  sortatom  = args["sortatom"]

  if not os.path.isfile(dealwith):
    sys.exit(f"Error: input file not found: {dealwith}")

  if shutil.which("obabel") is None:
    sys.exit("Error: 'obabel' not found on PATH. Please install Open Babel first.")

  ext = os.path.splitext(dealwith)[1]
  if ext == ".xyz":
    lines = open(dealwith).readlines()
    if not (int(lines[0].split()[0]) == len(lines)-1):
      xyz = dealwith
      dealwith = _replace_ext(dealwith, ".txyz")
      _run(["obabel", "-ixyz", xyz, "-otxyz", "-O", dealwith])
  elif ext == ".pdb":
    xyz = dealwith
    dealwith = _replace_ext(dealwith, ".txyz")
    _run(["obabel", "-ipdb", xyz, "-otxyz", "-O", dealwith])
  elif ext == ".sdf":
    xyz = dealwith
    dealwith = _replace_ext(dealwith, ".txyz")
    _run(["obabel", "-isdf", xyz, "-otxyz", "-O", dealwith])

  fname = dealwith + "_2"
  for template in templates:
    atoms1, coord1, _, types1, connections1 =  readTXYZ(dealwith)
    atoms2, coord2, _, types2, connections2 =  readTXYZ(template)
    if len(atoms1) == len(atoms2):
      g1 = txyz2graph(template)
      g2 = txyz2graph(dealwith)
      match, newidx = matchgraphs(g2, g1)

      if match:
        with open(fname, 'w') as f:
          f.write("%3s\n"%len(atoms1))
          for i in range(len(newidx)):
            idx = int(newidx[str(i+1)]) - 1
            f.write("%3s%3s%12.6f%12.6f%12.6f  %s   %s\n"%(i+1,atoms1[i], coord1[i][0], coord1[i][1], coord1[i][2], types2[idx], '  '.join(connections1[i])))
