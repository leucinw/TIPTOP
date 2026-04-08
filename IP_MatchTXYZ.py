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
    ~ Uses graph isomorphism with element-aware matching (VF2 algorithm via networkx)
    ~ make sure you can run obabel command before you use this program
'''

import os
import sys
import shutil
import subprocess
import argparse
import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism as iso


# Two-letter element symbols that may appear in TXYZ atom names
_TWO_LETTER_ELEMENTS = {
  'He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
  'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
  'Rb','Sr','Zr','Nb','Mo','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb',
  'Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy',
  'Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt','Au','Hg',
  'Tl','Pb','Bi',
}


def _run(cmd, check=True):
  ''' Run a command using subprocess; cmd is a list of arguments. '''
  result = subprocess.run(cmd, capture_output=True, text=True)
  if check and result.returncode != 0:
    sys.exit(f"Command failed: {' '.join(str(c) for c in cmd)}\n{result.stderr}")
  return result


def _replace_ext(filepath, new_ext):
  ''' Replace the file extension, e.g. .pdb -> .txyz '''
  base, _ = os.path.splitext(filepath)
  return base + new_ext


def get_element(atom_name):
  ''' Extract the element symbol from a TXYZ/PDB atom name.

  Handles single-letter elements (C, N, O, H, S, P, etc.) and
  two-letter elements (Cl, Br, Mg, Zn, Ca, Fe, etc.).
  TXYZ files may write elements in ALL CAPS (CA, MG, ZN) or mixed
  case (Ca, Mg, Zn) — both forms are recognized.
  '''
  name = atom_name.strip()
  if not name:
    return ''
  # Check for two-letter element: try title-case normalization
  if len(name) >= 2 and name[0].isalpha() and name[1].isalpha():
    candidate = name[0].upper() + name[1].lower()
    if candidate in _TWO_LETTER_ELEMENTS:
      return candidate
  return name[0].upper()


def read_txyz(filepath):
  ''' Read a TINKER .txyz file.

  Returns:
    natom:       number of atoms
    atoms:       list of atom name strings
    coords:      list of [x, y, z] floats
    types:       list of atom type strings
    connections: list of lists of connected atom index strings (1-based)
  '''
  if not os.path.isfile(filepath):
    sys.exit(f"Error: TXYZ file not found: {filepath}")

  with open(filepath) as f:
    lines = f.readlines()

  natom = int(lines[0].split()[0])
  atoms = []
  coords = []
  types = []
  connections = []

  # Atom lines are the last `natom` lines (handles title/header on line 1)
  atom_lines = lines[1:natom + 1] if len(lines) >= natom + 1 else lines[1:]

  for line in atom_lines:
    data = line.split()
    if len(data) < 6:
      continue
    atoms.append(data[1])
    coords.append([float(data[2]), float(data[3]), float(data[4])])
    types.append(data[5])
    connections.append(data[6:])

  return natom, atoms, coords, types, connections


def build_molecular_graph(atoms, connections):
  ''' Build a networkx Graph with element labels as node attributes.

  Nodes are integer-indexed (1-based, matching TXYZ convention).
  Each node has an 'element' attribute derived from the atom name,
  used during isomorphism matching to ensure chemically valid mappings.

  Args:
    atoms:       list of atom name strings (index i → atom i+1)
    connections: list of lists of connected atom index strings (1-based)

  Returns:
    nx.Graph with 'element' node attributes
  '''
  G = nx.Graph()

  for i, (atom, conns) in enumerate(zip(atoms, connections), start=1):
    element = get_element(atom)
    G.add_node(i, element=element)
    for c in conns:
      j = int(c)
      if j > 0:
        G.add_edge(i, j)

  return G


def find_isomorphic_mapping(g_template, g_target):
  ''' Find an atom mapping from target to template using VF2 with element matching.

  Uses element labels as constraints so that only atoms of the same element
  can be mapped to each other (C→C, N→N, etc.). This prevents incorrect
  mappings in molecules where different elements share similar connectivity.

  Args:
    g_template: nx.Graph of the template molecule (with 'element' attributes)
    g_target:   nx.Graph of the target molecule (with 'element' attributes)

  Returns:
    (success, mapping) where mapping is {target_node: template_node}
  '''
  node_match = iso.categorical_node_match('element', '')
  GM = iso.GraphMatcher(g_target, g_template, node_match=node_match)

  if GM.is_isomorphic():
    return True, GM.mapping

  return False, {}


def convert_to_txyz(input_file):
  ''' Convert input file to .txyz format using obabel if needed.

  Supports .xyz, .pdb, .sdf input formats. If the file is already
  in .txyz format (or is a .xyz file that matches TXYZ layout), it
  is returned as-is.

  Returns:
    Path to the .txyz file.
  '''
  ext = os.path.splitext(input_file)[1].lower()

  if ext == ".txyz":
    return input_file

  if shutil.which("obabel") is None:
    sys.exit("Error: 'obabel' not found on PATH. Please install Open Babel first.")

  txyz_file = _replace_ext(input_file, ".txyz")
  format_map = {".xyz": "xyz", ".pdb": "pdb", ".sdf": "sdf", ".mol2": "mol2"}

  if ext == ".xyz":
    # Check if it's already a TXYZ (atom_count == num_data_lines)
    with open(input_file) as f:
      lines = f.readlines()
    if int(lines[0].split()[0]) == len(lines) - 1:
      # Already TXYZ layout — use directly
      return input_file

  if ext in format_map:
    _run(["obabel", f"-i{format_map[ext]}", input_file, "-otxyz", "-O", txyz_file])
    return txyz_file

  sys.exit(f"Error: unsupported file format: {ext}")


def match_and_assign_types(template_file, target_file, output_file=None):
  ''' Assign template atom types to a target molecule via graph isomorphism.

  Builds molecular graphs from both files with element labels, then uses
  the VF2 algorithm (with element-aware node matching) to find a valid
  isomorphic mapping. The template's AMOEBA atom types are transferred
  to the target's coordinates and connectivity.

  Args:
    template_file: path to template .txyz (with known atom types)
    target_file:   path to target .txyz (types will be assigned)
    output_file:   output path (default: target_file + "_2")

  Returns:
    True if matching succeeded, False otherwise
  '''
  if output_file is None:
    output_file = target_file + "_2"

  # Read both files
  _, t_atoms, t_coords, t_types, t_conns = read_txyz(template_file)
  _, d_atoms, d_coords, d_types, d_conns = read_txyz(target_file)

  if len(t_atoms) != len(d_atoms):
    print(f"  Skipping template {template_file}: atom count mismatch "
          f"(template={len(t_atoms)}, target={len(d_atoms)})")
    return False

  # Build graphs with element attributes
  g_template = build_molecular_graph(t_atoms, t_conns)
  g_target = build_molecular_graph(d_atoms, d_conns)

  # Find isomorphic mapping: target_node → template_node
  success, mapping = find_isomorphic_mapping(g_template, g_target)

  if not success:
    print(f"  Warning: no element-consistent isomorphism found between "
          f"{template_file} and {target_file}")
    return False

  # Transfer types using the mapping
  new_types = list(d_types)
  for target_idx, template_idx in mapping.items():
    new_types[target_idx - 1] = t_types[template_idx - 1]

  # Write output in standard TXYZ format
  natom = len(d_atoms)
  with open(output_file, 'w') as f:
    f.write("%3s\n" % natom)
    for i in range(natom):
      conn_str = '  '.join(d_conns[i])
      if conn_str:
        conn_str = '   ' + conn_str
      f.write("%3s%3s%12.6f%12.6f%12.6f  %s%s\n" % (
        i + 1, d_atoms[i],
        d_coords[i][0], d_coords[i][1], d_coords[i][2],
        new_types[i], conn_str))

  return True


# --- Backward-compatible aliases for callers that import these ---
def readTXYZ(filepath):
  ''' Backward-compatible wrapper around read_txyz.
  Returns (atoms, coords, order, types, connections) to match old interface.
  '''
  _, atoms, coords, types, connections = read_txyz(filepath)
  order = [str(i + 1) for i in range(len(atoms))]
  return atoms, coords, order, types, connections

def txyz2graph(txyz_file):
  ''' Build a molecular graph from a .txyz file (backward-compatible). '''
  _, atoms, _, _, connections = read_txyz(txyz_file)
  return build_molecular_graph(atoms, connections)

def matchgraphs(G1, G2):
  ''' Match two graphs with element-aware isomorphism (backward-compatible). '''
  return list(find_isomorphic_mapping(G1, G2))


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Assign atom types from a template TXYZ to a target structure "
                "using element-aware graph isomorphism (VF2 algorithm).")
  parser.add_argument('-t', dest='template', nargs='+',
                      help="Template txyz file(s)", required=True)
  parser.add_argument('-d', dest='dealwith',
                      help="File to deal-with, can be .xyz/.txyz/.pdb/.sdf",
                      required=True)
  parser.add_argument('-s', dest='sortatom',
                      help="Sort atoms according to template file (not yet implemented)",
                      default=False, type=bool)
  args = vars(parser.parse_args())
  templates = args["template"]
  dealwith = args["dealwith"]

  if not os.path.isfile(dealwith):
    sys.exit(f"Error: input file not found: {dealwith}")

  # Convert to .txyz if necessary
  dealwith = convert_to_txyz(dealwith)

  # Try each template until one matches
  fname = dealwith + "_2"
  matched = False
  for template in templates:
    if not os.path.isfile(template):
      print(f"  Warning: template file not found: {template}")
      continue
    if match_and_assign_types(template, dealwith, fname):
      matched = True
      break

  if not matched:
    print(f"  Error: no matching template found for {dealwith}")
    sys.exit(1)
