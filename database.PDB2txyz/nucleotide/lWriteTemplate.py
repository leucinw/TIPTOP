''' write template xyz file from amoebabio18.prm for amino acid '''

''' usage: python lWriteTemplate.py '''

import math
import os

# Reuse internal phosphodiester atom types for terminal monophosphate atoms
# (amoebabio18.prm assigns -1 to the dedicated 5'/3'-Monophosphate biotypes,
# so the closest parameterized stand-ins are the phosphodiester P/OP).
TYPE_P  = 343
TYPE_OP = 344


def _parse_atom_line(line):
  dd = line.split()
  return {
    "idx":   int(dd[0]),
    "name":  dd[1],
    "x":     float(dd[2]),
    "y":     float(dd[3]),
    "z":     float(dd[4]),
    "type":  dd[5],
    "bonds": [int(b) for b in dd[6:]],
  }


def _format_atom(idx, element, x, y, z, atype, bonds):
  conn = '   '.join(str(b) for b in bonds)
  conn = ('   ' + conn) if bonds else ''
  return (f"{idx:>6d}{element:>3s}"
          f"{x:12.6f}{y:12.6f}{z:12.6f}"
          f"{str(atype):>6s}{conn}\n")


def _find_o3prime(atoms):
  ''' O3' is the only O in _2 with a single bond to a carbon (the C3'-O3'
      stub left after the cross-residue bond to next-P was stripped). '''
  for a in atoms:
    if not a["name"].startswith("O") or len(a["bonds"]) != 1:
      continue
    if atoms[a["bonds"][0] - 1]["name"].startswith("C"):
      return a
  raise RuntimeError("O3' not found in _2 template")


def _add_5prime_phosphate(residue):
  ''' Build _4 (5'-monophosphate) from _2 by adding a third non-bridging
      O on the existing terminal P. '''
  f2 = f"{residue}_2.txyz"
  f4 = f"{residue}_4.txyz"
  lines = open(f2).readlines()
  natom = int(lines[0].split()[0])
  atoms = [_parse_atom_line(lines[i]) for i in range(1, natom + 1)]
  p = atoms[0]   # P is always atom 1 in _2
  assert p["name"].startswith("P"), f"_2 atom 1 is not P in {residue}"
  new_idx = natom + 1
  # Place new OP opposite the centroid of existing P bonds (rough tetrahedral)
  cx = cy = cz = 0.0
  for b in p["bonds"]:
    nb = atoms[b - 1]
    cx += nb["x"] - p["x"]
    cy += nb["y"] - p["y"]
    cz += nb["z"] - p["z"]
  norm = math.sqrt(cx*cx + cy*cy + cz*cz) or 1.0
  ox = p["x"] - 1.5 * cx / norm
  oy = p["y"] - 1.5 * cy / norm
  oz = p["z"] - 1.5 * cz / norm
  with open(f4, "w") as f:
    f.write(f"{natom + 1:>6d} Template\n")
    for a in atoms:
      bonds = a["bonds"] + ([new_idx] if a["idx"] == p["idx"] else [])
      f.write(_format_atom(a["idx"], a["name"][0], a["x"], a["y"], a["z"],
                           a["type"], bonds))
    f.write(_format_atom(new_idx, "O", ox, oy, oz, TYPE_OP, [p["idx"]]))


def _add_3prime_phosphate(residue):
  ''' Build _5 (3'-monophosphate) from _2 by attaching a terminal PO3 group
      to the dangling O3'. '''
  f2 = f"{residue}_2.txyz"
  f5 = f"{residue}_5.txyz"
  lines = open(f2).readlines()
  natom = int(lines[0].split()[0])
  atoms = [_parse_atom_line(lines[i]) for i in range(1, natom + 1)]
  o3 = _find_o3prime(atoms)
  c3 = atoms[o3["bonds"][0] - 1]
  # Place P along the C3'-O3' bond axis, extended past O3'
  dx = o3["x"] - c3["x"]
  dy = o3["y"] - c3["y"]
  dz = o3["z"] - c3["z"]
  n = math.sqrt(dx*dx + dy*dy + dz*dz) or 1.0
  dx, dy, dz = dx / n, dy / n, dz / n
  px = o3["x"] + 1.6 * dx
  py = o3["y"] + 1.6 * dy
  pz = o3["z"] + 1.6 * dz
  # Two unit vectors perpendicular to the P-O3' axis for placing the 3 OPs
  ax = (1.0, 0.0, 0.0) if abs(dx) < 0.9 else (0.0, 1.0, 0.0)
  dot = ax[0]*dx + ax[1]*dy + ax[2]*dz
  e1 = (ax[0] - dot*dx, ax[1] - dot*dy, ax[2] - dot*dz)
  n1 = math.sqrt(sum(c*c for c in e1)) or 1.0
  e1 = (e1[0]/n1, e1[1]/n1, e1[2]/n1)
  e2 = (dy*e1[2] - dz*e1[1], dz*e1[0] - dx*e1[2], dx*e1[1] - dy*e1[0])
  # Three OP directions 120° apart in the plane perpendicular to P-O3'
  cos120, sin120 = -0.5, 0.8660254
  op_dirs = [
    e1,
    (cos120*e1[0] + sin120*e2[0], cos120*e1[1] + sin120*e2[1], cos120*e1[2] + sin120*e2[2]),
    (cos120*e1[0] - sin120*e2[0], cos120*e1[1] - sin120*e2[1], cos120*e1[2] - sin120*e2[2]),
  ]
  p_idx = natom + 1
  op_idx = [natom + 2, natom + 3, natom + 4]
  with open(f5, "w") as f:
    f.write(f"{natom + 4:>6d} Template\n")
    for a in atoms:
      bonds = a["bonds"] + ([p_idx] if a["idx"] == o3["idx"] else [])
      f.write(_format_atom(a["idx"], a["name"][0], a["x"], a["y"], a["z"],
                           a["type"], bonds))
    f.write(_format_atom(p_idx, "P", px, py, pz, TYPE_P,
                         [o3["idx"]] + op_idx))
    for idx, d in zip(op_idx, op_dirs):
      f.write(_format_atom(idx, "O",
                           px + 1.5 * d[0], py + 1.5 * d[1], pz + 1.5 * d[2],
                           TYPE_OP, [p_idx]))


def getXYZ(residue):
  residue = residue.strip()
  cmd = f"sh nucleic.sh {residue}"
  os.system(cmd)
  tripep = residue + "3.xyz"
  lines = open(tripep).readlines()
  natom = int(lines[0].split()[0])
  n2 = int((natom+1)/3)
  n1 = n2 - 2
  n3 = n2 + 1
  f1 = f"{residue}_1.txyz"
  f2 = f"{residue}_2.txyz"
  f3 = f"{residue}_3.txyz"
  with open(f1, "w") as f:
    f.write(f"{n1:>6d} Template\n")
    for i in range(1, n1+1):
      dd = lines[i].split()
      conlist = []
      for n in dd[6:]:
        if int(n) <= n1:
          conlist.append(n)
      linestr = f"{dd[0]:>6s}{dd[1][0]:>3s}{float(dd[2]):12.6f}{float(dd[3]):12.6f}{float(dd[4]):12.6f}{dd[5]:>6s}" + '   ' + '   '.join(conlist)
      f.write(f"{linestr}\n")
  with open(f2, "w") as f:
    f.write(f"{n2:>6d} Template\n")
    for i in range(n1+1, n2+n1+1):
      dd = lines[i].split()
      conlist = []
      for n in dd[6:]:
        if n1 < int(n) <= n1+n2:
          conlist.append(str(int(n)-n1))
      linestr = f"{int(dd[0])-n1:>6d}{dd[1][0]:>3s}{float(dd[2]):12.6f}{float(dd[3]):12.6f}{float(dd[4]):12.6f}{dd[5]:>6s}" + '   ' + '   '.join(conlist)
      f.write(f"{linestr}\n")
  with open(f3, "w") as f:
    f.write(f"{n3:>6d} Template\n")
    for i in range(n1+n2+1, n2+n3+n1+1):
      dd = lines[i].split()
      conlist = []
      for n in dd[6:]:
        if n1+n2 < int(n) <= n1+n2+n3:
          conlist.append(str(int(n)-n1-n2))
      linestr = f"{int(dd[0])-n1-n2:>6d}{dd[1][0]:>3s}{float(dd[2]):12.6f}{float(dd[3]):12.6f}{float(dd[4]):12.6f}{dd[5]:>6s}" + '   ' + '   '.join(conlist)
      f.write(f"{linestr}\n")
  _add_5prime_phosphate(residue)
  _add_3prime_phosphate(residue)
  os.system(f"rm {tripep}")
  return


def main():
  ## nucleic acids in amoebabio18.prm
  residues = \
  ['DA ', 'DT ', 'DC ', 'DG ', 'A ', 'C ', 'U ', 'G '] 
  for residue in residues[0:]:
    getXYZ(residue)
  return

if __name__ == "__main__":
  main()
