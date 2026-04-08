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


info = ''' convert PDB structure to TINKER xyz 
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
import shutil
import subprocess
import numpy as np
import concurrent.futures


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


def _check_tool(name):
  ''' Check that an external tool is available on PATH. '''
  if shutil.which(name) is None:
    sys.exit(f"Error: '{name}' not found on PATH. Please install it first.")


def _get_base_name(template_name):
  ''' Strip trailing _N numeric suffix from a template name.

  Examples: DA_1 -> DA, HIS_2 -> HIS, HISC_3 -> HISC, ALA -> ALA
  '''
  parts = template_name.rsplit('_', 1)
  if len(parts) == 2 and parts[1].isdigit():
    return parts[0]
  return template_name


def _matches_residue(template_name, resname):
  ''' Check if a template name belongs to a given residue.

  Handles both naming conventions:
    - Amino acid terminal variants: ALA matches ALA, ALAC, ALAN
      (single uppercase letter appended for C/N-terminal)
    - Numbered parts: DA matches DA_1, DA_2, DA_3
      (_N suffix for multi-part templates)
    - Combined: HIS matches HIS, HISC, HISN, HIS_1, HISC_1, etc.
  '''
  if template_name == resname:
    return True
  base = _get_base_name(template_name)
  if base == resname:
    return True
  # Amino acid terminal variants: base name has exactly one extra uppercase character
  if base.startswith(resname) and len(base) == len(resname) + 1 and base[-1].isupper():
    return True
  return False


def _find_templates(resname, template_names, template_dir):
  ''' Find all template file paths for a residue from a list of template names. '''
  matched = [name for name in template_names if _matches_residue(name, resname)]
  return sorted(os.path.join(template_dir, m + ".txyz") for m in matched)


''' write a pdb for obabel '''
def prepare(pdb):
  _check_tool("obabel")
  # Extract ATOM/HETATM lines
  with open(pdb) as fin, open(pdb + "_1", "w") as fout:
    for line in fin:
      if line.startswith("ATOM") or line.startswith("HETATM"):
        fout.write(line)

  lines = open(pdb + "_1").readlines()
  with open(pdb + "_2", "w") as f:
    for line in lines:
      atomline = line[12:16]
      if not (('MG' in atomline)):
        line = line[:76] + ' ' + atomline.split()[0][0] + "\n"
      else:
        line = line[:76] + ' ' + atomline.strip() + "\n"
      f.write(line)
  txyz_out = _replace_ext(pdb, '.txyz')
  _run(["obabel", "-ipdb", pdb, "-otxyz", "-O", txyz_out])
  return


def _load_ion_types(database, rootdir):
  ''' Read atom types from single-atom templates in the database.

  Returns a dict mapping residue name -> atom type string for all
  single-atom (ion) templates in both solvent and ion categories.
  '''
  ion_types = {}
  for category in ['ion', 'solvent']:
    for resname in database.get(category, []):
      tpath = os.path.join(rootdir, 'database.PDB2txyz', category, resname + '.txyz')
      if os.path.isfile(tpath):
        with open(tpath) as f:
          tlines = f.readlines()
        if len(tlines) >= 2:
          try:
            natom = int(tlines[0].split()[0])
          except (ValueError, IndexError):
            print(f"Warning: malformed template header in {tpath}, skipping")
            continue
          if natom == 1:
            parts = tlines[1].split()
            if len(parts) >= 6:
              ion_types[resname] = parts[5]
  return ion_types


''' split input pdb file into fragment pdbs '''
def splitpdb(pdb, database, rootdir):
  if not os.path.isfile(pdb):
    sys.exit(f"Error: PDB file not found: {pdb}")

  lines = open(pdb).readlines()
  pdb_lines = {}
  number_res = 0
  number_atm = 0
  pdbstrs = []
  pdbs = []
  tmp = []
  txyzstr = []
  solvent_ion_residues = set(database.get('solvent', []) + database.get('ion', []))
  ion_types = _load_ion_types(database, rootdir)
  for line in lines:
    number_atm += 1
    curresnm = line[17:21].strip()
    if ("TER" not in line) and ("END" not in line) and ("REMARK" not in line) and ("CRYST" not in line) and ("MODEL" not in line) and (curresnm not in solvent_ion_residues):
      pdbstrs.append(line)
      curresid = line[22:26].strip()
      name_id = line[17:26]
      if name_id not in tmp:
        number_res += 1
        tmp.append(name_id)
        pdbname = "%04d"%number_res + f"_{curresnm.split()[0]}.pdb"
        pdbs.append(pdbname)

      if pdbname not in pdb_lines.keys():
        pdb_lines[pdbname] = [line]
      else:
        pdb_lines[pdbname] = pdb_lines[pdbname] + [line]

    # write solvent and ions if there are
    if curresnm in solvent_ion_residues:
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      atom = line[12:16].strip()
      if curresnm in ion_types:
        # Single-atom ion: look up type from template
        line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} {ion_types[curresnm]}\n"
        txyzstr.append(line_s)
      elif atom in ['OH2', 'OW', 'O']:
        line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 349 {number_atm+1} {number_atm+2}\n"
        txyzstr.append(line_s)
      elif atom in ['H1', 'HW1']:
        line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 350 {number_atm-1} \n"
        txyzstr.append(line_s)
      elif atom in ['H2', 'HW2']:
        line_s = f"{number_atm:>8d}{atom:>5s}{x:12.4f}{y:12.4f}{z:12.4f} 350 {number_atm-2} \n"
        txyzstr.append(line_s)
      else:
        sys.exit(f'Could not recognize atom {atom} in residue {curresnm} (line {number_atm})')

  for pdbname in pdb_lines:
    with open(pdbname, 'w') as f:
      pdbstr = pdb_lines[pdbname]
      for s in pdbstr:
        f.write(s)

  with open('pdblist', 'w') as f:
    for s in pdbs:
      f.write(s + '\n')

  txyzname = "solvent.txyz"
  natom = len(txyzstr)
  if natom > 0:
    with open(txyzname, 'w') as f:
      for s in txyzstr:
        f.write(s)
  return


''' read the template database '''
def readdatabase(rootdir):
  database = {}
  dbdir = os.path.join(rootdir, 'database.PDB2txyz')
  if not os.path.isdir(dbdir):
    sys.exit(f"Error: database directory not found: {dbdir}")
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


def _get_glycan_residues(database):
  ''' Derive glycan residue names from the database instead of hardcoding. '''
  return database.get('glycan', [])


''' convert xyz to tinker xyz '''
def pdbtxyz(pdb, database, rootdir):
  resname = pdb.split(".pdb")[0].split("_")[-1]

  if not os.path.isfile(pdb):
    print(f"Warning: PDB file not found: {pdb}")
    return

  # glycan is special — derive from database
  glycan_residues = _get_glycan_residues(database)
  isGlycan = resname in glycan_residues

  for key, value in database.items():
    template_dir = os.path.join(rootdir, 'database.PDB2txyz', key)
    templates = _find_templates(resname, value, template_dir)
    if not templates:
      continue

    t = _replace_ext(pdb, '.txyz')
    if not os.path.isfile(t):
      print(resname)
      pdb_lines = open(pdb).readlines()
      if len(pdb_lines) > 1:
        if isGlycan:
          cmd = [sys.executable, os.path.join(rootdir, "IP_MatchTXYZ_Glycan.py")] + templates + [pdb.split('.')[0]]
          print(f"Running: {' '.join(cmd)}")
          _run(cmd, check=False)
        else:
          cmd = [sys.executable, os.path.join(rootdir, "IP_MatchTXYZ.py"), "-t"] + templates + ["-d", pdb]
          _run(cmd, check=False)
      else:
        # ion
        if resname not in ['HOH']:
          line = pdb_lines[0]
          x = float(line[30:38])
          y = float(line[38:46])
          z = float(line[46:54])
          with open(_replace_ext(pdb, '.txyz'), 'w') as f:
            f.write(f"1 {pdb}\n")
            f.write(f"1 {line[12:16]}{x:10.3f}{y:10.3f}{z:10.3f} 0\n")
        # water
    else:
      pass
  return


''' check the correctness of pdb and txyz mapping'''
def check_pdb_xyz(pdblist, xyzlist):
  if not os.path.isfile(pdblist):
    sys.exit(f"Error: file not found: {pdblist}")
  if not os.path.isfile(xyzlist):
    sys.exit(f"Error: file not found: {xyzlist}")
  pdbs = list(np.loadtxt(pdblist, usecols=(0), dtype='str', unpack=True))
  xyzs = list(np.loadtxt(xyzlist, usecols=(0), dtype='str', unpack=True))
  npdb = len(pdbs)
  nxyz = len(xyzs)
  if npdb == nxyz:
    for pdb, xyz in zip(pdbs, xyzs):
      plines = open(pdb).readlines()
      xlines = open(xyz).readlines()
      if len(plines) != len(xlines)-1:
        sys.exit(f"Error: {pdb} not the same as {xyz}")
  return


''' generate the final txyz '''
def connect(txyz, txyzs):
  check_pdb_xyz('pdblist', 'txyzlist')
  fname = txyz + "_2"
  alltypes = []
  for t in txyzs:
    types = np.loadtxt(t, usecols=(5,), dtype="str", unpack=True, skiprows=1)
    if types.ndim == 0:
      alltypes += [str(types), ]
    else:
      alltypes += list(types)
  lines = open(txyz).readlines()

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


def _apply_charmm_fixups(pdb):
  ''' Apply CHARMM-GUI residue name fixups using in-memory replacement. '''
  ress_in_pdbs = {"CYT": "DC ", "GUA": "DG ", "ADE": "DA ", "THY": "DT "}
  with open(pdb) as f:
    content = f.read()
  for res, res_ in ress_in_pdbs.items():
    content = content.replace(res, res_)
  with open(pdb, 'w') as f:
    f.write(content)


def _pdbtxyz_wrapper(args):
  ''' Wrapper for pdbtxyz to unpack arguments for ProcessPoolExecutor. '''
  pdb, database, rootdir = args
  return pdbtxyz(pdb, database, rootdir)


if __name__ == "__main__":
  if len(sys.argv) != 3:
    sys.exit(info)

  pdb = sys.argv[1]

  if not os.path.isfile(pdb):
    sys.exit(f"Error: PDB file not found: {pdb}")

  # special names in CHARMM-GUI file
  _apply_charmm_fixups(pdb)

  mode = sys.argv[2].upper()
  rootdir = os.path.join(os.path.split(__file__)[0])
  database = readdatabase(rootdir)


  if 'A' in mode:
    print("=" * 60)
    print("Step A: Preparing PDB and splitting into residue fragments")
    print("=" * 60)
    print(f"  Input PDB: {pdb}")
    prepare(pdb)
    print(f"  Generated TXYZ: {_replace_ext(pdb, '.txyz')}")
    splitpdb(pdb + '_2', database, rootdir)
    residue_fragments = np.atleast_1d(np.loadtxt("pdblist", dtype='str'))
    print(f"  Split into {len(residue_fragments)} residue fragment(s)")
    print("Step A completed.\n")

  if 'B' in mode:
    print("=" * 60)
    print("Step B: Converting residue PDB fragments to TINKER XYZ")
    print("=" * 60)
    pdbs = np.loadtxt("pdblist", dtype='str')
    print(f"  Processing {len(pdbs)} residue fragment(s)...")
    jobs = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
      args_list = [(p, database, rootdir) for p in pdbs]
      results = [executor.submit(_pdbtxyz_wrapper, args) for args in args_list]
      for f in concurrent.futures.as_completed(results):
        jobs.append(f.result())

    # Check that all expected txyz_2 files were generated
    txyzs = [_replace_ext(p, '.txyz_2') for p in pdbs]
    missing = []
    for p, t in zip(pdbs, txyzs):
      if not os.path.isfile(t):
        resname = str(p).split(".pdb")[0].split("_")[-1]
        missing.append((str(p), t, resname))

    if missing:
      dbdir = os.path.join(rootdir, 'database.PDB2txyz')
      print("\n" + "!" * 60)
      print("WARNING: Some residue TXYZ files were NOT generated.")
      print("!" * 60)
      print(f"  {len(missing)} of {len(pdbs)} residue(s) failed:\n")
      for pdb_frag, txyz_frag, resname in missing:
        print(f"    - {pdb_frag}  (residue: {resname})")
        print(f"      Expected output: {txyz_frag}")
      print("\n  To add a missing residue to the database:")
      print(f"    1. Create a TXYZ template file named <RESNAME>.txyz")
      print(f"    2. Place it in the appropriate subdirectory under:")
      print(f"         {dbdir}")
      categories = sorted(os.listdir(dbdir)) if os.path.isdir(dbdir) else []
      print(f"       Available categories: {', '.join(categories)}")
      print(f"    3. The template must follow TINKER XYZ format:")
      print(f"         Line 1: <natom> <title>")
      print(f"         Lines 2+: <index> <atom_name> <x> <y> <z> <type> <connected_indices...>")
      print(f"    4. Re-run Step B after adding the template.\n")

    with open('txyzlist', 'w') as f:
      for s in txyzs:
        f.write(s + '\n')
    n_success = len(pdbs) - len(missing)
    print(f"  Successfully generated {n_success} of {len(pdbs)} TXYZ file(s).")
    print("Step B completed.\n")

  if 'C' in mode:
    print("=" * 60)
    print("Step C: Assembling final TINKER XYZ file")
    print("=" * 60)
    _check_tool("obabel")
    if not os.path.isfile('bio.pdb'):
      txyz = _replace_ext(pdb, '.txyz')
    else:
      txyz = 'bio.txyz'
      pdb = 'bio.pdb'
    print(f"  Converting {pdb} -> {txyz} via obabel")
    _run(["obabel", "-ipdb", pdb, "-otxyz", "-O", txyz])
    txyzs = np.loadtxt("txyzlist",dtype='str')
    print(f"  Connecting {len(txyzs)} fragment TXYZ file(s)...")
    connect(txyz,txyzs)
    if os.path.isfile('solvent.txyz'):
      print("  Including solvent atoms from solvent.txyz")
      with open(txyz + "_2") as src, open("solvent.txyz") as sol, open("final.txyz", "w") as out:
        out.write(src.read())
        out.write(sol.read())
    else:
      shutil.copy(txyz + "_2", "final.txyz")
    print(f"  Final output: final.txyz")
    print("Step C completed.\n")
