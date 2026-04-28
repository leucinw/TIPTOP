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
#        Chengwen Liu               #
#    liuchw2010@gmail.com           #
#   University of Texas at Austin   #
#===================================

info = ''' convert PDB structure to TINKER xyz

  1. Amino acids and nucleotides are from amoebabio18.prm
  2. Glycans are generated using Poltype
  3. Lipids are taken from literature and generated using Poltype
  4. Can add more support by adding `residue` in database directory

  Usage: python IP_PDB2txyz.py your.pdb mode [--residue LIG.txyz MOL.txyz ...]
  ~ mode can be chosen from [A,B,C]
  ~ suggest run a,b,c sequencially, to make sure every step is good
  ~ can also run ABC at the same time if you are confident
  ~ use --residue (-r) to supply TXYZ templates for non-standard residues;
    the residue name is taken from the filename stem (LIG.txyz -> residue LIG)

  !!! your.pdb MUST have hydrogens and protonation state correctly!!!
'''

import os
import sys
import shutil
import subprocess
import concurrent.futures

# ---------------------------------------------------------------------------
#  Solvent / ion residue name handling
#  Atom types are read from database.PDB2txyz/solvent/ template files.
# ---------------------------------------------------------------------------

# Maps variant residue names to the canonical template name in the database.
# Add entries here when the PDB uses a name that differs from the template file.
SOLVENT_ALIASES = {
    'WAT':  'HOH', 'TIP3': 'HOH', 'TIP4': 'HOH', 'SPC': 'HOH', 'SOL': 'HOH',
    'SOD':  'NA',
    'CL':   'CLA',
    'POT':  'K',
    'MG2':  'MG',
}

# Residue names treated as solvent / ion (skipped in biomolecule pipeline).
# Includes both canonical template names and their aliases above.
SOLVENT_RESIDUE_NAMES = {
    'HOH', 'WAT', 'TIP3', 'TIP4', 'SPC', 'SOL',   # water
    'NA', 'SOD',                                     # sodium
    'CLA', 'CL',                                     # chloride
    'K', 'POT',                                      # potassium
    'MG', 'MG2',                                     # magnesium
}


def _read_lines(path):
    ''' Read non-empty stripped lines from a text file. '''
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


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

    Handles multiple naming conventions:
      - Numbered parts:  DA  matches DA_1, DA_2, DA_3
      - AA terminals:    ALA matches ALAC, ALAN  (single uppercase letter suffix)
      - Future-proof for digit suffixes (e.g. DA5, DA3) if ever added
    '''
    if template_name == resname:
        return True
    base = _get_base_name(template_name)
    if base == resname:
        return True
    # Terminal / variant suffix: one extra alphanumeric character
    #   - uppercase letter for amino acids (C = C-term, N = N-term)
    #   - digit for nucleotides (5 = 5'-term, 3 = 3'-term)
    if (base.startswith(resname)
            and len(base) == len(resname) + 1
            and base[-1].isalnum()):
        return True
    return False


def _find_templates(resname, template_names, template_dir):
    ''' Find all template file paths for a residue. '''
    matched = [name for name in template_names if _matches_residue(name, resname)]
    return sorted(os.path.join(template_dir, m + ".txyz") for m in matched)


def _is_solvent_or_ion(resname):
    ''' Check whether a residue name belongs to solvent or ion. '''
    return resname in SOLVENT_RESIDUE_NAMES


def _atom_element(atom_name):
    ''' Return the element symbol from a PDB atom name (e.g. OH2 -> O, HW1 -> H). '''
    for ch in atom_name:
        if ch.isalpha():
            return ch.upper()
    return atom_name[0].upper()


def _read_txyz_templates(directory):
    ''' Read atom types and connectivity from all .txyz files in a directory.

    Returns dict: resname -> list of (atom_name, atom_type, bonded_indices)
    where bonded_indices are 0-based indices into the template atom list.
    '''
    templates = {}
    if not os.path.isdir(directory):
        return templates
    for fname in os.listdir(directory):
        if not fname.endswith('.txyz'):
            continue
        resname = fname[:-5]
        with open(os.path.join(directory, fname)) as fh:
            lines = fh.readlines()
        atoms = []
        for line in lines[1:]:
            parts = line.split()
            if len(parts) < 6:
                continue
            name = parts[1]
            atype = int(parts[5])
            bonds = [int(b) - 1 for b in parts[6:]]  # convert to 0-based
            atoms.append((name, atype, bonds))
        if atoms:
            templates[resname] = atoms
    return templates


def _read_solvent_templates(rootdir):
    ''' Read atom types and connectivity from database.PDB2txyz/solvent/ templates. '''
    solvent_dir = os.path.join(rootdir, 'database.PDB2txyz', 'solvent')
    return _read_txyz_templates(solvent_dir)


def _read_ion_templates(rootdir):
    ''' Read atom types and connectivity from database.PDB2txyz/ion/ templates. '''
    ion_dir = os.path.join(rootdir, 'database.PDB2txyz', 'ion')
    return _read_txyz_templates(ion_dir)


# -----------------------------------------------------------------------
#  Step A helpers
# -----------------------------------------------------------------------

def _split_pdb_biomol_solvent(pdb):
    ''' Split a PDB into biomolecule-only and solvent/ion-only PDB files.

    Returns:
        (biomol_pdb, solvent_pdb)  -- paths written; solvent_pdb is None
                                      if no solvent/ion lines found.
    '''
    biomol_path = _replace_ext(pdb, '.biomol.pdb')
    solvent_path = _replace_ext(pdb, '.solvent.pdb')

    biomol_lines = []
    solvent_lines = []

    with open(pdb) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            resname = line[17:21].strip()
            if _is_solvent_or_ion(resname):
                solvent_lines.append(line)
            else:
                biomol_lines.append(line)

    # Fix element column for obabel (column 76-77)
    fixed_biomol = []
    for line in biomol_lines:
        atomname = line[12:16]
        if 'MG' not in atomname:
            line = line[:76] + ' ' + atomname.split()[0][0] + "\n"
        else:
            line = line[:76] + ' ' + atomname.strip() + "\n"
        fixed_biomol.append(line)

    with open(biomol_path, 'w') as f:
        for line in fixed_biomol:
            f.write(line)

    solvent_pdb = None
    if solvent_lines:
        with open(solvent_path, 'w') as f:
            for line in solvent_lines:
                f.write(line)
        solvent_pdb = solvent_path

    return biomol_path, solvent_pdb


def _build_solvent_txyz(solvent_pdb, solvent_templates, ion_templates):
    ''' Build a solvent/ion TXYZ using atom types and connectivity from
    database templates (database.PDB2txyz/solvent/ and database.PDB2txyz/ion/).

    PDB atoms are matched to template atoms by element symbol.
    Template connectivity is preserved: O bonds to H atoms for water,
    ions have no bonds.  Returns path to solvent.txyz, or None if no atoms.
    '''
    if solvent_pdb is None:
        return None

    with open(solvent_pdb) as f:
        pdb_lines = [l for l in f if l.startswith("ATOM") or l.startswith("HETATM")]

    if not pdb_lines:
        return None

    # Group PDB atoms by residue, preserving order of first appearance.
    residue_groups = {}   # (chain, resid, resname) -> [(atom_name, x, y, z)]
    residue_order = []

    for line in pdb_lines:
        atom_name = line[12:16].strip()
        resname   = line[17:21].strip()
        chain     = line[21:22]
        resid     = line[22:27].strip()
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        key = (chain, resid, resname)
        if key not in residue_groups:
            residue_groups[key] = []
            residue_order.append(key)
        residue_groups[key].append((atom_name, x, y, z))

    txyz_lines = []
    idx = 1  # 1-based sequential atom index

    for key in residue_order:
        chain, resid, resname = key
        pdb_atoms = residue_groups[key]

        # Resolve to the canonical template name.
        canonical = SOLVENT_ALIASES.get(resname, resname)
        template = solvent_templates.get(canonical) or ion_templates.get(canonical)
        if template is None:
            sys.exit(
                f"Error: no template found for solvent/ion residue '{resname}' "
                f"(canonical: '{canonical}'). "
                f"Add '{canonical}.txyz' to database.PDB2txyz/solvent/ or "
                f"database.PDB2txyz/ion/."
            )

        n_template = len(template)

        if n_template == 1:
            # ---- Ion: single atom, no bonds ----
            if len(pdb_atoms) != 1:
                sys.exit(f"Error: ion residue {key} has {len(pdb_atoms)} atoms, expected 1")
            atom_name, x, y, z = pdb_atoms[0]
            _, atype, _ = template[0]
            txyz_lines.append(
                f"{idx:>8d}{atom_name:>5s}{x:12.4f}{y:12.4f}{z:12.4f}{atype:>5d}\n"
            )
            idx += 1

        else:
            # ---- Multi-atom molecule (e.g. water): match by element ----
            # Build per-element queues from PDB atoms.
            elem_queue = {}
            for a, x, y, z in pdb_atoms:
                e = _atom_element(a)
                elem_queue.setdefault(e, []).append((a, x, y, z))

            # Assign PDB atoms to template slots in template order.
            elem_cursors = {}
            assigned = []  # (pdb_atom_name, x, y, z, atype, bonds_0based)
            for tname, atype, bonds in template:
                e = _atom_element(tname)
                i = elem_cursors.get(e, 0)
                pool = elem_queue.get(e, [])
                if i >= len(pool):
                    sys.exit(
                        f"Error: residue {key} has too few '{e}' atoms "
                        f"to match template '{canonical}'"
                    )
                aname, x, y, z = pool[i]
                assigned.append((aname, x, y, z, atype, bonds))
                elem_cursors[e] = i + 1

            # Write TXYZ lines with connectivity remapped to absolute indices.
            start_idx = idx
            for i, (aname, x, y, z, atype, bonds) in enumerate(assigned):
                atom_idx = start_idx + i
                abs_bonds = [start_idx + b for b in bonds]
                bond_str = ''.join(f"{b:>8d}" for b in abs_bonds)
                txyz_lines.append(
                    f"{atom_idx:>8d}{aname:>5s}{x:12.4f}{y:12.4f}{z:12.4f}"
                    f"{atype:>5d}{bond_str}\n"
                )
            idx += len(assigned)

    natom = idx - 1
    outpath = "solvent.txyz"
    with open(outpath, 'w') as f:
        f.write(f"{natom:>8d} Solvent and ions\n")
        for line in txyz_lines:
            f.write(line)

    return outpath


def prepare(pdb, database, rootdir):
    ''' Step A: split PDB into biomolecule vs solvent/ion, convert
    biomolecule-only PDB to TXYZ via obabel (no spurious bonds),
    build solvent TXYZ from database templates, and split biomolecule
    PDB into per-residue fragments.
    '''
    _check_tool("obabel")

    # 1. Split into biomolecule-only and solvent-only PDBs
    biomol_pdb, solvent_pdb = _split_pdb_biomol_solvent(pdb)
    n_biomol = sum(1 for _ in open(biomol_pdb))
    n_solvent = sum(1 for _ in open(solvent_pdb)) if solvent_pdb else 0
    print(f"  Biomolecule atoms: {n_biomol}")
    print(f"  Solvent/ion atoms: {n_solvent}")

    # 2. Convert biomolecule-only PDB to TXYZ via obabel
    #    This avoids spurious bonds between water/ions and biomolecule
    biomol_txyz = _replace_ext(biomol_pdb, '.txyz')
    _run(["obabel", "-ipdb", biomol_pdb, "-otxyz", "-O", biomol_txyz])
    print(f"  Biomolecule TXYZ (obabel): {biomol_txyz}")

    # 3. Build solvent TXYZ using atom types from database templates
    solvent_templates = _read_solvent_templates(rootdir)
    ion_templates = _read_ion_templates(rootdir)
    solvent_txyz = _build_solvent_txyz(solvent_pdb, solvent_templates, ion_templates)
    if solvent_txyz:
        print(f"  Solvent TXYZ: {solvent_txyz}")

    # 4. Split biomolecule PDB into per-residue fragments
    splitpdb(biomol_pdb, database)

    return biomol_pdb, biomol_txyz, solvent_txyz


def splitpdb(pdb, database):
    ''' Split a biomolecule-only PDB into per-residue fragment PDB files.
    Writes pdblist with the fragment filenames.
    '''
    if not os.path.isfile(pdb):
        sys.exit(f"Error: PDB file not found: {pdb}")

    with open(pdb) as f:
        lines = [l for l in f if l.startswith("ATOM") or l.startswith("HETATM")]

    pdb_lines = {}
    number_res = 0
    pdbs = []
    seen_residues = []

    for line in lines:
        resname = line[17:21].strip()
        name_id = line[17:26]  # residue name + sequence number

        if name_id not in seen_residues:
            number_res += 1
            seen_residues.append(name_id)
            pdbname = f"{number_res:04d}_{resname}.pdb"
            pdbs.append(pdbname)

        pdbname = f"{number_res:04d}_{resname}.pdb"
        if pdbname not in pdb_lines:
            pdb_lines[pdbname] = []
        pdb_lines[pdbname].append(line)

    for pdbname, pdbstr in pdb_lines.items():
        with open(pdbname, 'w') as f:
            for s in pdbstr:
                f.write(s)

    with open('pdblist', 'w') as f:
        for s in pdbs:
            f.write(s + '\n')

    return


def readdatabase(rootdir):
    ''' Read the template database directory structure. '''
    database = {}
    dbdir = os.path.join(rootdir, 'database.PDB2txyz')
    if not os.path.isdir(dbdir):
        sys.exit(f"Error: database directory not found: {dbdir}")
    for ff in os.listdir(dbdir):
        dirname = os.path.join(dbdir, ff)
        if os.path.isdir(dirname):
            templates = []
            for f in os.listdir(dirname):
                if f.endswith(".txyz"):
                    templates.append(f.split(".txyz")[0])
            if templates:
                database[ff] = templates
    return database


def _get_glycan_residues(database):
    return database.get('glycan', [])


def pdbtxyz(pdb, database, rootdir, user_residues=None):
    ''' Convert a single residue PDB fragment to TINKER XYZ using templates.

    Looks up templates in the database first; if none match, falls back to
    user_residues (a dict of resname -> [txyz_path, ...] supplied via --residue).
    '''
    resname = pdb.split(".pdb")[0].split("_")[-1]

    if not os.path.isfile(pdb):
        print(f"Warning: PDB file not found: {pdb}")
        return

    glycan_residues = _get_glycan_residues(database)
    isGlycan = resname in glycan_residues

    def _run_match(templates):
        ''' Run IP_MatchTXYZ (or Glycan variant) for this fragment. '''
        t = _replace_ext(pdb, '.txyz')
        if os.path.isfile(t):
            return
        print(resname)
        with open(pdb) as f:
            pdb_lines = f.readlines()
        if len(pdb_lines) > 1:
            if isGlycan:
                cmd = ([sys.executable,
                        os.path.join(rootdir, "IP_MatchTXYZ_Glycan.py")]
                       + templates + [pdb.split('.')[0]])
                print(f"Running: {' '.join(cmd)}")
                _run(cmd, check=False)
            else:
                cmd = ([sys.executable,
                        os.path.join(rootdir, "IP_MatchTXYZ.py"),
                        "-t"] + templates + ["-d", pdb])
                _run(cmd, check=False)
        else:
            # single-atom residue fallback
            if resname not in ('HOH', 'WAT'):
                line = pdb_lines[0]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                with open(_replace_ext(pdb, '.txyz'), 'w') as fout:
                    fout.write(f"1 {pdb}\n")
                    fout.write(f"1 {line[12:16]}{x:10.3f}{y:10.3f}{z:10.3f} 0\n")

    # 1. Search the database (built-in residues).
    for key, value in database.items():
        template_dir = os.path.join(rootdir, 'database.PDB2txyz', key)
        templates = _find_templates(resname, value, template_dir)
        if templates:
            _run_match(templates)
            return

    # 2. Fall back to user-supplied templates (--residue flag).
    if user_residues:
        templates = user_residues.get(resname)
        if templates:
            _run_match(templates)


def _pdbtxyz_wrapper(args):
    pdb, database, rootdir, user_residues = args
    return pdbtxyz(pdb, database, rootdir, user_residues)


def check_pdb_xyz(pdblist, xyzlist):
    ''' Validate that each fragment PDB matches its TXYZ in atom count. '''
    pdbs = _read_lines(pdblist)
    xyzs = _read_lines(xyzlist)

    if len(pdbs) != len(xyzs):
        sys.exit(f"Error: pdblist has {len(pdbs)} entries but txyzlist has {len(xyzs)}")

    for pdb, xyz in zip(pdbs, xyzs):
        if not os.path.isfile(pdb) or not os.path.isfile(xyz):
            continue
        with open(pdb) as f:
            npdb = sum(1 for _ in f)
        with open(xyz) as f:
            nxyz = sum(1 for _ in f) - 1  # header line
        if npdb != nxyz:
            sys.exit(f"Error: {pdb} ({npdb} atoms) does not match {xyz} ({nxyz} atoms)")


def connect(biomol_txyz, fragment_txyzs):
    ''' Replace atom types in the obabel-generated biomolecule TXYZ with
    the correct types from the per-residue template-matched TXYZ files.

    Only operates on biomolecule atoms (no solvent/ions).
    '''
    check_pdb_xyz('pdblist', 'txyzlist')

    # Collect atom types from all fragment txyz files
    alltypes = []
    for t in fragment_txyzs:
        if not os.path.isfile(t):
            sys.exit(f"Error: fragment TXYZ not found: {t}")
        with open(t) as f:
            for i, line in enumerate(f):
                if i == 0:
                    continue  # skip header
                parts = line.split()
                if len(parts) >= 6:
                    alltypes.append(parts[5])

    with open(biomol_txyz) as f:
        lines = f.readlines()

    natom = int(lines[0].split()[0])
    if len(alltypes) != natom:
        sys.exit(f"Error: type count mismatch — {len(alltypes)} types from fragments "
                 f"vs {natom} atoms in {biomol_txyz}")

    outpath = biomol_txyz + "_typed"
    with open(outpath, "w") as f:
        f.write(f"{natom:>10d} Generated by TIPTOP/IP_PDB2txyz.py\n")
        for i in range(1, len(lines)):
            d = lines[i].split()
            if len(d) < 6:
                continue
            d[5] = alltypes[i - 1]
            constr = ''.join([f"{s:>10s}" for s in d[6:]])
            line = (f"{d[0]:>10s}{d[1]:>3s}{float(d[2]):12.4f}"
                    f"{float(d[3]):12.4f}{float(d[4]):12.4f}{d[5]:>5s}{constr}\n")
            f.write(line)

    return outpath


def assemble_final(biomol_typed_txyz, solvent_txyz):
    ''' Merge biomolecule and solvent TXYZ files into final.txyz with
    correct sequential atom numbering and shifted connectivity.
    '''
    # Read biomolecule typed txyz
    with open(biomol_typed_txyz) as f:
        biomol_lines = f.readlines()

    biomol_natom = int(biomol_lines[0].split()[0])

    # Read solvent txyz if present
    solvent_lines = []
    solvent_natom = 0
    if solvent_txyz and os.path.isfile(solvent_txyz):
        with open(solvent_txyz) as f:
            solvent_lines = f.readlines()
        solvent_natom = int(solvent_lines[0].split()[0])

    total_natom = biomol_natom + solvent_natom

    with open("final.txyz", "w") as f:
        # Header
        f.write(f"{total_natom:>10d} Generated by TIPTOP/IP_PDB2txyz.py\n")

        # Write biomolecule atoms (indices unchanged, already 1-based)
        for line in biomol_lines[1:]:
            if line.strip():
                f.write(line)

        # Write solvent atoms with indices shifted by biomol_natom
        if solvent_natom > 0:
            for line in solvent_lines[1:]:
                if not line.strip():
                    continue
                parts = line.split()
                # parts[0] = atom index, parts[1] = name, [2-4] = xyz,
                # parts[5] = type, parts[6:] = connectivity indices
                old_idx = int(parts[0])
                new_idx = old_idx + biomol_natom
                new_conn = []
                for c in parts[6:]:
                    new_conn.append(str(int(c) + biomol_natom))

                constr = ''.join([f"{c:>8s}" for c in new_conn])
                line_out = (f"{new_idx:>8d}{parts[1]:>5s}"
                            f"{float(parts[2]):12.4f}{float(parts[3]):12.4f}"
                            f"{float(parts[4]):12.4f}{int(parts[5]):>5d}{constr}\n")
                f.write(line_out)

    print(f"  Final output: final.txyz  ({total_natom} atoms: "
          f"{biomol_natom} biomolecule + {solvent_natom} solvent/ion)")
    return "final.txyz"


def _apply_charmm_fixups(pdb):
    ''' Apply CHARMM-GUI residue name fixups using in-memory replacement. '''
    ress_in_pdbs = {"CYT": "DC ", "GUA": "DG ", "ADE": "DA ", "THY": "DT "}

    with open(pdb) as f:
        content = f.read()

    for res, res_ in ress_in_pdbs.items():
        content = content.replace(res, res_)

    with open(pdb, 'w') as f:
        f.write(content)


# =======================================================================
#  Main
# =======================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert PDB structure to TINKER XYZ format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=info,
    )
    parser.add_argument("pdb",  help="Input PDB file (must have H and correct protonation)")
    parser.add_argument("mode", help="Steps to run: A, B, C, or combinations like ABC")
    parser.add_argument(
        "--residue", "-r",
        metavar="TXYZ",
        nargs="+",
        default=[],
        help=(
            "TXYZ template file(s) for non-standard residues not in the database. "
            "The residue name is taken from the filename stem "
            "(e.g. LIG.txyz -> residue LIG, LIG_1.txyz -> residue LIG). "
            "Multiple files may be given; files sharing the same base name "
            "are all used as templates for that residue."
        ),
    )
    args = parser.parse_args()

    pdb = args.pdb
    if not os.path.isfile(pdb):
        parser.error(f"PDB file not found: {pdb}")

    # Build user_residues: resname -> [absolute txyz paths]
    user_residues = {}
    for path in args.residue:
        if not os.path.isfile(path):
            parser.error(f"Residue TXYZ file not found: {path}")
        stem = os.path.splitext(os.path.basename(path))[0]
        base = _get_base_name(stem)          # strips _N suffix (LIG_1 -> LIG)
        user_residues.setdefault(base, []).append(os.path.abspath(path))

    if user_residues:
        print(f"User-supplied residues: {', '.join(sorted(user_residues))}")

    _apply_charmm_fixups(pdb)

    mode = args.mode.upper()
    rootdir = os.path.join(os.path.split(__file__)[0])
    database = readdatabase(rootdir)

    # These are written/read across steps — persist as files so steps
    # can be run independently.
    biomol_txyz_path = _replace_ext(pdb, '.biomol.txyz')
    solvent_txyz_path = "solvent.txyz"

    # -----------------------------------------------------------------
    if 'A' in mode:
        print("=" * 60)
        print("Step A: Split PDB -> biomolecule + solvent, convert via obabel")
        print("=" * 60)
        print(f"  Input PDB: {pdb}")

        biomol_pdb, biomol_txyz, solvent_txyz = prepare(pdb, database, rootdir)

        residue_fragments = _read_lines("pdblist")
        print(f"  Split into {len(residue_fragments)} residue fragment(s)")
        print("Step A completed.\n")

    # -----------------------------------------------------------------
    if 'B' in mode:
        print("=" * 60)
        print("Step B: Converting residue PDB fragments to TINKER XYZ")
        print("=" * 60)

        pdbs = _read_lines("pdblist")
        print(f"  Processing {len(pdbs)} residue fragment(s)...")

        with concurrent.futures.ProcessPoolExecutor() as executor:
            args_list = [(p, database, rootdir, user_residues) for p in pdbs]
            futures = [executor.submit(_pdbtxyz_wrapper, args) for args in args_list]
            for fut in concurrent.futures.as_completed(futures):
                fut.result()  # propagate exceptions

        # Check results
        txyzs = [_replace_ext(p, '.txyz_2') for p in pdbs]
        missing = []
        for p, t in zip(pdbs, txyzs):
            if not os.path.isfile(t):
                resname = p.split(".pdb")[0].split("_")[-1]
                missing.append((p, t, resname))

        if missing:
            dbdir = os.path.join(rootdir, 'database.PDB2txyz')
            print("\n" + "!" * 60)
            print("WARNING: Some residue TXYZ files were NOT generated.")
            print("!" * 60)
            print(f"  {len(missing)} of {len(pdbs)} residue(s) failed:\n")
            for pdb_frag, txyz_frag, resname in missing:
                print(f"    - {pdb_frag} (residue: {resname})")
                print(f"      Expected output: {txyz_frag}")
            print(f"\n  Options to resolve missing residues:")
            print(f"    1. Re-run with --residue <RESNAME>.txyz for each missing residue")
            print(f"    2. Add <RESNAME>.txyz to the database:")
            print(f"         {dbdir}/<category>/")
            categories = sorted(os.listdir(dbdir)) if os.path.isdir(dbdir) else []
            print(f"       Available categories: {', '.join(categories)}")
            print(f"  Then re-run Step B.\n")

        with open('txyzlist', 'w') as f:
            for s in txyzs:
                f.write(s + '\n')

        n_success = len(pdbs) - len(missing)
        print(f"  Successfully generated {n_success} of {len(pdbs)} TXYZ file(s).")
        print("Step B completed.\n")

    # -----------------------------------------------------------------
    if 'C' in mode:
        print("=" * 60)
        print("Step C: Assembling final TINKER XYZ file")
        print("=" * 60)

        # Use biomolecule-only txyz (no water/ion bonds from obabel)
        if not os.path.isfile(biomol_txyz_path):
            sys.exit(f"Error: biomolecule TXYZ not found: {biomol_txyz_path}\n"
                     f"  Did you run Step A first?")

        fragment_txyzs = _read_lines("txyzlist")
        print(f"  Connecting {len(fragment_txyzs)} fragment TXYZ file(s)...")

        biomol_typed = connect(biomol_txyz_path, fragment_txyzs)

        # Merge biomolecule + solvent with correct renumbering
        solvent = solvent_txyz_path if os.path.isfile(solvent_txyz_path) else None
        assemble_final(biomol_typed, solvent)

        print("Step C completed.\n")
