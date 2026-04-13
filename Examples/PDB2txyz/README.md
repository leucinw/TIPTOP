# IP_PDB2txyz.py — PDB to TINKER XYZ Converter

Convert a standard PDB structure into a TINKER XYZ file (`final.txyz`) with
correct AMOEBA atom types, ready for molecular dynamics or energy calculations
with the Tinker software package.

---

## Overview

`IP_PDB2txyz.py` performs a three-step pipeline (A → B → C) to produce a
fully typed TINKER XYZ file from a PDB structure:

| Step | What it does |
|------|--------------|
| **A** | Splits the PDB into a biomolecule part and a solvent/ion part; converts the biomolecule to a bare TXYZ via `obabel`; splits the biomolecule into per-residue fragment PDB files |
| **B** | Converts each residue fragment to a typed TXYZ by matching it against template files in `database.PDB2txyz/` |
| **C** | Assembles all typed fragment TXYZs back into the biomolecule TXYZ, merges with the solvent TXYZ, and writes `final.txyz` |

Steps can be run individually or all at once (`ABC`).

---

## Dependencies

| Tool | Purpose |
|------|---------|
| Python ≥ 3.8 | Runtime |
| [Open Babel](https://openbabel.org) (`obabel`) | PDB → bare TXYZ conversion (Step A) |
| [PDBFixer](https://github.com/openmm/pdbfixer) *(optional)* | Add missing hydrogens / fix PDB before running |

> **Important:** The input PDB **must already have all hydrogens** and the
> correct protonation state.  Use `pdb_fixer.py` (provided in this directory)
> or another tool to prepare the PDB first.

---

## Supported Residues

Atom types are drawn from the **`amoebabio18`** force field and are organised
into the database categories below.

| Category | Residues |
|----------|---------|
| `aminoacid` | Standard amino acids (ALA, ARG, ASN, …) plus N-/C-terminal variants (ALAN, ALAC, …) and capping group ACE |
| `nucleotide` | RNA (A, C, G, U) and DNA (DA, DC, DG, DT) with 5′/3′ terminal variants |
| `glycan` | Common *N*- and *O*-linked glycan residues (BGAL, AMAN, ANE5, …) |
| `aaconjugate` | Amino-acid conjugates (CGL, NGL, SGL, TGL, CYP, THP2) |
| `lipid` | POPC, POPE, POPI, POPS, CHL1 |
| `solvent` | HOH (water) and aliases WAT, TIP3, TIP4, SPC, SOL |
| `ion` | NA (SOD), CLA (CL), K (POT), MG (MG2) |

Non-standard residues not in the database can be supplied at run time with
the `--residue` flag (see below).

---

## Usage

```
python IP_PDB2txyz.py your.pdb <mode> [--residue LIG.txyz ...]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `your.pdb` | Input PDB file — **must contain hydrogens and correct protonation** |
| `mode` | One or more of `A`, `B`, `C` (case-insensitive). Run sequentially (`A` then `B` then `C`) or all at once (`ABC`) |
| `--residue / -r` | One or more TXYZ template files for non-standard residues not in the database. The residue name is taken from the filename stem (e.g. `LIG.txyz` → residue `LIG`) |

### Modes

```
A   Split PDB → biomolecule + solvent, convert biomolecule via obabel,
    generate per-residue fragment PDB files

B   Convert each residue fragment PDB to a typed TXYZ using database templates
    (runs fragments in parallel)

C   Assemble all fragment TXYZs → final.txyz (biomolecule + solvent/ions)

ABC Run all three steps in sequence
```

> **Tip:** Run A, B, and C separately when working with a new system so you
> can inspect intermediate results and catch problems early.

---

## Step-by-step Example

```bash
# 0. Prepare the PDB (add missing atoms and hydrogens)
python pdb_fixer.py          # produces 1EHZ_fixed.pdb (edit script for your PDB ID)

# 1. Step A — split and convert
python /path/to/TIPTOP/IP_PDB2txyz.py 1EHZ_fixed.pdb A

# 2. Step B — type each residue fragment
python /path/to/TIPTOP/IP_PDB2txyz.py 1EHZ_fixed.pdb B

# 3. Step C — assemble final.txyz
python /path/to/TIPTOP/IP_PDB2txyz.py 1EHZ_fixed.pdb C

# Or run everything at once (if confident):
python /path/to/TIPTOP/IP_PDB2txyz.py 1EHZ_fixed.pdb ABC
```

Intermediate files produced:

| File | Description |
|------|-------------|
| `*.biomol.pdb` | Biomolecule-only PDB (no water/ions) |
| `*.solvent.pdb` | Solvent/ion-only PDB |
| `*.biomol.txyz` | Raw TXYZ from obabel (no AMOEBA types yet) |
| `solvent.txyz` | Solvent/ion TXYZ with AMOEBA types |
| `pdblist` | List of per-residue fragment PDB filenames |
| `txyzlist` | List of per-residue fragment TXYZ filenames |
| `NNNN_<RES>.pdb` | Individual residue fragment PDB files |
| `NNNN_<RES>.txyz_2` | Individual residue fragment typed TXYZ files |
| `final.txyz` | **Final output** — fully typed TINKER XYZ |

---

## Non-standard Residues

If your PDB contains a residue not in the database (e.g. a small-molecule
ligand), generate a TXYZ template with `IP_ParmGen.py` or another tool, then
supply it at run time:

```bash
# Single non-standard residue
python IP_PDB2txyz.py complex.pdb ABC --residue LIG.txyz

# Multiple non-standard residues
python IP_PDB2txyz.py complex.pdb ABC --residue LIG.txyz COF.txyz

# Multiple template fragments for the same residue
python IP_PDB2txyz.py complex.pdb ABC --residue LIG_1.txyz LIG_2.txyz
```

The residue name matched in the PDB is the **stem** of the filename
(e.g. `LIG_1.txyz` → residue `LIG`).

If Step B reports missing residues, the output lists the affected fragments
and suggests how to resolve them:

```
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
WARNING: Some residue TXYZ files were NOT generated.
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  Options to resolve missing residues:
    1. Re-run with --residue <RESNAME>.txyz for each missing residue
    2. Add <RESNAME>.txyz to the database: database.PDB2txyz/<category>/
  Then re-run Step B.
```

---

## Extending the Database

To permanently support a new residue, place its TXYZ template in the
appropriate subdirectory of `database.PDB2txyz/`:

```
database.PDB2txyz/
├── aminoacid/      ← standard amino acids
├── nucleotide/     ← DNA / RNA
├── glycan/         ← glycan residues
├── aaconjugate/    ← amino-acid conjugates
├── lipid/          ← lipid molecules
├── solvent/        ← water (HOH) and variants
└── ion/            ← monatomic ions (NA, CLA, K, MG …)
```

Naming conventions for template files:

- Standard residue: `ALA.txyz`
- N-/C-terminal variant: `ALAN.txyz`, `ALAC.txyz`
- Multi-fragment residue: `DA_1.txyz`, `DA_2.txyz`, `DA_3.txyz`

---

## CHARMM-GUI PDB Compatibility

The script automatically converts CHARMM-GUI nucleotide residue names to
the internal naming convention used by the database:

| CHARMM-GUI name | Internal name |
|-----------------|---------------|
| `CYT` | `DC` |
| `GUA` | `DG` |
| `ADE` | `DA` |
| `THY` | `DT` |

---

## Solvent / Ion Aliases

The script recognises common PDB variant names for water and ions:

| PDB name(s) | Canonical template |
|-------------|-------------------|
| `WAT`, `TIP3`, `TIP4`, `SPC`, `SOL` | `HOH` |
| `SOD` | `NA` |
| `CL` | `CLA` |
| `POT` | `K` |
| `MG2` | `MG` |

---

## License

MIT License — Copyright (c) 2021 Chengwen Liu
