# TIPTOP: From Tinker Input Preparing To Output Processing
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5227347.svg)](https://doi.org/10.5281/zenodo.5227347)

## Introduction

This software package provides several light tools for Tinker simulations from input preparing to output processing. Each program can also be used individually for specific purposes. 

* IP_AtomTyper.py: modify atom types in tinker xyz and key files 
* IP_MatchTXYZ.py: assign atom types for an interested molecule by refering a reference molecule
* IP_MatchTXYZ_Glycan.py: assign atom types for an interested glycan fragment by refering the reference glycan (with all hydrogens)
* IP_ParmGen.py: generate the AMOEBA parameters for a small molecule (using Gaussian and Tinker)
* IP_PDB2txyz.py: generate a tinker xyz file based on a standard PDB structure
* IP_TorDriver.py: do classical torsion scan at given QM level and fit AMOEBA torsional parameters
* IP_TorParmGen.py: quickly assign AMOEBA torsional parameters without doing QM torsion scan
* IP_ValenceInit.py: write the necessary bond/angle/strbnd/opbend/torsion parameters basd on a tinker xyz file 
* OP_ARC2PDB.py: convert Tinker trajectory file (arc) to a PDB file
* OP_PDBeditor.py: Edit a PDB file, supported function: trim, average

## Cite TIPTOP
The author hopes TIPTOP can make your Tinker simulations a little easier. It is appreciated if you cite it in your publication.
```
Chengwen Liu. (2021). leucinw/TIPTOP: TIPTOP: From Tinker Input Preparing To Output Processing (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5227347
```

## Tutorial
Run the individual program with `-h` option to get help information. Detailed tutorials for running each program is coming (not) soon.

## Note
Some of the functionality is not perfect yet.
