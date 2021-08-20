# TIPTOP

**TIPTOP: From Tinker Input Preparing To Output Processing**

This software package is providing a bunch of light tools for Tinker simulations, including generation of force field parameters and tinker input structures. Each Python program can also be used individually to fulfill certain goals. Below is brief summaries of the programs:

* IP_AtomTyper.py: modify atom types in tinker xyz and key files 
* IP_MatchTXYZ.py: assign atom types for an interested molecule by refering a reference molecule
* IP_MatchTXYZ_Glycan.py: assign atom types for an interested glycan fragment by refering the reference glycan (with all hydrogens)
* IP_ParmGen.py: generate the AMOEBA parameters for a small molecule (using Gaussian and Tinker)
* IP_PDB2txyz.py: generate a tinker xyz file based on a standard PDB structure
* IP_TorDriver.py: do classical torsion scan at given QM level and fit AMOEBA torsional parameters
* IP_TorParmGen.py: quickly assign AMOEBA torsional parameters without doing QM torsion scan
