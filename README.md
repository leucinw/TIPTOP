# TIPTOP
**T**inker **I**nput **P**reparing & **T**o **O**utput **P**rocessing 

This software package is providing a bunch of light tools for Tinker simulations, including generation of force field parameters and tinker input structures. Each Python program can also be used individually to fulfill certain goals. Below is brief summaries of the programs:

	* IP_AtomTyper.py: modify atom types in xyz and key file 
	* IP_MatchTXYZ.py: assign atom types for an interested molecule by looking at the reference molecule
	* IP_MatchTXYZ_Glycan.py: assign atom types for an interested glycan fragment by looking at the reference glycan (with all hydrogens)
	* IP_ParmGen.py: generate the AMOEBA parameters for a small molecule 
	* IP_PDB2txyz.py: generate a tinker xyz file based on a PDB input
	* IP_TorDriver.py: do torsion scan at QM level and fit AMOEBA torsional parameters
	* IP_TorParmGen.py: quickly assign AMOEBA torsional parameters without doing QM torsion scan
