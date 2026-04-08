from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(pdbid='1EHZ')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms() # This adds the hydrogens based on pH 7 by default
fixer.addMissingHydrogens(7.0)

with open('1EHZ_fixed.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)