# importing modules
from Bio.PDB.PDBParser import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

# Let's get our two structures; call them native (true) and model
p = PDBParser(QUIET=True)
native = p.get_structure("native", "1hso.pdb")
model = p.get_structure("model", "1ht0.pdb")

chain_a = native[0]['A']
residues_a = [r for r in chain_a]

