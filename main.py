# importing modules
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Structure import Structure
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

# Let's get our two structures; call them native (true) and model
p = PDBParser(QUIET=True)
native = p.get_structure("native", "1hso.pdb")
model = p.get_structure("model", "1ht0.pdb")

chain_a = native[0]['A']
residues_a = [r for r in chain_a]

AA = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
      "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

def align(native: Structure, model:Structure, atom_types = ["CA", "N", "C", "O"]) -> SVDSuperimposer:
      """Alignts a model structure onto a native structure using the atom types listed in 'atom_types'."""

      # A long one-liner that gets the one-letter amino acid representation for each residue in a structure,
      # then joins those letters into one long string.
      native_seq = "".join([three_to_one(r.resname) for r in native[0].get_residues() if r.resname in AA])