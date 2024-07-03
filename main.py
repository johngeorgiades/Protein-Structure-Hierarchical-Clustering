# importing modules
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, index_to_one, three_to_index
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

# Let's get our two structures; call them native (true) and model
#struc_list = ["1hso", "1ht0"]
# PDBList().retrieve_pdb_file(pdb_code=struc_list[0], file_format="pdb")
# PDBList().retrieve_pdb_file(pdb_code=struc_list[1], file_format="pdb")

#PDBList().download_pdb_files(struc_list,file_format="pdb")
# for struc in struc_list:
#     print(struc)
#     PDBList().retrieve_pdb_file(pdb_code=struc, file_format="pdb")
p = PDBParser(QUIET=True)
native = p.get_structure("native", "1hso.pdb")
model = p.get_structure("model", "1ht0.pdb")

chain_a = native[0]['A']
residues_a = [r for r in chain_a]

# AA = ["ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN",
#       "ARG", "SER", "THR", "VAL", "TRP", "TYR"]

def align(native: Structure, model: Structure, atom_types=["CA", "N", "C", "O"]) -> SVDSuperimposer:
    """Aligns a model structure onto a native structure using the atom types listed in 'atom_types'."""

    # A long one-liner that gets the one-letter amino acid representation for each residue in a structure,
    # then joins those letters into one long string.
    native_seq = "".join([index_to_one(three_to_index((r.get_resname()))) for r in native[0].get_residues() if is_aa(r)])

# for r in native[0].get_residues():
#     if is_aa(r):
#         index_to_one, three_to_index
#         print(index_to_one(three_to_index((r.get_resname()))))

# native_seq = "".join([index_to_one(three_to_index((r.get_resname()))) for r in native[0].get_residues() if is_aa(r)])
#
# # native_seq = "".join([index_to_one(r) for r in native[0].get_residues() if is_aa(r)])
# print(native_seq)
