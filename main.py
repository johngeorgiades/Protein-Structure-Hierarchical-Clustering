# importing modules
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, index_to_one, three_to_index
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

# Let's get our two structures; call them native (true) and model
# struc_list = ["1hso", "1ht0"]
# PDBList().retrieve_pdb_file(pdb_code=struc_list[0], file_format="pdb")
# PDBList().retrieve_pdb_file(pdb_code=struc_list[1], file_format="pdb")

# PDBList().download_pdb_files(struc_list,file_format="pdb")
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
    native_seq = "".join(
        [index_to_one(three_to_index((r.get_resname()))) for r in native[0].get_residues() if is_aa(r)])
    model_seq = "".join([index_to_one(three_to_index((r.get_resname()))) for r in model[0].get_residues() if is_aa(r)])

    # Some assertions that can be used
    assert len(model_seq) == len(native_seq), "The sequences should be of identical length."

    # Get the coordinates of the Atom object if the Atom is from an amino acid residue,
    # and the atom type is what's specified in atom_types.
    # Traditionally RMSD is calculated for either:
    # Only the alpha-carbon atoms (CA), or
    # The "protein backbone" atoms (CA, N, C, O), or
    # All atoms
    native_coords = [a.get_coord() for a in native[0].get_atoms() if
                     is_aa(a.parent.get_resname()) and a.get_id() in atom_types]
    model_coords = [a.get_coord() for a in model[0].get_atoms() if
                    is_aa(a.parent.get_resname()) and a.get_id() in atom_types]

    si = SVDSuperimposer()
    si.set(np.array(native_coords), np.array(model_coords))
    si.run()  # Run the SVD alignment

    return si


align_structures = align(native, model)
print("RMSD before alignment: {:.2f} angstroms; full-backbone RMSD after alignment: {:.2f} angstroms".format(
    align_structures.get_init_rms(),
    align_structures.get_rms()))

#################
# Distance matrix
#################

# Make distance_matrix 11 x 11 array, fill with zeros
template = np.zeros((11, 11))

# Make 1D array with each axis (remember that in range() and np.arange() the stop parameter is NOT included)
x_axis = np.arange(0, 11, 1)
y_axis = x_axis


# Make a function to call in the loop
def multiply(mult1, mult2):
    return mult1 * mult2


# Change each element of the array to be the product of the corresponding entry in the two axes
for row in y_axis:
    for column in x_axis:
        template[row, column] = multiply(y_axis[row], x_axis[column])

print(template.reshape(11, 11))