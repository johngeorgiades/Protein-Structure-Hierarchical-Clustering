# import necessary modules
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, index_to_one, three_to_index
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import os
import urllib.request
import sys

##################################
# retrieve structures from the PDB
##################################
# I wanted to use the PDBList.retrieve_pdb_file() method but on Nov 1 2024 the FTP protocol that this method uses to
# download PDB files goes offline.
# Instead, I'm going to use this download_pdb function I found on stack overflow
# https://stackoverflow.com/questions/37335759/using-python-to-download-specific-pdb-files-from-protein-data-bank


# download_pdb() function
def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)

    if os.path.exists(outfnm):
        print(f"The file {pdbfn} was not downloaded because it already exists in path {outfnm}")
    else:
        try:
            urllib.request.urlretrieve(url, outfnm)
            return outfnm
        except Exception as err:
            print(str(err), file=sys.stderr)
            return None


# Check if the /pdbFiles folder exists in the current directory. If not, create it.
pdbFileDir = "./pdbFiles"

if not os.path.exists(pdbFileDir):
    os.makedirs(pdbFileDir)
    print("Created file directory " + pdbFileDir)

# import pdb entries to /pdbFiles folder
pdbEntries = ["1hso", '1ht0']

for entry in pdbEntries:
    download_pdb(entry, pdbFileDir)

# set up alignment
p = PDBParser(QUIET=True)
native = p.get_structure("native", "pdbFiles/1hso.pdb")
model = p.get_structure("model", "pdbFiles/1ht0.pdb")

chain_a = native[0]['A']
residues_a = [r for r in chain_a]


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
