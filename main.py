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

####################################
# Import structures file as an array
####################################

# Import the pdb entry names where each row is an individual chain, column [0] is the entry name, and column [1] is
# the chain identifier. Entries that do not have a chain identifier only have one chain.
pdbEntries = np.genfromtxt("pdbEntries.csv", dtype=str, encoding="utf-8-sig", delimiter=",", usemask=True)

# Make a list (actually, a numpy array) that has the PDB code (+ chain id if there is one) for each structure
structureList = np.empty(np.ma.shape(pdbEntries)[0], dtype=np.dtype('U100'))

for struc in range(np.ma.shape(pdbEntries)[0]):
    if np.ma.getmask(pdbEntries[struc, 1]):
        structureList[struc] = pdbEntries[struc, 0]
    else:
        structureList[struc] = f"{pdbEntries[struc, 0]}_{pdbEntries[struc, 1]}"


##################################
# Retrieve Structures from the PDB
##################################

# I wanted to use the PDBList.retrieve_pdb_file() method but on Nov 1 2024 the FTP protocol that this method uses to
# download PDB files goes offline.
# Instead, I'm going to use this download_pdb function I found on stack overflow
# https://stackoverflow.com/questions/37335759/using-python-to-download-specific-pdb-files-from-protein-data-bank


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

for struc in range(np.ma.shape(pdbEntries)[0]):
    download_pdb(pdbEntries[struc, 0], pdbFileDir)

###########
# Alignment
###########

# set up alignment
p = PDBParser(QUIET=True)
templateStruc = p.get_structure("templateStruc", "pdbFiles/1hso.pdb")
mobileStruc = p.get_structure("mobileStruc", "pdbFiles/1ht0.pdb")

chain_a = templateStruc[0]['A']
residues_a = [r for r in chain_a]


def align(template: Structure, mobile: Structure, atom_types=["CA", "N", "C", "O"]) -> SVDSuperimposer:
    """Aligns a mobile structure onto a template structure using the atom types listed in 'atom_types'."""

    # A long one-liner that gets the one-letter amino acid representation for each residue in a structure,
    # then joins those letters into one long string.
    template_seq = "".join([index_to_one(three_to_index((r.get_resname())))
                            for r in template[0].get_residues() if is_aa(r)])
    mobile_seq = "".join([index_to_one(three_to_index((r.get_resname())))
                          for r in mobile[0].get_residues() if is_aa(r)])

    # Some assertions that can be used
    assert len(mobile_seq) == len(template_seq), "The sequences should be of identical length."

    # Get the coordinates of the Atom object if the Atom is from an amino acid residue,
    # and the atom type is what's specified in atom_types.
    # Traditionally RMSD is calculated for either:
    # Only the alpha-carbon atoms (CA), or
    # The "protein backbone" atoms (CA, N, C, O), or
    # All atoms
    template_coords = [a.get_coord() for a in template[0].get_atoms() if
                       is_aa(a.parent.get_resname()) and a.get_id() in atom_types]
    mobile_coords = [a.get_coord() for a in mobile[0].get_atoms() if
                     is_aa(a.parent.get_resname()) and a.get_id() in atom_types]

    si = SVDSuperimposer()
    si.set(np.array(template_coords), np.array(mobile_coords))
    si.run()  # Run the SVD alignment

    return si


align_structures = align(templateStruc, mobileStruc)
print("RMSD before alignment: {:.2f} angstroms; full-backbone RMSD after alignment: {:.2f} angstroms".format(
    align_structures.get_init_rms(),
    align_structures.get_rms()))

#################
# Distance matrix
#################

# Make distance_matrix 11 x 11 array, fill with zeros
distance_matrix = np.zeros((11, 11))
# distance_matrix = np.empty((np.ma.shape(pdbEntries)[0], np.ma.shape(pdbEntries)[0]))


# Make 1D array with each axis (remember that in range() and np.arange() the stop parameter is NOT included)
x_axis = np.arange(0, 11, 1)
y_axis = x_axis


# Make a function to call in the loop
def multiply(mult1, mult2):
    return mult1 * mult2


# Change each element of the array to be the product of the corresponding entry in the two axes
for row in y_axis:
    for column in x_axis:
        distance_matrix[row, column] = multiply(y_axis[row], x_axis[column])

print(distance_matrix.reshape(11, 11))
