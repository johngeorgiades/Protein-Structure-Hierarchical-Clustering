import os
import sys
import urllib.request

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Structure import Structure
from Bio.SVDSuperimposer import SVDSuperimposer

####################################
# Import structures file as an array
####################################

# Import the pdb entry names where each row is an individual chain, column [0] is the entry name, and column [1] is
# the chain identifier. Entries that do not have a chain identifier only have one chain. LIST HAS BEEN FILTERED TO
# REMOVE ENTRIES WITH NON-IDENTICAL SEQUENCES. GOAL WOULD BE TO INCLUDE THOSE SOMEHOW, BUT CURRENT CODE CANNOT DO THAT

pdbEntries = np.genfromtxt("pdbEntriesFiltered.csv", dtype=str, encoding="utf-8-sig", delimiter=",", usemask=True)
numEntries = np.ma.shape(pdbEntries)[0]

# Make a list (actually, a numpy array) that has the PDB code (+ chain id if there is one) for each structure

structureList = np.empty(numEntries, dtype=np.dtype('U100'))

for struc in range(numEntries):
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

for struc in range(numEntries):
    download_pdb(pdbEntries[struc, 0], pdbFileDir)

###########
# Alignment
###########
p = PDBParser(QUIET=True)


def align(template: Structure, mobile: Structure, atom_types="CA") -> SVDSuperimposer:
    """Aligns a mobile structure onto a template structure using the atom types listed in 'atom_types'."""

    # A long one-liner that gets the one-letter amino acid representation for each residue in a structure,
    # then joins those letters into one long string.

    template_seq = "".join([r.get_resname()
                            for r in template.get_residues() if is_aa(r)])
    mobile_seq = "".join([r.get_resname()
                          for r in mobile.get_residues() if is_aa(r)])

    # Check if identical sequences. Throw error if not (since atom coordinates won't work).

    assert len(mobile_seq) == len(template_seq), "The sequences should be of identical length."

    # Get the coordinates of the Atom object if the Atom is from an amino acid residue,
    # and the atom type is what's specified in atom_types.

    template_coords = [a.get_coord() for a in template.get_atoms() if
                       is_aa(a.parent.get_resname()) and a.get_id() in atom_types]
    mobile_coords = [a.get_coord() for a in mobile.get_atoms() if
                     is_aa(a.parent.get_resname()) and a.get_id() in atom_types]

    si = SVDSuperimposer()
    si.set(np.array(template_coords), np.array(mobile_coords))
    si.run()  # Run the SVD alignment

    return si


#################
# Distance matrix
#################

distance_matrix = np.empty((numEntries, numEntries))

for row in range(numEntries):
    if np.ma.getmask(pdbEntries[row, 1]):
        templateStruc = p.get_structure("templateStruc", f"pdbFiles/{pdbEntries[row, 0]}.pdb")[0]["A"]
    else:
        templateStruc = p.get_structure("templateStruc", f"pdbFiles/{pdbEntries[row, 0]}.pdb")[0][
            pdbEntries[row, 1]]
    for column in range(numEntries):
        if np.ma.getmask(pdbEntries[column, 1]):
            mobileStruc = p.get_structure("mobileStruc", f"pdbFiles/{pdbEntries[column, 0]}.pdb")[0]["A"]
        else:
            mobileStruc = p.get_structure("mobileStruc", f"pdbFiles/{pdbEntries[column, 0]}.pdb")[0][
                pdbEntries[column, 1]]
        align_structures = align(templateStruc, mobileStruc)
        distance_matrix[row, column] = align_structures.get_rms()
    print(f"Alignment {100 * (row + 1) // numEntries} % complete. Starting iteration {row + 2}.")

print(distance_matrix.reshape(np.ma.shape(pdbEntries)[0], np.ma.shape(pdbEntries)[0]))
# np.savetxt(fname=distance_matrix, X=distance_matrix, delimiter=",")
