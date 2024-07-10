import os
import sys
import time
import urllib.request
import math

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Structure import Structure
from Bio.SVDSuperimposer import SVDSuperimposer
from scipy.cluster.hierarchy import dendrogram, linkage  # use for scipy clustering + visualization
from scipy.spatial.distance import squareform  # use to convert redundant distance matrix to condensed distance matrix
from matplotlib import pyplot as plt  # use to plot the scipy dendrogram

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


# RMSD calculation for a specific region


def rmsd_specific(template_coords, mobile_coords, resi_start: int, resi_end: int, rot, tran):
    """Calculates the RMSD between two protein atom coordinate sets 'template_coords' and 'mobile_coords' for a
    specific region of the amino acid sequence between 'resi_start' and 'resi_end'. """
    s = resi_start - 1
    e = resi_end
    template_coords_specific = template_coords[s:e]
    mobile_coords_rotated = np.dot(mobile_coords, rot) + tran
    mobile_coords_rotated_specific = mobile_coords_rotated[s:e]
    diff = template_coords_specific - mobile_coords_rotated_specific
    rmsd = np.sqrt(sum(sum(diff ** 2)) / template_coords_specific.shape[0])
    return rmsd


# alignment function

def align(template: Structure, mobile: Structure, resi_start=None, resi_end=None, atom_type="CA") -> list:
    """Aligns a mobile structure onto a template structure using the atom types listed in 'atom_types' and calculates
    the RMSD for a specific region of that alignment."""

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

    template_coords = np.array([a.get_coord() for a in template.get_atoms() if
                                is_aa(a.parent.get_resname()) and a.get_id() == atom_type])
    mobile_coords = np.array([a.get_coord() for a in mobile.get_atoms() if
                              is_aa(a.parent.get_resname()) and a.get_id() == atom_type])

    si = SVDSuperimposer()
    si.set(template_coords, mobile_coords)
    si.run()  # Run the SVD alignment

    if not ((resi_start is None) and (resi_end is None)):
        rmsd_spec = rmsd_specific(template_coords, mobile_coords, resi_start, resi_end, si.rot, si.tran)
    else:
        rmsd_spec = None

    return [si, rmsd_spec]


#################
# Distance matrix
#################

distance_matrix_global = np.empty((numEntries, numEntries))
distance_matrix_specific = np.empty((numEntries, numEntries))

for row in range(numEntries):
    start = time.time()
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
        global_alignment, specificRMSD = align(templateStruc, mobileStruc, 9, 24)
        distance_matrix_global[row, column] = global_alignment.get_rms()
        distance_matrix_specific[row, column] = specificRMSD
    end = time.time()
    elapsed = end - start
    remainingMin = elapsed * (numEntries - row + 1) // 60
    remainingSec = elapsed * (numEntries - row + 1) % 60
    print(f"Alignment {100 * (row + 1) // numEntries} % complete. Completed row {row + 1} in {elapsed} seconds. "
          f"Estimated time to completion: {math.floor(remainingMin)} minutes {math.floor(remainingSec)} seconds.")

print(distance_matrix_global.reshape(np.ma.shape(pdbEntries)[0], np.ma.shape(pdbEntries)[0]))
np.savetxt(fname="distance_matrix_global.csv", X=distance_matrix_global, delimiter=",")
print(distance_matrix_specific.reshape(np.ma.shape(pdbEntries)[0], np.ma.shape(pdbEntries)[0]))
np.savetxt(fname="distance_matrix_specific.csv", X=distance_matrix_specific, delimiter=",")

##############################################
# Hierarchical Clustering from Distance Matrix
##############################################

# Must condense the matrix for linkage() to read. checks=False because the matrix is essentially symmetrical and the
# diagonal elements are essentially zero.

distance_matrix_global_condensed = squareform(distance_matrix_global, checks=False)
distance_matrix_specific_condensed = squareform(distance_matrix_specific, checks=False)

# perform the hierarchical clustering using the average-linkage method

globalRMSDTree = linkage(distance_matrix_global_condensed, "average", optimal_ordering=True)
specificRMSDTree = linkage(distance_matrix_specific_condensed, "average", optimal_ordering=True)

# generate the dendrogram using the matplotlib package's pyplot module

globalRMSD_fig = plt.figure(figsize=(6.5, 10), dpi=600)
global_dn = dendrogram(globalRMSDTree, orientation="left", labels=structureList)
plt.savefig(fname="globalRMSDTree.pdf")

specificRMSD_fig = plt.figure(figsize=(6.5, 10), dpi=600)
specific_dn = dendrogram(specificRMSDTree, orientation="left", labels=structureList)
plt.savefig(fname="specificRMSDTree.pdf")
