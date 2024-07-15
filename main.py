import math
import os
import sys
import time
import urllib.request

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Structure import Structure
from Bio.SVDSuperimposer import SVDSuperimposer
from matplotlib import pyplot as plt  # use to plot the scipy dendrogram
from scipy.cluster.hierarchy import dendrogram, linkage  # use for scipy clustering + visualization
from scipy.spatial.distance import squareform  # use to convert redundant distance matrix to condensed distance matrix

####################################
# Import structures file as an array
####################################

# Import the pdb entry names where each row is an individual chain, column [0] is the entry name, and column [1] is
# the chain identifier. Entries that do not have a chain identifier only have one chain. LIST HAS BEEN FILTERED TO
# REMOVE ENTRIES WITH NON-IDENTICAL SEQUENCES. GOAL WOULD BE TO INCLUDE THOSE SOMEHOW, BUT CURRENT CODE CANNOT DO THAT

assert os.path.exists("./pdbEntries.csv"), ("The pdbEntries.csv file is missing from the root directory. Please add a "
                                            ".csv file that includes a list of PDB codes in column A and (optional) "
                                            "chains in column B that correspond to the PDB codes in column A. Each row "
                                            "will be read as an entry.")

pdbEntries = np.genfromtxt("pdbEntries.csv", dtype=str, encoding="utf-8-sig", delimiter=",", usemask=True)
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

# Download the pdbEntries into the pdbFileDir
# NOTE: if you want to include a local .pdb file, simply move it to the ./pdbFiles folder AND include it in your
# pdbEntires.csv file with the IDENTICAL name as the .pdb file.

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
    e = resi_end  # This is not - 1 because slicing does not include the last referenced number.
    template_coords_specific = template_coords[s:e]
    mobile_coords_rotated = np.dot(mobile_coords, rot) + tran
    mobile_coords_rotated_specific = mobile_coords_rotated[s:e]
    diff = template_coords_specific - mobile_coords_rotated_specific
    rmsd = np.sqrt(sum(sum(diff ** 2)) / template_coords_specific.shape[0])
    return rmsd


# alignment function


def align(template: Structure,
          mobile: Structure,
          align_resi_start=None,
          align_resi_end=None,
          rmsd_resi_start=None,
          rmsd_resi_end=None,
          atom_type="CA") -> list:
    """Aligns a mobile structure onto a template structure using the atom types listed in 'atom_types' and calculates
    the RMSD for a specific region of that alignment."""

    # Get the coordinates of the Atom object if the Atom is from an amino acid residue,
    # and the atom type is what's specified in atom_types. Known modified / non-canonical amino acids are supported.

    template_coords = np.array([a.get_coord() for a in template.get_atoms() if
                                is_aa(a.parent.get_resname()) and a.get_id() == atom_type])
    mobile_coords = np.array([a.get_coord() for a in mobile.get_atoms() if
                              is_aa(a.parent.get_resname()) and a.get_id() == atom_type])

    # Determine if alignment is global or a small region of the enzyme. If you have sequences that include things at N
    # or C term, there might be coordinates in there you don't want to consider.

    if align_resi_start is None:
        a = None
    else:
        a = align_resi_start - 1
    if align_resi_end is None:
        b = None
    else:
        b = align_resi_end  # This is not - 1 because slicing does not include the last referenced number.

    # Set up alignment

    si = SVDSuperimposer()
    si.set(template_coords[a:b], mobile_coords[a:b])
    si.run()  # Run the SVD alignment

    # Determine if going to calculate the rmsd for a specific region of the protein.

    if not ((rmsd_resi_start is None) and (rmsd_resi_end is None)):
        rmsd_spec = rmsd_specific(template_coords, mobile_coords, rmsd_resi_start, rmsd_resi_end, si.rot, si.tran)
    else:
        rmsd_spec = None

    return [si, rmsd_spec]


#################
# Distance matrix
#################

# Set up distance matrix arrays

distance_matrix_global = np.empty((numEntries, numEntries))
distance_matrix_specific = np.empty((numEntries, numEntries))

# Make the distance matrix by iterating through each row / column, performing the alignment, and recording the
# respective rmsd. If any distance matrix exists in the ./distance_matrices directory, this will be skipped.

if not (os.path.exists("./distance_matrices/distance_matrix_global.csv")
        or os.path.exists("./distance_matrices/distance_matrix_global.npy")
        or os.path.exists("./distance_matrices/distance_matrix_specific.csv")
        or os.path.exists("./distance_matrices/distance_matrix_specific.npy")):
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
            global_alignment, specificRMSD = align(templateStruc,
                                                   mobileStruc,
                                                   align_resi_start=1,
                                                   align_resi_end=159,
                                                   rmsd_resi_start=9,
                                                   rmsd_resi_end=24)
            distance_matrix_global[row, column] = global_alignment.get_rms()
            distance_matrix_specific[row, column] = specificRMSD
        end = time.time()
        elapsed = end - start
        remainingMin = elapsed * (numEntries - row + 1) // 60
        remainingSec = elapsed * (numEntries - row + 1) % 60
        print(f"Alignment {100 * (row + 1) // numEntries} % complete. Completed row {row + 1} in {elapsed} seconds. "
              f"Estimated time to completion: {math.floor(remainingMin)} minutes {math.floor(remainingSec)} seconds.")

    # print(distance_matrix_global.reshape(np.ma.shape(pdbEntries)[0], np.ma.shape(pdbEntries)[0]))
    # print(distance_matrix_specific.reshape(np.ma.shape(pdbEntries)[0], np.ma.shape(pdbEntries)[0]))

    # Check if the /distance_matrices folder exists in the current directory. If not, create it.

    distance_matrices_FileDir = "./distance_matrices"

    if not os.path.exists(distance_matrices_FileDir):
        os.makedirs(distance_matrices_FileDir)
        print("Created file directory " + distance_matrices_FileDir)

    # Save .csv and .npy files of the distance matrices.

    np.savetxt(fname="./distance_matrices/distance_matrix_global.csv", X=distance_matrix_global, delimiter=",")
    np.save("./distance_matrices/distance_matrix_global.npy", distance_matrix_global)

    np.savetxt(fname="./distance_matrices/distance_matrix_specific.csv", X=distance_matrix_specific, delimiter=",")
    np.save("./distance_matrices/distance_matrix_specific.npy", distance_matrix_specific)

    print(f"Distance matrix generated for the alignment of {numEntries} structures.")

    distance_matrix_ran = True
else:
    print("No distance matrix was calculated because there is already a distance matrix in ./distance_matrices")
    distance_matrix_ran = False

##############################################
# Hierarchical Clustering from Distance Matrix
##############################################

# Load the distance matrices from saved files in the ./distance_matrices folder

if not distance_matrix_ran:
    distance_matrix_global = np.load("./distance_matrices/distance_matrix_global.npy")
    distance_matrix_specific = np.load('./distance_matrices/distance_matrix_specific.npy')

# Must condense the matrix for linkage() to read. checks=False because the matrix is essentially symmetrical and the
# diagonal elements are essentially zero.

distance_matrix_global_condensed = squareform(distance_matrix_global, checks=False)
distance_matrix_specific_condensed = squareform(distance_matrix_specific, checks=False)

# perform the hierarchical clustering using the average-linkage method

globalRMSDTree = linkage(distance_matrix_global_condensed, "average", optimal_ordering=True)
specificRMSDTree = linkage(distance_matrix_specific_condensed, "average", optimal_ordering=True)

# Check if the /dendrograms folder exists in the current directory. If not, create it.

dendrogramsFileDir = "./dendrograms"

if not os.path.exists(dendrogramsFileDir):
    os.makedirs(dendrogramsFileDir)
    print("Created file directory " + dendrogramsFileDir)

# generate the dendrograms using the matplotlib package's pyplot module (overwrites any current files)

globalRMSD_fig = plt.figure(figsize=(6.5, 10), dpi=600)
global_dn = dendrogram(globalRMSDTree, color_threshold=0, orientation="left", labels=structureList, leaf_font_size=3,
                       above_threshold_color="k")
plt.xlabel("RMSD (Å) of Full-Length EcDHFR (Residues 1 to 159)")
plt.savefig(fname="./dendrograms/globalRMSDTree.pdf")

specificRMSD_fig = plt.figure(figsize=(6.5, 10), dpi=600)
specific_dn = dendrogram(specificRMSDTree, color_threshold=0, orientation="left", labels=structureList,
                         leaf_font_size=3, above_threshold_color="k")
plt.xlabel("RMSD (Å) of the Met20 Loop (Residues 9 to 24)")

plt.savefig(fname="./dendrograms/specificRMSDTree.pdf")

print(f"Hierarchical clustering performed and dendrograms generated for the alignment of {numEntries} structures.")
