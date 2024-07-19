# Protein Structure Hierarchical Clustering

## Python Packages

All necessary Python packages are specified in the requirements.txt folder.

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install these using:

```bash
pip install -r requirements.txt
```

## Required Files

### pdbEntries.csv
column A: PDB codes ("XXXX" format) or file names (without extension)

column B: (optional) chain identifier ("A", "B", "C", etc.; defaults to "A")

To include multiple chains from one PDB code, you must specify the PDB code and chain identifier for each
structure as separate rows.

The pdbEntries.csv file must be in your root directory.

## Calling the ```align()``` Function

```python
# Align mobileStruc to templateStruc using residues 1 to 159 and calculating an
# RMSD for residues 9 to 24.
align(templateStruc, mobileStruc, 
      align_resi_start=1, align_resi_end=159, 
      rmsd_resi_start=9,rmsd_resi_end=24)
```

The dimensions of the coordinates matrices for template and mobile structures must be identical.

Specify the residue range (```align_resi_start=``` and ```align_resi_end=```) you seek to use for the alignment and the residue 
range (```rmsd_resi_start=``` and ```rmsd_resi_end=```) you seek to calculate an RMSD for, in addition to the global RMSD. For both 
of these ranges, the "resi_end" residue IS included.

## Explanation of File Directories and Their Contents

### ```"./pdbFiles"```
The ```"./pdbFiles"``` directory will be in your root directory.

The main.py script creates this directory, if it does not already exist, downloads PDB files from the RCSB PDB, and
places them into this directory.

If you wish to include local PDB files for analysis, please place those files in this directory. The file name must
be IDENTICAL to the column A of the corresponding entry in pdbEntries.csv.

### ```"./distance_matrices"```
The ```"./distance_matrices"``` directory will be in your root directory.

The main.py script creates this directory, if it does not already exist, and places in it the calculated distance
matrices in .csv and .npy formats.

If any files exist in this directory, main.py will NOT calculate new distance matrices. This is done to save time,
as distance matrix calculation can take several minutes. Please delete the ```"./distance_matrices"``` directory if you wish
to calculate new distance matrices.

### ```"./dendrograms"```
The ```"./dendrograms"``` directory will be in your root directory.

The main.py script creates this directory, if it does not already exist, and places in it the dendrograms plotted
based on the hierarchical clustering in the .pdb format.

The main.py script WILL overwrite files in this directory.

## Citation

If this code helped you in your work, please consider supporting us by citing the following [paper](https://doi.org/).

## License

[MIT](https://choosealicense.com/licenses/mit/)