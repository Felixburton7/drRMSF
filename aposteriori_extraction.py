#Pip install pdbUtils
#From pdbUtils import pdbUtils
#Dataframe = pdbUtils.pdb2df
#df2pdb(path/)

#If your pdb is a dataframe and 
#Just to get labels. 


#!/usr/bin/env python3
"""
aposteriori_extraction.py

Extracts necessary structural information from mdCATH domain .h5 files and
outputs PDB files for compatibility with `aposteriori. Each output file
corresponds to a single domain.

Usage:
  python aposteriori_extraction.py \
      --input_dir /path/to/mdCATH \
      --output_dir /path/to/output_pdbs \
      --format pdb \
      [--verbose]

Author: Researcher / Senior Scientist
"""

import os
import re
import glob
import h5py
import argparse
import numpy as np
from tqdm import tqdm

def parse_domain_id(h5_filename):
    """
    Extracts the domain ID from an mdCATH .h5 filename.
    Example: 'mdcath_dataset_1r9lA02.h5' -> '1r9lA02'.
    """
    basename = os.path.basename(h5_filename)
    # Assumes files are named like 'mdcath_dataset_<domain_id>.h5'.
    match = re.match(r"mdcath_dataset_(.+)\.h5", basename)
    if not match:
        # Fall back to removing extension if pattern not matched
        return os.path.splitext(basename)[0]
    return match.group(1)

def extract_pdb(domain_file, output_dir, output_format="pdb", verbose=False):
    """
    Extracts the protein structure (as PDB) from the specified mdCATH
    domain HDF5 file and saves it to the output directory.

    Parameters:
    -----------
    domain_file : str
        Path to the .h5 file containing the domain data.
    output_dir : str
        Directory to which the extracted PDB file will be written.
    output_format : str
        Output format for the structural file (default='pdb').
    verbose : bool
        If True, prints additional processing information.
    """
    domain_id = parse_domain_id(domain_file)

    try:
        with h5py.File(domain_file, "r") as h5f:
            # The top-level group is typically the domain ID
            # but we shall retrieve it dynamically in case the
            # file layout changes.
            top_keys = list(h5f.keys())
            if len(top_keys) == 0:
                if verbose:
                    print(f"[WARNING] No groups found in {domain_file}. Skipping.")
                return

            # We handle each domain group in the file (commonly just one).
            for dkey in top_keys:
                # Check if this is truly a domain group (mdCATH usually has the same name).
                if verbose:
                    print(f"Processing domain group '{dkey}' in file '{domain_file}'.")
                # Retrieve the PDB from 'pdbProteinAtoms' or fall back to 'pdb'.
                if "pdbProteinAtoms" in h5f[dkey]:
                    pdb_string = h5f[dkey]["pdbProteinAtoms"][()].decode("utf-8")
                elif "pdb" in h5f[dkey]:
                    pdb_string = h5f[dkey]["pdb"][()].decode("utf-8")
                else:
                    if verbose:
                        print(f"[WARNING] No PDB field found in {domain_file}. Skipping domain '{dkey}'.")
                    continue

                # Write the file with the appropriate extension.
                out_filename = f"{domain_id}.{output_format}"
                out_path = os.path.join(output_dir, out_filename)
                with open(out_path, "w") as f_out:
                    f_out.write(pdb_string)
                if verbose:
                    print(f"Saved {out_path}")
    except Exception as e:
        if verbose:
            print(f"[ERROR] Unable to process '{domain_file}': {e}")

def aposteriori_extraction(input_path, output_dir, output_format="pdb", verbose=False):
    """
    Main extraction function. Handles both single .h5 file and directory input.

    Parameters:
    -----------
    input_path : str
        Path to a single .h5 file or a directory containing multiple .h5 files.
    output_dir : str
        Directory to which extracted files will be saved.
    output_format : str
        Desired output format (e.g., 'pdb').
    verbose : bool
        If True, prints additional processing information.
    """
    # Ensure output directory exists.
    os.makedirs(output_dir, exist_ok=True)

    # Determine if input is a directory or a single file.
    if os.path.isdir(input_path):
        # Recursively gather .h5 files.
        domain_files = glob.glob(os.path.join(input_path, "**/*.h5"), recursive=True)
    else:
        domain_files = [input_path]

    if verbose:
        print(f"Found {len(domain_files)} .h5 file(s). Beginning extraction.")

    # File-level progress bar.
    for domain_file in tqdm(domain_files, desc="Extracting PDB", unit="file"):
        extract_pdb(domain_file, output_dir, output_format=output_format, verbose=verbose)


def main():
    parser = argparse.ArgumentParser(description="Extract structural data from mdCATH domain .h5 files for aposteriori.")
    parser.add_argument("--input_dir", required=True, help="Path to a single .h5 file or a directory of .h5 files.")
    parser.add_argument("--output_dir", required=True, help="Directory where extracted PDB files will be saved.")
    parser.add_argument("--format", default="pdb", help="Output file format. Typically 'pdb'.")
    parser.add_argument("--verbose", action="store_true", help="Increase output verbosity.")
    args = parser.parse_args()

    aposteriori_extraction(
        input_path=args.input_dir,
        output_dir=args.output_dir,
        output_format=args.format,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()


#Somewhere they will have the raw trajectory files - we can use a program like mdAnalysis or pdbframes. 
#Pretend you want to find the trajectory file - dcd 


#Load the the coordinates and trajectory x3 x3 
#Every period just replace the coordinate columns
#Have pdb open as a dataframe and then update coordinates based on the array. 
#Coordinates will be indexed to each atom in the file. 
# Sampole from trajectory at 5 different points - each will have the same RMSF labels
# BUT
# These will produce very subtly differnet aposteriori frames
# - So this model will not hyper fixate one particular molecule
# - This is protein A in multiple different ways
# - So 5 snapshots, 
# - Crystal structure is at a high minima 
# - Sequence prediction. 

# - From a really big model - start within one - then model two - then model three - 

#Make a funciton take an (n) paramater so we can edit the dataframes 
#Every proper run as configs 
#Anything that ever needs to change 
#Anything change dictionary 
# 
