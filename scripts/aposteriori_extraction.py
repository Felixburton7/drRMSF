# OLD CODE THAT WORKS

# #!/usr/bin/env python3
# """
# aposteriori_extraction.py

# Extracts necessary structural information from mdCATH domain .h5 files and
# outputs cleaned PDB files for compatibility with `aposteriori`. 
# Each output file corresponds to a single domain and is processed using pdbUtils.

# Usage:
# python aposteriori_extraction.py \
# --input_dir /path/to/mdCATH \
# --output_dir /path/to/output_pdbs \
# --format pdb \
# [--verbose]
# [--clean]

# Author: Felix
# """

# import os
# import re
# import glob
# import h5py
# import argparse
# import numpy as np
# from tqdm import tqdm
# from pdbUtils import pdbUtils

# def parse_domain_id(h5_filename):
#     """
#     Extracts the domain ID from an mdCATH .h5 filename.
#     Example: 'mdcath_dataset_1r9lA02.h5' -> '1r9lA02'.
#     """
#     basename = os.path.basename(h5_filename) #basename = takes in path as input, extracts the filename from path provided, takes away everything except the the filename. 
#     # Assumes files are named like 'mdcath_dataset_<domain_id>.h5'.
#     match = re.match(r"mdcath_dataset_(.+)\.h5", basename) # re (regular expression) - extracts the domain id. '()' capturing group, '.' matches any character except new line, '+) One or more of any character in preceding group. 
#     if not match:
#         domain_id = os.path.splitext(basename)[0] # returns a TUPLE datatype which you can then take the first '0' from. 
#         if domain_id.startswith("mdcath_dataset_"): 
#             domain_id =  domain_id.replace("mdcath_dataset_", "")
#         return domain_id
#         # Fall back to removing extension if pattern not matched
#         # return os.path.splitext(basename)[0] # This will return 'mdcath_dataset_1r91A02' which isn't ideal
#     return match.group(1) # If match successful, it extracts the 'capture group' AKA domain_id which is part of the string that matched the (.+) pattern. 

# def extract_pdb(domain_file, output_dir, output_format="pdb", verbose=False, clean=False):
#     """
#     Extracts the protein structure (as PDB) from the specified mdCATH
#     domain HDF5 file and saves it to the output directory.

#     Parameters:
#     -----------
#     domain_file : str
#         Path to the .h5 file containing the domain data.
#     output_dir : str
#         Directory to which the extracted PDB file will be written.
#     output_format : str
#         Output format for the structural file (default='pdb').
#     verbose : bool
#         If True, prints additional processing information.
#     clean : bool
#         If True, uses pdbUtils to clean the PDB file.
#     """
#     domain_id = parse_domain_id(domain_file) # This line calls the `parse_domain_id` function and assigns the result (the extracted domain ID) to the variable `domain_id`.
#     # The `parse_domain_id` function takes the `domain_file` (which is the path to the current .h5 file being processed) as input.
#     # This function is responsible for extracting the domain ID from the filename using regular expressions or a fallback mechanism.


#     try:
#         with h5py.File(domain_file, "r") as h5f: # open in read mode. 
#             # The top-level group is typically the domain ID
#             # but we shall retrieve it dynamically in case the
#             # file layout changes.
#             top_keys = list(h5f.keys()) 
#             if len(top_keys) == 0: 
#                 if verbose:
#                     print(f"[WARNING] No groups found in {domain_file}. Skipping.")
#                 return

#             # We handle each domain group in the file (commonly just one).
#             for dkey in top_keys:
#                 # Check if this is truly a domain group (mdCATH usually has the same name).
#                 if verbose:
#                     print(f"Processing domain group '{dkey}' in file '{domain_file}'.")
                
#                 # Retrieve the PDB from 'pdbProteinAtoms' or fall back to 'pdb'.
#                 if "pdbProteinAtoms" in h5f[dkey]:
#                     pdb_string = h5f[dkey]["pdbProteinAtoms"][()].decode("utf-8")
#                 elif "pdb" in h5f[dkey]:
#                     pdb_string = h5f[dkey]["pdb"][()].decode("utf-8")
#                 else:
#                     if verbose:
#                         print(f"[WARNING] No PDB field found in {domain_file}. Skipping domain '{dkey}'.")
#                     continue

#                 # Write the initial file
#                 out_filename = f"{domain_id}.{output_format}" #Gives a domain_id.pdb
#                 out_path = os.path.join(output_dir, out_filename)
#                 with open(out_path, "w") as f_out:
#                     f_out.write(pdb_string)

#                 # Clean the PDB if requested
#                 if clean:
#                     cleaned_filename = f"{domain_id}_cleaned.{output_format}"
#                     cleaned_path = os.path.join(output_dir, cleaned_filename)
                    
#                     # Convert to DataFrame and back to cleaned PDB
#                     pdb_df = pdbUtils.pdb2df(out_path)
#                     pdbUtils.df2pdb(pdb_df, cleaned_path)

#                     if verbose:
#                         print(f"Saved cleaned PDB to {cleaned_path}")
#                 else:
#                     if verbose:
#                         print(f"Saved {out_path}")

#     except Exception as e:
#         if verbose:
#             print(f"[ERROR] Unable to process '{domain_file}': {e}")

# def aposteriori_extraction(input_path, output_dir, output_format="pdb", verbose=False, clean=False):
#     """
#     Main extraction function. Handles both single .h5 file and directory input.

#     Parameters:
#     -----------
#     input_path : str
#         Path to a single .h5 file or a directory containing multiple .h5 files.
#     output_dir : str
#         Directory to which extracted files will be saved.
#     output_format : str
#         Desired output format (e.g., 'pdb').
#     verbose : bool
#         If True, prints additional processing information.
#     clean : bool
#         If True, cleans PDB files using pdbUtils.
#     """
#     # Ensure output directory exists.
#     os.makedirs(output_dir, exist_ok=True)

#     # Determine if input is a directory or a single file.
#     if os.path.isdir(input_path):
#         # Recursively gather .h5 files.
#         domain_files = glob.glob(os.path.join(input_path, "**/*.h5"), recursive=True) #Join to construct a path pattern, 
#     else:
#         domain_files = [input_path]

#     if verbose:
#         print(f"Found {len(domain_files)} .h5 file(s). Beginning extraction.")

#     # File-level progress bar.
#     for domain_file in tqdm(domain_files, desc="Extracting PDB", unit="file"):
#         extract_pdb(domain_file, output_dir, output_format=output_format, verbose=verbose, clean=clean)

# def main():
#     #Argspase library used for parsing CLI arguments.
#     parser = argparse.ArgumentParser(description="Extract structural data from mdCATH domain .h5 files for aposteriori.") #Create parser object which can be used to parse CLI argument.
#     parser.add_argument("--input_dir", required=True, help="Path to a single .h5 file or a directory of .h5 files.")
#     parser.add_argument("--output_dir", required=True, help="Directory where extracted PDB files will be saved.")
#     parser.add_argument("--format", default="pdb", help="Output file format. Typically 'pdb'.")
#     parser.add_argument("--verbose", action="store_true", help="Increase output verbosity.") 
#     parser.add_argument("--clean", action="store_true", help="Clean PDB files using pdbUtils.")

#     args = parser.parse_args() #Call parse_args method on parser to store in an object called args. 

#     aposteriori_extraction(
#         input_path=args.input_dir, #Explicitly pass in the values - The input_dir argument is the input_path in aposteriori script. 
#         output_dir=args.output_dir,
#         output_format=args.format,
#         verbose=args.verbose,
#         clean=args.clean
#     )

# if __name__ == "__main__":
#     main()

# '''

# Task: 
# This currently extract the PDB file and cleans it then give it an output. 
# I want to change this function. It should be able to look at the .H5 file, look at a particular at the temperatures, look at the the repeats nested in that temperature. 
# Look at the at the coordinates in this repeat. Match these coordinates to the specific residue/name number etc in the extracted PDB file. 
# So if I change 'n' to '5' that means number of frames '5' or if I choose 'n' then I want '10 frames from that specific temperature. It can randomly select them from different repeats within that temperature. I want a pdb file identical to the currently extracted (post pdbUtils one) BUT it has different coordinates. The coordinates related to that frame. 
# So it creates 5 files which are identical but with different coordinates in this column related column '-3.090  -3.392 -18.053' 
# ATOM      1 CAY  ALA 0  92      -3.090  -3.392 -18.053  0.00  0.00           C  
# ATOM      2 HY1  ALA 0  92      -3.921  -3.067 -18.503  0.00  0.00           H  
# ATOM      3 HY2  ALA 0  92      -3.233  -4.316 -17.698  0.00  0.00           H  
# ATOM      4 HY3  ALA 0  92      -2.328  -3.389 -18.700  0.00  0.00           H  
# ATOM      5 CY   ALA 0  92      -2.873  -2.781 -17.291  0.00  0.00           C  
# ATOM      6 OY   ALA 0  92      -2.041  -3.106 -16.840  0.00  0.00           O  
# ATOM      7 N    ALA 0  92      -3.563  -2.059 -17.246  1.00  0.00           N  
# ATOM      8 HN   ALA 0  92      -4.395  -1.734 -17.697  0.00  0.00           H  
# ATOM      9 CA   ALA 0  92      -3.248  -1.173 -16.141  1.00  0.00           C  
# ATOM     10 HA   ALA 0  92      -2.335  -0.801 -16.267  1.00  0.00           H  
# ATOM     11 CB   ALA 0  92      -4.235  -0.026 -16.068  1.00  0.00           C 

# '''