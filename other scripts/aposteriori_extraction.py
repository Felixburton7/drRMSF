





# # OLD AND WORKS 
# #!/usr/bin/env python3
# """
# aposteriori_extraction_revised.py

# Revised script to extract and sample coordinate frames from mdCATH domain .h5 files,
# matching them to specific residues in the original PDB file, and outputting updated PDB files.

# Usage:
# python aposteriori_extraction_revised.py \
#     --h5_file mdcath_dataset_1r9lA02.h5 \
#     --output_dir aposteriori_ready \
#     --n 10 \
#     --verbose \
#     [--clean]

# Author: Felix
# """

# import os  # For interacting with the operating system, like file paths
# import re  # For regular expression operations
# import h5py  # To handle HDF5 files
# import argparse  # For parsing command-line arguments
# import numpy as np  # For numerical operations, especially with arrays
# from tqdm import tqdm  # For displaying progress bars
# from pdbUtils import pdbUtils  # Custom utility for handling PDB files
# import random  # For random sampling
# import shutil  # For high-level file operations like copying files

# def parse_domain_id(h5_filename):
#     """
#     Extracts the domain ID from an mdCATH .h5 filename.
#     Example: 'mdcath_dataset_1r9lA02.h5' -> '1r9lA02'.
#     """
#     basename = os.path.basename(h5_filename)  # Gets the filename from the full path
#     # Uses regex to match the pattern and extract the domain ID
#     match = re.match(r"mdcath_dataset_(.+)\.h5", basename)
#     if not match:
#         # If the pattern doesn't match, remove the extension and the prefix manually
#         domain_id = os.path.splitext(basename)[0]
#         if domain_id.startswith("mdcath_dataset_"):
#             domain_id = domain_id.replace("mdcath_dataset_", "")
#         return domain_id
#     return match.group(1)  # Returns the captured domain ID

# def extract_pdb_from_h5(h5_file, domain_id, verbose=False):
#     """
#     Extracts the PDB string from the H5 file using pdbUtils.

#     Parameters:
#     -----------
#     h5_file : h5py.File
#         Opened HDF5 file object.
#     domain_id : str
#         Domain identifier extracted from the filename.
#     verbose : bool
#         If True, prints additional processing information.

#     Returns:
#     --------
#     str
#         The extracted PDB string.
#     """
#     # Access the top-level group corresponding to the domain ID
#     domain_group = h5_file[domain_id]
    
#     # Retrieve the PDB string from 'pdbProteinAtoms' or fallback to 'pdb'
#     if "pdbProteinAtoms" in domain_group:
#         pdb_string = domain_group["pdbProteinAtoms"][()].decode("utf-8")
#         if verbose:
#             print(f"Extracted PDB from 'pdbProteinAtoms'.")
#     elif "pdb" in domain_group:
#         pdb_string = domain_group["pdb"][()].decode("utf-8")
#         if verbose:
#             print(f"Extracted PDB from 'pdb'.")
#     else:
#         raise KeyError(f"No PDB field found in domain '{domain_id}'.")

#     return pdb_string  # Returns the PDB string

# def load_pdb(pdb_path):
#     """
#     Loads the PDB file and returns it as a list of lines.
#     """
#     with open(pdb_path, 'r') as pdb_file:
#         pdb_lines = pdb_file.readlines()  # Reads all lines from the PDB file
#     return pdb_lines  # Returns the list of PDB lines

# def write_pdb(pdb_lines, coords, output_path):
#     """
#     Writes a new PDB file with updated coordinates.

#     Parameters:
#     -----------
#     pdb_lines : list of str
#         Original PDB file lines.
#     coords : np.ndarray
#         Array of shape (numAtoms, 3) containing updated coordinates.
#     output_path : str
#         Path to save the new PDB file.
#     """
#     new_pdb_lines = []  # Initializes a list to hold the new PDB lines
#     coord_idx = 0  # Index to keep track of the current coordinate

#     for line in pdb_lines:
#         if line.startswith("ATOM") or line.startswith("HETATM"):
#             # Ensure we have enough coordinates to update
#             if coord_idx >= coords.shape[0]:
#                 raise ValueError(f"Insufficient coordinates provided. Expected at least {coord_idx + 1}, got {coords.shape[0]}")

#             # Extract the x, y, z coordinates for the current atom
#             x, y, z = coords[coord_idx]
#             # Create a new line with updated coordinates, maintaining the original formatting
#             new_line = (
#                 f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}"
#                 f"{line[54:]}"
#             )
#             new_pdb_lines.append(new_line)  # Add the updated line to the list
#             coord_idx += 1  # Move to the next coordinate
#         else:
#             new_pdb_lines.append(line)  # Non-ATOM/HETATM lines are added unchanged

#     # Check if all coordinates have been used
#     if coord_idx != coords.shape[0]:
#         raise ValueError(f"Unused coordinates detected. Provided {coords.shape[0]}, used {coord_idx}")

#     # Write all the new lines to the output PDB file
#     with open(output_path, 'w') as out_file:
#         out_file.writelines(new_pdb_lines)

# def sample_frames(h5_file, domain_id, temperature, n, verbose=False):
#     """
#     Samples 'n' random frames from any repeat within the specified temperature.

#     Parameters:
#     -----------
#     h5_file : h5py.File
#         Opened HDF5 file object.
#     domain_id : str
#         Domain identifier.
#     temperature : str
#         Temperature to process.
#     n : int
#         Number of frames to sample.
#     verbose : bool
#         If True, prints detailed logs.

#     Returns:
#     --------
#     List of tuples: [(repeat_id, frame_index), ...]
#     """
#     available_replicas = list(h5_file[domain_id][temperature].keys())  # Lists all replicas for the given temperature
#     if verbose:
#         print(f"Available replicas ({temperature}K): {available_replicas}")

#     sampled_frames = []  # List to store sampled (replica, frame) tuples
#     total_frames = 0  # Counter for total available frames

#     # Iterate over each replica to count total frames
#     for repl in available_replicas:
#         num_frames = h5_file[domain_id][temperature][repl]['coords'].shape[0]  # Number of frames in this replica - This accesses the coordinates in teh dataset for a speicifc. (450, 2165, 3) 
#         total_frames += num_frames  # Update total frame count, 
#         if verbose:
#             print(f"Replica {repl} has {num_frames} frames.")

#     # Check if there are enough frames to sample
#     if total_frames < n:
#         raise ValueError(f"Not enough frames ({total_frames}) available to sample {n} frames for temperature {temperature}.")

#     # Randomly sample frames until we have 'n' unique samples
#     while len(sampled_frames) < n:
#         repl = random.choice(available_replicas)  # Randomly select a replica
#         num_frames = h5_file[domain_id][temperature][repl]['coords'].shape[0]  # Frames in the selected replica
#         frame_idx = random.randint(0, num_frames - 1)  # Random frame index within the replica
#         if (repl, frame_idx) not in sampled_frames:
#             sampled_frames.append((repl, frame_idx))  # Add unique sample
#             if verbose:
#                 print(f"Sampled frame {frame_idx} from replica {repl} for temperature {temperature}.")
#     return sampled_frames  # Return the list of sampled frames

# def validate_coordinates(h5_coords_shape, pdb_num_atoms):
#     """
#     Validates that the number of atoms in the H5 coordinates matches the PDB file.

#     Parameters:
#     -----------
#     h5_coords_shape : tuple
#         Shape of the coordinates array from H5 (numAtoms, 3).
#     pdb_num_atoms : int
#         Number of atoms in the PDB file.

#     Raises:
#     -------
#     ValueError if the number of atoms does not match.
#     """
#     if h5_coords_shape[0] != pdb_num_atoms:
#         raise ValueError(f"Number of atoms in H5 file ({h5_coords_shape[0]}) does not match PDB file ({pdb_num_atoms}).")
#     # If the counts match, the function simply passes without raising an error

# def extract_and_sample(h5_file_path, output_dir, temperature=None, n=10, verbose=False, clean=False):
#     """
#     Main function to extract PDB from H5, sample frames, and generate updated PDB files.

#     Parameters:
#     -----------
#     h5_file_path : str
#         Path to the .h5 file.
#     output_dir : str
#         Directory to save the generated PDB files.
#     temperature : str or None
#         Specific temperature to process. If None, process all available temperatures.
#     n : int
#         Number of frames to sample per temperature.
#     verbose : bool
#         If True, prints detailed logs.
#     clean : bool
#         If True, cleans the PDB files using pdbUtils.
#     """
#     domain_id = parse_domain_id(h5_file_path)  # Extract the domain ID from the H5 filename
#     if verbose:
#         print(f"Domain ID: {domain_id}")

#     # Ensure the output directory exists; create it if it doesn't
#     os.makedirs(output_dir, exist_ok=True)

#     # Open the H5 file for reading
#     try:
#         with h5py.File(h5_file_path, 'r') as h5f:
#             # Extract the PDB string from the H5 file
#             pdb_string = extract_pdb_from_h5(h5f, domain_id, verbose=verbose)

#             # Define the initial PDB file path
#             initial_pdb_path = os.path.join(output_dir, f"{domain_id}.pdb")
            
#             # Write the extracted PDB string to the initial PDB file
#             with open(initial_pdb_path, 'w') as pdb_file:
#                 pdb_file.write(pdb_string)
#             if verbose:
#                 print(f"Extracted PDB written to {initial_pdb_path}")

#             # Clean the PDB if the clean flag is set
#             if clean:
#                 cleaned_pdb_path = os.path.join(output_dir, f"{domain_id}_cleaned.pdb")
#                 # Use pdbUtils to clean the PDB file
#                 pdb_df = pdbUtils.pdb2df(initial_pdb_path)  # Convert PDB to DataFrame
#                 pdbUtils.df2pdb(pdb_df, cleaned_pdb_path)  # Convert DataFrame back to PDB
#                 if verbose:
#                     print(f"Cleaned PDB written to {cleaned_pdb_path}")
#                 pdb_file_path = cleaned_pdb_path  # Use the cleaned PDB for further processing
#             else:
#                 pdb_file_path = initial_pdb_path  # Use the initial PDB for further processing

#             # Load the PDB file lines
#             original_pdb_lines = load_pdb(pdb_file_path)  # Get PDB lines as a list
#             # Count the number of atoms in the PDB file by counting ATOM and HETATM lines
#             pdb_num_atoms = sum(1 for line in original_pdb_lines if line.startswith("ATOM") or line.startswith("HETATM"))
#             if verbose:
#                 print(f"Number of atoms in PDB file: {pdb_num_atoms}")

#             # List all temperature keys that are purely digits (e.g., '320', '348')
#             available_temperatures = [temp for temp in h5f[domain_id].keys() if temp.isdigit()]
#             if temperature:
#                 # If a specific temperature is provided, check its availability
#                 if temperature not in available_temperatures:
#                     raise ValueError(f"Specified temperature {temperature}K not found in the H5 file.")
#                 temperatures_to_process = [temperature]  # Only process the specified temperature
#             else:
#                 temperatures_to_process = available_temperatures  # Process all available temperatures
#                 if verbose:
#                     print(f"No specific temperature specified. Processing all available temperatures: {temperatures_to_process}")

#             # Iterate over each temperature to process
#             for temp in temperatures_to_process:
#                 if verbose:
#                     print(f"\nProcessing temperature: {temp}K")

#                 # Sample 'n' frames from the current temperature
#                 sampled = sample_frames(h5f, domain_id, temp, n, verbose=verbose)

#                 # Iterate over each sampled frame to generate updated PDB files
#                 for idx, (repl, frame_idx) in enumerate(sampled, 1):
#                     coords_dataset = h5f[domain_id][temp][repl]['coords']  # Access the 'coords' dataset

#                     # If verbose, print the shape of the 'coords' dataset for this replica and temperature
#                     if verbose:
#                         print(f"Coords dataset shape for {temp}K, replica {repl}: {coords_dataset.shape}")

#                     # Extract the specific frame's coordinates; shape should be (2165, 3)
#                     coords = coords_dataset[frame_idx, :, :]  # Slice to get all atoms for the selected frame

#                     # Validate that the number of atoms matches between H5 and PDB
#                     validate_coordinates(coords.shape, pdb_num_atoms)

#                     # Create a filename for the new PDB file, indicating domain, temperature, and frame number
#                     frame_filename = f"{domain_id}_{temp}_frame{idx}.pdb"
#                     output_path = os.path.join(output_dir, frame_filename)  # Full path for the new PDB file

#                     # Write the new PDB file with updated coordinates
#                     write_pdb(original_pdb_lines, coords, output_path)

#                     if verbose:
#                         print(f"Generated PDB file: {output_path}")

#     except Exception as e:
#         # Catch and print any errors that occur during processing
#         print(f"[ERROR] {e}")

# def main():
#     # Initialise the argument parser to handle command-line inputs
#     parser = argparse.ArgumentParser(description="Extract and sample frames from mdCATH domain .h5 files to generate PDB files with updated coordinates.")
#     # Define required and optional arguments
#     parser.add_argument("--h5_file", required=True, help="Path to the .h5 file.")
#     parser.add_argument("--output_dir", required=True, help="Directory to save the generated PDB files.")
#     parser.add_argument("--temperature", type=str, default=None, help="Specific temperature to process (e.g., '348'). If not specified, all temperatures will be processed.")
#     parser.add_argument("--n", type=int, default=10, help="Number of frames to sample per temperature.")
#     parser.add_argument("--verbose", action="store_true", help="Increase output verbosity.")
#     parser.add_argument("--clean", action="store_true", help="Clean PDB files using pdbUtils.")

#     args = parser.parse_args()  # Parse the provided command-line arguments

#     # Call the main extraction and sampling function with the parsed arguments
#     extract_and_sample(
#         h5_file_path=args.h5_file,
#         output_dir=args.output_dir,
#         temperature=args.temperature,
#         n=args.n,
#         verbose=args.verbose,
#         clean=args.clean
#     )

# if __name__ == "__main__":
#     main()  # Entry point of the script


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