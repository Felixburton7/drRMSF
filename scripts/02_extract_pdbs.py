#!/usr/bin/env python3
"""
aposteriori_extraction_revised.py

Revised script to extract and sample coordinate frames from mdCATH domain .h5 files,
matching them to specific residues in the original PDB file, and outputting updated PDB files.

Usage:
python aposteriori_extraction_revised.py \
    --h5_file mdcath_dataset_1r9lA02.h5 \
    --output_dir aposteriori_ready \
    --n 10 \
    --verbose \
    [--clean]

Author: Felix
"""

import os  # For interacting with the operating system, like file paths
import re  # For regular expression operations
import h5py  # To handle HDF5 files
import argparse  # For parsing command-line arguments
import numpy as np  # For numerical operations, especially with arrays
from tqdm import tqdm  # For displaying progress bars
from pdbUtils import pdbUtils  # Custom utility for handling PDB files
import random  # For random sampling
import shutil  # For high-level file operations like copying files

def parse_domain_id(h5_filename):
    """
    Extracts the domain ID from an mdCATH .h5 filename.
    Example: 'mdcath_dataset_1r9lA02.h5' -> '1r9lA02'.
    
    Parameters:
    -----------
    h5_filename : str
        The filename of the .h5 file.
    
    Returns:
    --------
    str
        The extracted domain ID.
    """
    basename = os.path.basename(h5_filename)  # Extracts the filename from the full path
    # Example: '/path/to/mdcath_dataset_1r9lA02.h5' -> 'mdcath_dataset_1r9lA02.h5'
    
    # Uses regex to match the pattern and extract the domain ID
    match = re.match(r"mdcath_dataset_(.+)\.h5", basename)
    if not match:
        # If the pattern doesn't match, remove the extension and the prefix manually
        domain_id = os.path.splitext(basename)[0]  # Removes '.h5' extension
        # Example: 'mdcath_dataset_1r9lA02.h5' -> 'mdcath_dataset_1r9lA02'
        
        if domain_id.startswith("mdcath_dataset_"):
            domain_id = domain_id.replace("mdcath_dataset_", "")  # Removes 'mdcath_dataset_' prefix
            # Example: 'mdcath_dataset_1r9lA02' -> '1r9lA02'
        return domain_id  # Returns the cleaned domain ID
    
    return match.group(1)  # Returns the captured domain ID (e.g., '1r9lA02')

def extract_pdb_from_h5(h5_file, domain_id, verbose=False):
    """
    Extracts the PDB string from the H5 file using pdbUtils.

    Parameters:
    -----------
    h5_file : h5py.File
        Opened HDF5 file object.
    domain_id : str
        Domain identifier extracted from the filename.
    verbose : bool
        If True, prints additional processing information.

    Returns:
    --------
    str
        The extracted PDB string.
    
    Raises:
    -------
    KeyError:
        If neither 'pdbProteinAtoms' nor 'pdb' datasets are found within the domain group.
    """
    # Access the top-level group corresponding to the domain ID
    domain_group = h5_file[domain_id]
    # File Structure Example:
    # /1r9lA02/
    # ├── 320/
    # ├── 348/
    # ├── pdbProteinAtoms
    # └── pdb

    # Retrieve the PDB string from 'pdbProteinAtoms' or fallback to 'pdb'
    if "pdbProteinAtoms" in domain_group:
        pdb_string = domain_group["pdbProteinAtoms"][()].decode("utf-8")
        # Example:
        # h5_file['1r9lA02']['pdbProteinAtoms'] contains the PDB data as bytes
        if verbose:
            print(f"Extracted PDB from 'pdbProteinAtoms'.")
    elif "pdb" in domain_group:
        pdb_string = domain_group["pdb"][()].decode("utf-8")
        # Fallback if 'pdbProteinAtoms' is not available
        if verbose:
            print(f"Extracted PDB from 'pdb'.")
    else:
        # If neither PDB field is found, raise an error
        raise KeyError(f"No PDB field found in domain '{domain_id}'.")

    return pdb_string  # Returns the PDB string

def load_pdb(pdb_path):
    """
    Loads the PDB file and returns it as a list of lines.
    
    Parameters:
    -----------
    pdb_path : str
        Path to the PDB file.
    
    Returns:
    --------
    list of str
        List containing each line of the PDB file.
    
    Example:
    --------
    If the PDB file contains:
        ATOM      1  N   ALA A   1      -3.090  -3.392 -18.053  1.00  0.00           N  
        ATOM      2  CA  ALA A   1      -3.248  -1.173 -16.141  1.00  0.00           C  
    The function will return:
        [
            "ATOM      1  N   ALA A   1      -3.090  -3.392 -18.053  1.00  0.00           N  \n",
            "ATOM      2  CA  ALA A   1      -3.248  -1.173 -16.141  1.00  0.00           C  \n",
            ...
        ]
    """
    with open(pdb_path, 'r') as pdb_file:
        pdb_lines = pdb_file.readlines()  # Reads all lines from the PDB file
    return pdb_lines  # Returns the list of PDB lines

def write_pdb(pdb_lines, coords, output_path):
    """
    Writes a new PDB file with updated coordinates.

    Parameters:
    -----------
    pdb_lines : list of str
        Original PDB file lines.
    coords : np.ndarray
        Array of shape (numAtoms, 3) containing updated coordinates.
        Example shape: (2165, 3)
    output_path : str
        Path to save the new PDB file.

    Raises:
    -------
    ValueError:
        If the number of provided coordinates does not match the number of atoms in the PDB.
    
    Example:
    --------
    Original PDB line:
        "ATOM      1  N   ALA A   1      -3.090  -3.392 -18.053  1.00  0.00           N  \n"
    Updated PDB line:
        "ATOM      1  N   ALA A   1      -4.500  -4.500 -19.000  1.00  0.00           N  \n"
    """
    new_pdb_lines = []  # Initializes a list to hold the new PDB lines
    coord_idx = 0  # Index to keep track of the current coordinate

    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Ensure we have enough coordinates to update
            if coord_idx >= coords.shape[0]:
                raise ValueError(f"Insufficient coordinates provided. Expected at least {coord_idx + 1}, got {coords.shape[0]}")

            # Extract the x, y, z coordinates for the current atom
            x, y, z = coords[coord_idx]
            # Create a new line with updated coordinates, maintaining the original formatting
            new_line = (
                f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{line[54:]}"
            )
            new_pdb_lines.append(new_line)  # Add the updated line to the list
            coord_idx += 1  # Move to the next coordinate
        else:
            new_pdb_lines.append(line)  # Non-ATOM/HETATM lines are added unchanged

    # Check if all coordinates have been used
    if coord_idx != coords.shape[0]:
        raise ValueError(f"Unused coordinates detected. Provided {coords.shape[0]}, used {coord_idx}")

    # Write all the new lines to the output PDB file
    with open(output_path, 'w') as out_file:
        out_file.writelines(new_pdb_lines)

def sample_frames(h5_file, domain_id, temperature, n, verbose=False):
    """
    Samples 'n' random frames from any repeat within the specified temperature.

    Parameters:
    -----------
    h5_file : h5py.File
        Opened HDF5 file object.
    domain_id : str
        Domain identifier.
    temperature : str
        Temperature to process.
    n : int
        Number of frames to sample.
    verbose : bool
        If True, prints detailed logs.

    Returns:
    --------
    List of tuples: [(repeat_id, frame_index), ...]
    
    Example:
    --------
    If temperature '320K' has replicas '0', '1', '2', '3', each with 450 frames,
    and n=2, the function might return [('1', 123), ('3', 456)]
    """
    available_replicas = list(h5_file[domain_id][temperature].keys())  # Lists all replicas for the given temperature
    # Example:
    # h5_file['1r9lA02']['320']['0']['coords'].shape -> (450, 2165, 3)
    if verbose:
        print(f"Available replicas ({temperature}K): {available_replicas}")

    sampled_frames = []  # List to store sampled (replica, frame) tuples
    total_frames = 0  # Counter for total available frames

    # Iterate over each replica to count total frames
    for repl in available_replicas:
        num_frames = h5_file[domain_id][temperature][repl]['coords'].shape[0]  # Number of frames in this replica
        # Example: (450, 2165, 3) -> shape[0] = 450 frames
        total_frames += num_frames  # Update total frame count
        if verbose:
            print(f"Replica {repl} has {num_frames} frames.")

    # Check if there are enough frames to sample
    if total_frames < n:
        raise ValueError(f"Not enough frames ({total_frames}) available to sample {n} frames for temperature {temperature}.")

    # Randomly sample frames until we have 'n' unique samples
    while len(sampled_frames) < n:
        repl = random.choice(available_replicas)  # Randomly select a replica
        num_frames = h5_file[domain_id][temperature][repl]['coords'].shape[0]  # Frames in the selected replica
        frame_idx = random.randint(0, num_frames - 1)  # Random frame index within the replica
        if (repl, frame_idx) not in sampled_frames:
            sampled_frames.append((repl, frame_idx))  # Add unique sample
            if verbose:
                print(f"Sampled frame {frame_idx} from replica {repl} for temperature {temperature}.")
    
    # Example sampled_frames after sampling:
    # [('1', 123), ('3', 456), ...]
    return sampled_frames  # Return the list of sampled frames

def validate_coordinates(h5_coords_shape, pdb_num_atoms):
    """
    Validates that the number of atoms in the H5 coordinates matches the PDB file.

    Parameters:
    -----------
    h5_coords_shape : tuple
        Shape of the coordinates array from H5 (numAtoms, 3).
        Example: (2165, 3)
    pdb_num_atoms : int
        Number of atoms in the PDB file.
        Example: 2165

    Raises:
    -------
    ValueError:
        If the number of atoms does not match.
    """
    if h5_coords_shape[0] != pdb_num_atoms:
        raise ValueError(f"Number of atoms in H5 file ({h5_coords_shape[0]}) does not match PDB file ({pdb_num_atoms}).")
    # If the counts match, the function simply passes without raising an error

def extract_and_sample(h5_file_path, output_dir, temperature=None, n=10, verbose=False, clean=False):
    """
    Main function to extract PDB from H5, sample frames, and generate updated PDB files.

    Parameters:
    -----------
    h5_file_path : str
        Path to the .h5 file.
        Example: 'mdcath_dataset_1r9lA02.h5'
    output_dir : str
        Directory to save the generated PDB files.
        Example: 'aposteriori_ready'
    temperature : str or None
        Specific temperature to process. If None, process all available temperatures.
        Example: '320'
    n : int
        Number of frames to sample per temperature.
        Example: 10
    verbose : bool
        If True, prints detailed logs.
    clean : bool
        If True, cleans the PDB files using pdbUtils.
    
    Process Overview:
    -----------------
    1. Parse the domain ID from the .h5 filename.
    2. Extract the PDB string from the .h5 file and write it to an initial PDB file.
    3. Optionally clean the PDB file using pdbUtils.
    4. Load the PDB file and count the number of atoms.
    5. Identify temperatures to process (all or specific).
    6. For each temperature:
        a. Sample 'n' frames from available replicas.
        b. For each sampled frame:
            i. Retrieve the corresponding coordinates.
            ii. Validate atom count.
            iii. Write a new PDB file with updated coordinates.
    
    """
    domain_id = parse_domain_id(h5_file_path)  # Extract the domain ID from the H5 filename
    if verbose:
        print(f"Domain ID: {domain_id}")

    # Ensure the output directory exists; create it if it doesn't
    os.makedirs(output_dir, exist_ok=True)

    # Open the H5 file for reading
    try:
        with h5py.File(h5_file_path, 'r') as h5f:
            # Extract the PDB string from the H5 file
            pdb_string = extract_pdb_from_h5(h5f, domain_id, verbose=verbose)

            # Define the initial PDB file path
            initial_pdb_path = os.path.join(output_dir, f"{domain_id}.pdb")
            # Example:
            # output_dir = 'aposteriori_ready'
            # domain_id = '1r9lA02'
            # initial_pdb_path = 'aposteriori_ready/1r9lA02.pdb'

            # Write the extracted PDB string to the initial PDB file
            with open(initial_pdb_path, 'w') as pdb_file:
                pdb_file.write(pdb_string)
            if verbose:
                print(f"Extracted PDB written to {initial_pdb_path}")
                # Example Output:
                # "Extracted PDB written to aposteriori_ready/1r9lA02.pdb"

            # Clean the PDB if the clean flag is set
            if clean:
                cleaned_pdb_path = os.path.join(output_dir, f"{domain_id}_cleaned.pdb")
                # Example: 'aposteriori_ready/1r9lA02_cleaned.pdb'

                # Use pdbUtils to clean the PDB file
                pdb_df = pdbUtils.pdb2df(initial_pdb_path)  # Convert PDB to DataFrame
                # Example:
                # pdb_df = pdbUtils.pdb2df('aposteriori_ready/1r9lA02.pdb') 
                
                pdbUtils.df2pdb(pdb_df, cleaned_pdb_path)  # Convert DataFrame back to PDB
                # Example:
                # pdbUtils.df2pdb(pdb_df, 'aposteriori_ready/1r9lA02_cleaned.pdb')
                
                if verbose:
                    print(f"Cleaned PDB written to {cleaned_pdb_path}")
                    # Example Output:
                    # "Cleaned PDB written to aposteriori_ready/1r9lA02_cleaned.pdb"
                
                pdb_file_path = cleaned_pdb_path  # Use the cleaned PDB for further processing
            else:
                pdb_file_path = initial_pdb_path  # Use the initial PDB for further processing

            # Load the PDB file lines
            original_pdb_lines = load_pdb(pdb_file_path)  # Get PDB lines as a list
            # Example:
            # original_pdb_lines = [
            #     "ATOM      1  N   ALA A   1      -3.090  -3.392 -18.053  1.00  0.00           N  \n",
            #     "ATOM      2  CA  ALA A   1      -3.248  -1.173 -16.141  1.00  0.00           C  \n",
            #     ...
            # ]

            # Count the number of atoms in the PDB file by counting ATOM and HETATM lines
            pdb_num_atoms = sum(1 for line in original_pdb_lines if line.startswith("ATOM") or line.startswith("HETATM"))
            # Example:
            # If there are 2165 ATOM/HETATM lines, pdb_num_atoms = 2165
            if verbose:
                print(f"Number of atoms in PDB file: {pdb_num_atoms}")
                # Example Output:
                # "Number of atoms in PDB file: 2165"

            # List all temperature keys that are purely digits (e.g., '320', '348')
            available_temperatures = [temp for temp in h5f[domain_id].keys() if temp.isdigit()]
            # Example:
            # available_temperatures = ['320', '348', '379', '413', '450']

            if temperature:
                # If a specific temperature is provided, check its availability
                if temperature not in available_temperatures:
                    raise ValueError(f"Specified temperature {temperature}K not found in the H5 file.")
                temperatures_to_process = [temperature]  # Only process the specified temperature
                # Example:
                # temperature = '320'
                # temperatures_to_process = ['320']
            else:
                temperatures_to_process = available_temperatures  # Process all available temperatures
                if verbose:
                    print(f"No specific temperature specified. Processing all available temperatures: {temperatures_to_process}")
                    # Example Output:
                    # "No specific temperature specified. Processing all available temperatures: ['320', '348', '379', '413', '450']"

            # Iterate over each temperature to process
            for temp in temperatures_to_process:
                if verbose:
                    print(f"\nProcessing temperature: {temp}K")
                    # Example Output:
                    # "\nProcessing temperature: 320K"

                # Sample 'n' frames from the current temperature
                sampled = sample_frames(h5f, domain_id, temp, n, verbose=verbose)
                # Example:
                # sampled = [('1', 123), ('3', 456), ...]

                # Iterate over each sampled frame to generate updated PDB files
                for idx, (repl, frame_idx) in enumerate(sampled, 1):
                    coords_dataset = h5f[domain_id][temp][repl]['coords']  # Access the 'coords' dataset
                    # Example:
                    # coords_dataset = h5_file['1r9lA02']['320']['1']['coords']
                    # coords_dataset.shape -> (450, 2165, 3)
                    # Meaning:
                    # - 450 frames
                    # - 2165 atoms per frame
                    # - 3 coordinates (x, y, z) per atom

                    # If verbose, print the shape of the 'coords' dataset for this replica and temperature
                    if verbose:
                        print(f"Coords dataset shape for {temp}K, replica {repl}: {coords_dataset.shape}")
                        # Example Output:
                        # "Coords dataset shape for 320K, replica 1: (450, 2165, 3)"

                    # Extract the specific frame's coordinates; shape should be (2165, 3)
                    coords = coords_dataset[frame_idx, :, :]  # Slice to get all atoms for the selected frame
                    # Example:
                    # coords.shape -> (2165, 3)
                    # coords[0] -> [x1, y1, z1]
                    # coords[1] -> [x2, y2, z2]
                    # ...

                    # Validate that the number of atoms matches between H5 and PDB
                    validate_coordinates(coords.shape, pdb_num_atoms)
                    # Ensures that coords.shape[0] == pdb_num_atoms

                    # Create a filename for the new PDB file, indicating domain, temperature, and frame number
                    frame_filename = f"{domain_id}_{temp}_frame{idx}.pdb"
                    # Example:
                    # frame_filename = '1r9lA02_320_frame1.pdb'

                    output_path = os.path.join(output_dir, frame_filename)  # Full path for the new PDB file
                    # Example:
                    # output_path = 'aposteriori_ready/1r9lA02_320_frame1.pdb'

                    # Write the new PDB file with updated coordinates
                    write_pdb(original_pdb_lines, coords, output_path)
                    # The write_pdb function updates the coordinates in the PDB lines and writes to output_path

                    if verbose:
                        print(f"Generated PDB file: {output_path}")
                        # Example Output:
                        # "Generated PDB file: aposteriori_ready/1r9lA02_320_frame1.pdb"

    except Exception as e:
        # Catch and print any errors that occur during processing
        print(f"[ERROR] {e}")

def main():
    # Initialise the argument parser to handle command-line inputs
    parser = argparse.ArgumentParser(description="Extract and sample frames from mdCATH domain .h5 files to generate PDB files with updated coordinates.")
    # Define required and optional arguments
    parser.add_argument("--h5_file", required=True, help="Path to the .h5 file.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the generated PDB files.")
    parser.add_argument("--temperature", type=str, default=None, help="Specific temperature to process (e.g., '348'). If not specified, all temperatures will be processed.")
    parser.add_argument("--n", type=int, default=10, help="Number of frames to sample per temperature.")
    parser.add_argument("--verbose", action="store_true", help="Increase output verbosity.")
    parser.add_argument("--clean", action="store_true", help="Clean PDB files using pdbUtils.")

    args = parser.parse_args()  # Parse the provided command-line arguments

    # Call the main extraction and sampling function with the parsed arguments
    extract_and_sample(
        h5_file_path=args.h5_file,
        output_dir=args.output_dir,
        temperature=args.temperature,
        n=args.n,
        verbose=args.verbose,
        clean=args.clean
    )

if __name__ == "__main__":
    main()  # Entry point of the script



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