
#!/usr/bin/env python3
"""
01_extract_RMSF.py

Extract residue-level RMSF data from an mdCATH .h5 file, across all
temperature/replica groups. If the domain stores 'resid' and 'resname'
per-atom (e.g., 2165 entries) rather than per-residue (141), we group
by residue to produce a unique mapping of length == numResidues.

Alternatively, if grouping fails or is undesirable, we can parse the PDB
'pdbProteinAtoms' to reconstruct (resid, resname) in PDB order,
matching the domain's 'numResidues' attribute.

**Enhancements:**
1. Added 'protein_id' column to CSVs.
2. Computed per-temperature averages across replicates.
3. Adjusted folder structure to include domain-specific directories.
4. Modified script to process all domains within the .h5 file.

Workflow:
  1) Identify the domain's numResidues attribute (e.g., 141).
  2) Attempt to build residue-level arrays by grouping or by parsing the PDB.
  3) Loop over each temperature (e.g., '320', '348', etc.) and each replica.
  4) Extract 'rmsf' (must be length == numResidues).
  5) Save CSV in a temperature_<temp> folder.
  6) Compute per-temperature and overall average RMSF, saving to 'average_RMSF'.

Example usage:
  python extract_RMSF.py \
    --h5_file mdcath_dataset.h5 \
    --output_dir ./rmsf_output \
    --verbose

Author: Senior Data Scientist
"""

import os
import re
import h5py
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

def parse_arguments():
    """
    Command-line arguments:
      --h5_file   : Path to the mdCATH .h5 file (e.g., 'mdcath_dataset.h5').
      --output_dir: Base directory for storing RMSF CSV files.
      --verbose   : If set, prints additional logs and warnings.
    """
    parser = argparse.ArgumentParser(
        description="Extract RMSF data from an mdCATH .h5 file, ensuring correct residue names and indices."
    )
    parser.add_argument(
        "--h5_file",
        required=True,
        help="Path to the mdCATH .h5 file (e.g., 'mdcath_dataset.h5')."
    )
    parser.add_argument(
        "--output_dir",
        default="rmsf_output",
        help="Base directory for storing RMSF CSV files (default: 'rmsf_output')."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="If set, prints logs and warnings."
    )
    return parser.parse_args()

def create_directory(path):
    """Create an output directory if it doesn't already exist."""
    os.makedirs(path, exist_ok=True)

def parse_pdb_for_resid_resname(pdb_string, expected_num, verbose=False):
    """
    Parse 'pdbProteinAtoms' lines to extract a unique (resid, resname) list
    in the order of first appearance. Returns arrays of shape (expected_num,)
    if successful, else None.

    Parameters
    ----------
    pdb_string : str
        Full PDB text from 'pdbProteinAtoms' or 'pdb' dataset.
    expected_num : int
        The domain's numResidues attribute (e.g., 141).
    verbose : bool

    Returns
    -------
    (resid_array, resname_array) or None if we can't produce exactly 'expected_num'.
    """
    lines = pdb_string.splitlines()
    unique_map = []  # list of (resid_int, resname_str)
    seen = set()
    for line in lines:
        # Typical ATOM/HETATM line format:
        # Columns:
        #  1-6   Record name "ATOM  " or "HETATM"
        #  7-11  Atom serial number
        #  13-16 Atom name
        #  17    Alternate location indicator
        # 18-20 Residue name
        # 22    Chain identifier
        # 23-26 Residue sequence number
        # 27    Code for insertion of residues
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            try:
                resid = int(line[22:26].strip())
            except ValueError:
                # Malformed line, skip
                continue
            if resid not in seen:
                seen.add(resid)
                unique_map.append((resid, resname))
    if len(unique_map) == expected_num:
        if verbose:
            print(f"[INFO] Extracted {len(unique_map)} unique residues from PDB lines.")
        resid_array = np.array([p[0] for p in unique_map], dtype=int)
        resname_array = [p[1] for p in unique_map]
        return resid_array, resname_array
    else:
        if verbose:
            print(f"[WARNING] Found {len(unique_map)} unique residues in PDB, but expected {expected_num}.")
        return None

def group_atomwise_resid(domain_group, expected_num, verbose=False):
    """
    If 'resid' and 'resname' each have e.g., 2165 entries (per-atom),
    group them by residue ID. Each unique residue ID => one entry
    in the final arrays. The order is the order of first appearance.

    Returns
    -------
    (resid_array, resname_array) of length 'expected_num' or None
    """
    if "resid" not in domain_group or "resname" not in domain_group:
        if verbose:
            print("[WARNING] 'resid' or 'resname' dataset missing.")
        return None
    resid_full = domain_group["resid"][:]      # e.g., shape (2165,)
    resname_full = domain_group["resname"][:]
    if resid_full.shape[0] <= expected_num:
        # This function only relevant if we have more than 'expected_num'
        return None

    # We'll do stable grouping by the order of appearance
    unique_resids = []
    unique_resnames = []
    seen = set()

    for i in range(len(resid_full)):
        r = resid_full[i]
        # decode resname from bytes
        rn = resname_full[i].decode("utf-8") if hasattr(resname_full[i], 'decode') else str(resname_full[i])
        if r not in seen:
            seen.add(r)
            unique_resids.append(r)
            unique_resnames.append(rn)

    if len(unique_resids) == expected_num:
        if verbose:
            print(f"[INFO] Grouped per-atom arrays => {len(unique_resids)} unique residues.")
        return np.array(unique_resids, dtype=int), unique_resnames
    else:
        if verbose:
            print(f"[WARNING] Grouping found {len(unique_resids)} unique residues, not {expected_num}.")
        return None

def determine_resid_resname(domain_group, verbose=False):
    """
    Attempt to produce residue-level arrays (resid, resname) of length
    domain_group.attrs['numResidues'].

    The logic:
      1) If domain_group['resid'].shape == numResidues => use them directly.
      2) Else if 'pdbProteinAtoms' is present => parse PDB lines for unique residue info
         of length numResidues.
      3) Else attempt grouping from the 2165-atom arrays.
      4) If all fails, raise an error.

    Returns
    -------
    resid_array : np.ndarray
    resname_array : list of str
    """
    if "numResidues" not in domain_group.attrs:
        raise ValueError("[ERROR] 'numResidues' attribute missing. Cannot proceed.")

    num_residues = domain_group.attrs["numResidues"]

    # 1) If we already have domain_group['resid'].shape == (num_residues,)
    if "resid" in domain_group and domain_group["resid"].shape[0] == num_residues:
        # Perfect case
        resid_array = domain_group["resid"][:]
        # decode resname
        if "resname" in domain_group and domain_group["resname"].shape[0] == num_residues:
            resname_raw = domain_group["resname"][:]
            resname_array = [x.decode("utf-8") for x in resname_raw]
            if verbose:
                print("[INFO] Found residue-level arrays directly matching numResidues.")
            return resid_array, resname_array

    # 2) Try parsing from PDB if present
    pdb_string = None
    if "pdbProteinAtoms" in domain_group:
        pdb_string = domain_group["pdbProteinAtoms"][()].decode("utf-8")
    elif "pdb" in domain_group:
        pdb_string = domain_group["pdb"][()].decode("utf-8")

    if pdb_string is not None:
        parsed = parse_pdb_for_resid_resname(pdb_string, num_residues, verbose=verbose)
        if parsed is not None:
            if verbose:
                print("[INFO] Successfully derived residue arrays from PDB.")
            return parsed

    # 3) Attempt grouping from 2165-atom arrays
    grouped = group_atomwise_resid(domain_group, num_residues, verbose=verbose)
    if grouped is not None:
        if verbose:
            print("[INFO] Successfully grouped per-atom data into residue-level arrays.")
        return grouped

    # 4) If we reach here, nothing worked
    raise ValueError("[ERROR] Unable to construct residue-level arrays of length == numResidues.")

def main():
    args = parse_arguments()

    # Base output directory
    create_directory(args.output_dir)

    # Open HDF5
    if not os.path.isfile(args.h5_file):
        raise FileNotFoundError(f"File not found: {args.h5_file}")

    with h5py.File(args.h5_file, "r") as h5f:
        # Identify all domains in the .h5 file
        domains = [d for d in h5f.keys()]
        if args.verbose:
            print(f"[INFO] Found {len(domains)} domain(s) in the .h5 file.")

        for domain_id in domains:
            if args.verbose:
                print(f"\n[INFO] Processing domain '{domain_id}'.")

            domain_group = h5f[domain_id]

            try:
                # Determine (resid, resname)
                resid_data, resname_data = determine_resid_resname(domain_group, verbose=args.verbose)
                num_residues = domain_group.attrs["numResidues"]
            except (KeyError, ValueError) as e:
                if args.verbose:
                    print(f"[ERROR] {e} Skipping domain '{domain_id}'.")
                continue  # Skip this domain due to inconsistency

            # Create domain-specific output directory
            domain_output_dir = os.path.join(args.output_dir, f"rmsf_{domain_id}_data")
            create_directory(domain_output_dir)

            # Subdirectories
            temp_data_dir = os.path.join(domain_output_dir, "temperature_data")
            create_directory(temp_data_dir)

            avg_rmsf_dir = os.path.join(domain_output_dir, "average_RMSF")
            create_directory(avg_rmsf_dir)

            # Identify numeric temperature groups
            temp_groups = sorted([k for k in domain_group.keys() if k.isdigit()])
            if len(temp_groups) == 0 and args.verbose:
                print(f"[WARNING] No numeric temperature groups found in domain '{domain_id}'.")

            # Gather all (temp, replica) combinations
            combos = []
            for temp in temp_groups:
                # Sort replicas numerically/alphabetically
                replica_keys = sorted(domain_group[temp].keys(), key=lambda x: (int(x) if x.isdigit() else x))
                for rep in replica_keys:
                    combos.append((temp, rep))

            if args.verbose:
                print(f"[INFO] Found {len(combos)} temperature–replica combos in domain '{domain_id}'.")

            # Initialize structures to hold per-temperature RMSF data
            temp_rmsf_dict = {temp: [] for temp in temp_groups}

            # Progress bar for extraction
            pbar = tqdm(
                combos,
                total=len(combos),
                desc=f"Extracting RMSF for {domain_id}",
                unit="set",
                bar_format="{l_bar}{bar}| {percentage:5.1f}% [{elapsed}<{remaining}]",
                dynamic_ncols=True
            )

            for (temp, rep) in pbar:
                replica_group = domain_group[temp][rep]

                # Check for 'rmsf'
                if "rmsf" not in replica_group:
                    if args.verbose:
                        print(f"[WARNING] 'rmsf' missing at {domain_id}/{temp}/{rep}. Skipping.")
                    continue

                # Extract 'rmsf' and validate shape
                rmsf_data = replica_group["rmsf"][:]
                if rmsf_data.shape[0] != num_residues:
                    if args.verbose:
                        print(f"[WARNING] 'rmsf' length ({rmsf_data.shape[0]}) != numResidues ({num_residues}) at {temp}/{rep}. Skipping.")
                    continue

                # Build DataFrame with 'protein_id' column
                df_rmsf = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    "rmsf": rmsf_data
                })

                # Save to temperature_data folder
                out_filename = f"{domain_id}_temperature_{temp}_replica_{rep}_rmsf.csv"
                out_path = os.path.join(temp_data_dir, out_filename)
                df_rmsf.to_csv(out_path, index=False)

                if args.verbose:
                    print(f"[INFO] Saved '{out_filename}' to '{temp_data_dir}'.")

                # Accumulate RMSF for per-temperature averaging
                temp_rmsf_dict[temp].append(rmsf_data)

                   # After processing all replicates, compute per-temperature averages
            for temp, rmsf_list in temp_rmsf_dict.items():
                if len(rmsf_list) == 0:
                    if args.verbose:
                        print(f"[WARNING] No valid RMSF data found for temperature {temp} in domain '{domain_id}'.")
                    continue
                avg_rmsf = np.mean(np.vstack(rmsf_list), axis=0)
                df_avg_temp = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    f"rmsf_{temp}": avg_rmsf  # Renamed column
                })
                avg_temp_filename = f"{domain_id}_temperature_{temp}_average_rmsf.csv"
                avg_temp_path = os.path.join(avg_rmsf_dir, avg_temp_filename)
                df_avg_temp.to_csv(avg_temp_path, index=False)

                if args.verbose:
                    print(f"[INFO] Saved per-temperature average RMSF '{avg_temp_filename}' to '{avg_rmsf_dir}'.")

            # Compute overall average RMSF across all temperatures and replicates
            all_rmsf = []
            for rmsf_list in temp_rmsf_dict.values():
                all_rmsf.extend(rmsf_list)

            if len(all_rmsf) > 0:
                stacked_all = np.vstack(all_rmsf)  # Shape: (total_replicates, num_residues)
                overall_avg_rmsf = np.mean(stacked_all, axis=0)
                df_overall_avg = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    "rmsf_total": overall_avg_rmsf  # Renamed column
                })
                overall_avg_filename = f"{domain_id}_total_average_rmsf.csv"
                overall_avg_path = os.path.join(avg_rmsf_dir, overall_avg_filename)
                df_overall_avg.to_csv(overall_avg_path, index=False)

                if args.verbose:
                    print(f"[INFO] Saved overall average RMSF '{overall_avg_filename}' to '{avg_rmsf_dir}'.")
            else:
                if args.verbose:
                    print(f"[WARNING] No valid RMSF data found for domain '{domain_id}'; overall average not computed.")

            # Compute overall average RMSF across all temperatures and replicates
            all_rmsf = []
            for rmsf_list in temp_rmsf_dict.values():
                all_rmsf.extend(rmsf_list)

            if len(all_rmsf) > 0:
                stacked_all = np.vstack(all_rmsf)  # Shape: (total_replicates, num_residues)
                overall_avg_rmsf = np.mean(stacked_all, axis=0)
                df_overall_avg = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    "average_rmsf": overall_avg_rmsf
                })
                overall_avg_filename = f"{domain_id}_total_average_rmsf.csv"
                overall_avg_path = os.path.join(avg_rmsf_dir, overall_avg_filename)
                df_overall_avg.to_csv(overall_avg_path, index=False)

                if args.verbose:
                    print(f"[INFO] Saved overall average RMSF '{overall_avg_filename}' to '{avg_rmsf_dir}'.")
            else:
                if args.verbose:
                    print(f"[WARNING] No valid RMSF data found for domain '{domain_id}'; overall average not computed.")

if __name__ == "__main__":
    main()

''' 

5 Histograms of the RMSFs 
- One label to predict. 
- Visuallize the RMSF - which temperature we want as the Y-value for the dataset. 
- Choose the temperature that we want. n frames = 1 to begin with. 
- 

'''


''' 

Null Hypothesis: 
- Is more subtle 
- 

Critism: 
- Yeah I know - Highly flexible regions. 
- N

Response: 
- Simple ML algorithm shouldn't cover it. 
- Featurize it very simply
- Conventional logic to produce features 
- Extreme gradient boosting to see whether we can predict it. 
- CNN's are better. 

Show: 
- Progression of bar charts 
- This one is the worst (random forest - boosted)
- Ours is better! 
- 

NOW: 
- Multi-processing.pool
- For each worker. Make it output a csv. 
- Don't bother with result = funciton. 
- Get it to return nothing and just return a file. 
- If output file exists return STOP. 
- 
'''


'''

- Logical FLOW: 
1) tmux - run the scripts overnight to get the csv files. 
2) and the other script as a PDB file. 
3) 


'''