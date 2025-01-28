#!/usr/bin/env python3
"""
03_extract_RMSF.py

Extract residue-level RMSF data from an mdCATH .h5 file, mirroring
the exact logic in 01_extract_RMSF.py. This script:
  1) Identifies each domain in the .h5 file,
  2) Determines the correct (resid, resname) arrays,
  3) Iterates over all temperature–replica groups,
  4) Extracts and saves RMSF CSVs in 'temperature_data',
  5) Computes per-temperature and overall average RMSF,
  6) Saves averages in 'average_RMSF'.

Usage:
  python 03_extract_RMSF.py \
      --h5_file /path/to/mdcath_dataset.h5 \
      --output_dir /home/s_felix/RMSF_data \
      --verbose
"""

import os
import re
import h5py
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract RMSF data from a single mdCATH .h5 file, fully replicating the 01_extract_RMSF.py logic."
    )
    parser.add_argument("--h5_file", required=True, help="Path to the mdCATH .h5 file.")
    parser.add_argument(
        "--output_dir",
        default="rmsf_output",
        help="Directory for storing RMSF output (default: 'rmsf_output')."
    )
    parser.add_argument("--verbose", action="store_true", help="Print detailed logs.")
    return parser.parse_args()

def create_directory(path):
    """Create a directory if it doesn't already exist."""
    os.makedirs(path, exist_ok=True)

def parse_pdb_for_resid_resname(pdb_string, expected_num, verbose=False):
    """
    Parse PDB lines to create (resid, resname) arrays of length == expected_num.
    """
    lines = pdb_string.splitlines()
    unique_map = []
    seen = set()

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resname = line[17:20].strip()
            try:
                resid = int(line[22:26].strip())
            except ValueError:
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
            print(f"[WARNING] Found {len(unique_map)} unique residues, expected {expected_num}.")
        return None, None

def group_atomwise_resid(domain_group, expected_num, verbose=False):
    """
    If 'resid'/'resname' are per-atom, group them into residue-level arrays.
    """
    if "resid" not in domain_group or "resname" not in domain_group:
        return None, None

    resid_full = domain_group["resid"][:]
    resname_full = domain_group["resname"][:]
    if resid_full.shape[0] <= expected_num:
        return None, None

    unique_resids = []
    unique_resnames = []
    seen = set()

    for i in range(len(resid_full)):
        r = resid_full[i]
        rn = resname_full[i].decode("utf-8")
        if r not in seen:
            seen.add(r)
            unique_resids.append(r)
            unique_resnames.append(rn)

    if len(unique_resids) == expected_num:
        if verbose:
            print(f"[INFO] Grouped {len(unique_resids)} residues from atomwise data.")
        return np.array(unique_resids, dtype=int), unique_resnames
    else:
        if verbose:
            print(f"[WARNING] Grouped {len(unique_resids)} residues, expected {expected_num}.")
        return None, None

def determine_resid_resname(domain_group, verbose=False):
    """
    Replicates the logic from 01_extract_RMSF.py to obtain residue-level arrays.
    """
    if "numResidues" not in domain_group.attrs:
        raise ValueError("[ERROR] 'numResidues' attribute missing.")

    num_residues = domain_group.attrs["numResidues"]

    # Case 1: Direct match
    if "resid" in domain_group and domain_group["resid"].shape[0] == num_residues:
        resid_array = domain_group["resid"][:]
        if "resname" in domain_group and domain_group["resname"].shape[0] == num_residues:
            resname_raw = domain_group["resname"][:]
            resname_array = [x.decode("utf-8") for x in resname_raw]
            if verbose:
                print("[INFO] Found direct residue-level arrays.")
            return resid_array, resname_array

    # Case 2: Parse from PDB
    pdb_string = None
    if "pdbProteinAtoms" in domain_group:
        pdb_string = domain_group["pdbProteinAtoms"][()].decode("utf-8")
    elif "pdb" in domain_group:
        pdb_string = domain_group["pdb"][()].decode("utf-8")

    if pdb_string:
        resid, resname = parse_pdb_for_resid_resname(pdb_string, num_residues, verbose=verbose)
        if resid is not None and resname is not None:
            if verbose:
                print("[INFO] Derived residue arrays from PDB.")
            return resid, resname

    # Case 3: Group from per-atom data
    resid, resname = group_atomwise_resid(domain_group, num_residues, verbose=verbose)
    if resid is not None and resname is not None:
        if verbose:
            print("[INFO] Derived residue arrays by grouping per-atom data.")
        return resid, resname

    raise ValueError("[ERROR] Unable to construct residue-level arrays for domain.")

def main():
    args = parse_arguments()
    create_directory(args.output_dir)

    if not os.path.isfile(args.h5_file):
        raise FileNotFoundError(f"File not found: {args.h5_file}")

    with h5py.File(args.h5_file, "r") as h5f:
        # Identify domains in the .h5 file
        domains = sorted(h5f.keys())
        if args.verbose:
            print(f"[INFO] Found {len(domains)} domain(s) in '{args.h5_file}'.")

        for domain_id in domains:
            if args.verbose:
                print(f"\n[INFO] Processing domain '{domain_id}'.")
            domain_group = h5f[domain_id]

            # Attempt to determine resid, resname
            try:
                resid_data, resname_data = determine_resid_resname(domain_group, verbose=args.verbose)
                num_residues = domain_group.attrs["numResidues"]
            except ValueError as e:
                if args.verbose:
                    print(f"[ERROR] {e} Skipping domain '{domain_id}'.")
                continue

            # Prepare output folders for this domain
            domain_output_dir = os.path.join(args.output_dir, f"rmsf_{domain_id}_data")
            create_directory(domain_output_dir)

            temp_data_dir = os.path.join(domain_output_dir, "temperature_data")
            create_directory(temp_data_dir)

            avg_rmsf_dir = os.path.join(domain_output_dir, "average_RMSF")
            create_directory(avg_rmsf_dir)

            # Collect numeric temperature groups
            temp_groups = sorted([k for k in domain_group.keys() if k.isdigit()])
            if args.verbose:
                print(f"[INFO] Found {len(temp_groups)} temperature group(s).")

            combos = []
            for temp in temp_groups:
                # Sort replicas
                replica_keys = sorted(domain_group[temp].keys(), key=lambda x: (int(x) if x.isdigit() else x))
                for rep in replica_keys:
                    combos.append((temp, rep))

            if args.verbose:
                print(f"[INFO] Found {len(combos)} total temperature–replica combo(s).")

            # Store RMSF data for subsequent averaging
            temp_rmsf_dict = {temp: [] for temp in temp_groups}

            pbar = tqdm(combos, total=len(combos), desc=f"Extracting RMSF for {domain_id}", unit="set")
            for (temp, rep) in pbar:
                replica_group = domain_group[temp][rep]
                if "rmsf" not in replica_group:
                    if args.verbose:
                        print(f"[WARNING] 'rmsf' missing at {domain_id}/{temp}/{rep}. Skipping.")
                    continue

                rmsf_data = replica_group["rmsf"][:]
                if rmsf_data.shape[0] != num_residues:
                    if args.verbose:
                        print(f"[WARNING] RMSF length mismatch at {temp}/{rep}. Skipping.")
                    continue

                # Build DataFrame
                df_rmsf = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    "rmsf": rmsf_data
                })

                # Save per-replica CSV
                out_filename = f"{domain_id}_temperature_{temp}_replica_{rep}_rmsf.csv"
                out_path = os.path.join(temp_data_dir, out_filename)
                df_rmsf.to_csv(out_path, index=False)
                if args.verbose:
                    print(f"[INFO] Saved '{out_filename}' to '{temp_data_dir}'.")

                temp_rmsf_dict[temp].append(rmsf_data)

            # Compute averages per temperature
            for temp, all_rmsf_arrays in temp_rmsf_dict.items():
                if not all_rmsf_arrays:
                    continue
                avg_rmsf = np.mean(np.vstack(all_rmsf_arrays), axis=0)
                df_avg = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    f"rmsf_{temp}": avg_rmsf
                })
                filename_avg = f"{domain_id}_temperature_{temp}_average_rmsf.csv"
                df_avg.to_csv(os.path.join(avg_rmsf_dir, filename_avg), index=False)
                if args.verbose:
                    print(f"[INFO] Saved per-temperature average RMSF '{filename_avg}'.")

            # Compute overall average across all temperatures/replicas
            all_rmsf = []
            for temp, arrs in temp_rmsf_dict.items():
                all_rmsf.extend(arrs)

            if len(all_rmsf) > 0:
                stacked = np.vstack(all_rmsf)
                overall_avg = np.mean(stacked, axis=0)
                df_overall = pd.DataFrame({
                    "protein_id": [domain_id] * num_residues,
                    "resid": resid_data,
                    "resname": resname_data,
                    "average_rmsf": overall_avg
                })
                filename_overall = f"{domain_id}_total_average_rmsf.csv"
                df_overall.to_csv(os.path.join(avg_rmsf_dir, filename_overall), index=False)
                if args.verbose:
                    print(f"[INFO] Saved overall average RMSF '{filename_overall}'.")

if __name__ == "__main__":
    main()
