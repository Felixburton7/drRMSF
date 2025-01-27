#!/usr/bin/env python3
"""
validate_sample_extraction.py

Generates a representative sample of structural and RMSF data from mdCATH
domain .h5 files to verify that the extraction functions operate correctly.
This helps to avoid processing the entire dataset when doing a quick check.

Usage:
  python validate_sample_extraction.py \
      --input_dir /path/to/mdCATH \
      --output_dir /path/to/sampled_data \
      --sample_size 100 \
      --sample_fraction_coords 0.01 \
      [--verbose]

Author: Researcher / Senior Scientist
"""

import os
import re
import glob
import h5py
import yaml
import json
import argparse
import random
import numpy as np
from tqdm import tqdm
import pandas as pd

def parse_domain_id(h5_filename):
    """
    Extracts the domain ID from an mdCATH .h5 filename.
    Example: 'mdcath_dataset_1r9lA02.h5' -> '1r9lA02'.
    """
    basename = os.path.basename(h5_filename)
    match = re.match(r"mdcath_dataset_(.+)\.h5", basename)
    if not match:
        return os.path.splitext(basename)[0]
    return match.group(1)

def sample_coords(coords, fraction):
    """
    Randomly samples a fraction of the coordinate frames (F x N x 3).
    Returns the downsampled coordinate array.
    """
    if coords is None:
        return None
    if fraction <= 0 or fraction >= 1:
        return coords  # no sampling
    F = coords.shape[0]
    sample_count = max(1, int(F * fraction))
    sampled_indices = sorted(random.sample(range(F), sample_count))
    return coords[sampled_indices, :, :]

def sample_array(array, fraction):
    """
    Randomly samples a fraction of frames (1D array shape: F).
    Returns the downsampled array.
    """
    if array is None:
        return None
    if fraction <= 0 or fraction >= 1:
        return array
    F = array.shape[0]
    sample_count = max(1, int(F * fraction))
    sampled_indices = sorted(random.sample(range(F), sample_count))
    return array[sampled_indices]

def validate_aposteriori_extraction(
    input_path,
    output_dir,
    sample_size=100,
    sample_fraction_coords=0.01,
    verbose=False
):
    """
    Main function for extracting a validation sample from mdCATH domain .h5 files.

    Parameters:
    -----------
    input_path : str
        Path to a single .h5 file or a directory containing .h5 files.
    output_dir : str
        Directory to which the sampled data will be saved.
    sample_size : int
        Number of domain files to sample from.
    sample_fraction_coords : float
        Fraction of coordinates and frames to sample within each selected domain.
    verbose : bool
        If True, prints additional information.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Gather files
    if os.path.isdir(input_path):
        domain_files = glob.glob(os.path.join(input_path, "**/*.h5"), recursive=True)
    else:
        domain_files = [input_path]

    if len(domain_files) == 0:
        if verbose:
            print("[WARNING] No .h5 files found. Nothing to do.")
        return

    # Select a random sample of domain_files
    if sample_size < len(domain_files):
        selected_files = random.sample(domain_files, sample_size)
    else:
        selected_files = domain_files  # if sample_size >= total files

    if verbose:
        print(f"Selected {len(selected_files)} domain file(s) for validation sample.")

    # Process each selected file
    for domain_file in tqdm(selected_files, desc="Sampling Data", unit="file"):
        domain_id = parse_domain_id(domain_file)

        try:
            with h5py.File(domain_file, "r") as h5f:
                top_keys = list(h5f.keys())
                if len(top_keys) == 0:
                    if verbose:
                        print(f"[WARNING] No groups found in {domain_file}. Skipping.")
                    continue

                sample_data = {
                    "domain_id": domain_id,
                    "sampled_structures": []
                }

                for dkey in top_keys:
                    group = h5f[dkey]

                    # Decode 'chain' array if present (avoid bytes in JSON)
                    if "chain" in group:
                        chain_bytes = group["chain"][:]
                        chain_list = chain_bytes.tolist()
                        # If any item is bytes, decode it
                        if chain_list and isinstance(chain_list[0], bytes):
                            chain_list = [c.decode("utf-8") for c in chain_list]
                    else:
                        chain_list = []

                    # Decode 'resname' if present
                    if "resname" in group:
                        resname_list = [r.decode("utf-8") for r in group["resname"][:]]
                    else:
                        resname_list = []

                    # Minimal info about domain
                    dom_info = {
                        "dkey": dkey,
                        "numResidues": int(group.attrs["numResidues"]) if "numResidues" in group.attrs else None,
                        "chains": chain_list,
                        "resid": group["resid"][:].tolist() if "resid" in group else [],
                        "resname": resname_list,
                        "temp_groups": []
                    }

                    # Identify temperature groups
                    temp_groups = [tg for tg in group.keys() if tg.isdigit()]
                    for tg in temp_groups:
                        temp_group_data = {
                            "temperature": tg,
                            "replicas": []
                        }
                        # For each replica
                        rep_groups = list(group[tg].keys())
                        for rg in rep_groups:
                            replica_group = group[tg][rg]
                            # Load arrays if present
                            coords = replica_group["coords"][:] if "coords" in replica_group else None
                            dssp = replica_group["dssp"][:] if "dssp" in replica_group else None
                            gyration = replica_group["gyrationRadius"][:] if "gyrationRadius" in replica_group else None
                            rmsd = replica_group["rmsd"][:] if "rmsd" in replica_group else None
                            rmsf = replica_group["rmsf"][:] if "rmsf" in replica_group else None

                            # Sample data to keep the file size small
                            sampled_coords = sample_coords(coords, sample_fraction_coords)
                            sampled_dssp = sample_array(dssp, sample_fraction_coords)
                            sampled_gyration = sample_array(gyration, sample_fraction_coords)
                            sampled_rmsd = sample_array(rmsd, sample_fraction_coords)

                            rep_info = {
                                "replica": rg,
                                "coords_sample_shape": sampled_coords.shape if sampled_coords is not None else None,
                                "dssp_sample_shape": sampled_dssp.shape if sampled_dssp is not None else None,
                                "gyration_sample_shape": sampled_gyration.shape if sampled_gyration is not None else None,
                                "rmsd_sample_shape": sampled_rmsd.shape if sampled_rmsd is not None else None,
                                "rmsf": rmsf.tolist() if rmsf is not None else None
                            }
                            temp_group_data["replicas"].append(rep_info)

                        dom_info["temp_groups"].append(temp_group_data)
                    sample_data["sampled_structures"].append(dom_info)

                # Save the sampled data as JSON
                out_filename = f"{domain_id}_sample.json"
                out_path = os.path.join(output_dir, out_filename)
                with open(out_path, "w") as jsout:
                    json.dump(sample_data, jsout, indent=2)

                if verbose:
                    print(f"Saved sampled data to {out_path}")

        except Exception as e:
            if verbose:
                print(f"[ERROR] Unable to sample from '{domain_file}': {e}")

def main():
    parser = argparse.ArgumentParser(description="Validate mdCATH data extraction by sampling a fraction of the dataset.")
    parser.add_argument("--input_dir", required=True, help="Path to a single .h5 file or a directory of .h5 files.")
    parser.add_argument("--output_dir", required=True, help="Directory where sampled data files will be saved.")
    parser.add_argument("--sample_size", type=int, default=100, help="Number of domains to sample.")
    parser.add_argument("--sample_fraction_coords", type=float, default=0.01, help="Fraction of frames to sample.")
    parser.add_argument("--verbose", action="store_true", help="Increase output verbosity.")
    args = parser.parse_args()

    # Ensure reproducible sampling by fixing a seed if desired.
    random.seed(42)

    validate_aposteriori_extraction(
        input_path=args.input_dir,
        output_dir=args.output_dir,
        sample_size=args.sample_size,
        sample_fraction_coords=args.sample_fraction_coords,
        verbose=args.verbose
    )

if __name__ == "__main__":
    main()
