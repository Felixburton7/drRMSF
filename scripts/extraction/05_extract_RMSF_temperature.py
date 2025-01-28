#!/usr/bin/env python3
"""
extract_RMSF_413.py

This script scans all rmsf_*_data directories under ~/RMSF_data/ and copies only 
the CSV files corresponding to temperature_413_average_rmsf.csv into 
~/RMSF_data/RMSF_413_data/average_RMSF/.

Usage:
  python extract_RMSF_413.py \
    --rmsf_root_dir ~/RMSF_data \
    --output_dir ~/RMSF_data/RMSF_413_data/average_RMSF
    [--verbose]

Author: YourName
"""

import os
import shutil
import argparse
from glob import glob

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Extract and copy RMSF CSV files for temperature 413."
    )
    parser.add_argument(
        "--rmsf_root_dir",
        default="~/RMSF_data",
        help="Root directory containing rmsf_*_data subdirectories. Default: ~/RMSF_data"
    )
    parser.add_argument(
        "--output_dir",
        default="~/RMSF_data/RMSF_413_data/average_RMSF",
        help="Destination directory for copied 413 CSV files."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed output."
    )
    return parser.parse_args()

def main():
    args = parse_arguments()
    rmsf_root_dir = os.path.expanduser(args.rmsf_root_dir)
    output_dir = os.path.expanduser(args.output_dir)

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Find all directories of the form rmsf_*_data under RMSF root
    all_data_dirs = glob(os.path.join(rmsf_root_dir, "rmsf_*_data"))
    
    if args.verbose:
        print(f"[INFO] Found {len(all_data_dirs)} directories matching 'rmsf_*_data' under {rmsf_root_dir}.")

    copied_files_count = 0

    # Loop over each domain-specific directory
    for data_dir in all_data_dirs:
        avg_rmsf_subdir = os.path.join(data_dir, "average_RMSF")
        if not os.path.isdir(avg_rmsf_subdir):
            # Some domains might not have "average_RMSF" if they failed extraction
            if args.verbose:
                print(f"[WARNING] No 'average_RMSF' directory in {data_dir}. Skipping.")
            continue
        
        # Find files that end with "temperature_413_average_rmsf.csv"
        pattern_413 = os.path.join(avg_rmsf_subdir, "*temperature_413_average_rmsf.csv")
        csv_files_413 = glob(pattern_413)
        
        for csv_file in csv_files_413:
            filename = os.path.basename(csv_file)
            dest_path = os.path.join(output_dir, filename)
            shutil.copy2(csv_file, dest_path)
            copied_files_count += 1
            if args.verbose:
                print(f"[INFO] Copied {csv_file} -> {dest_path}")

    if args.verbose:
        print(f"[INFO] Completed. {copied_files_count} CSV file(s) copied to {output_dir}.")

if __name__ == "__main__":
    main()