#!/usr/bin/env python3
"""
04_extract_RMSF_Master.py

Master script to parallelise 03_extract_RMSF.py across all .h5 files
in a specified directory. Each .h5 is handled in one job, which
replicates the logic of 01_extract_RMSF.py for all domains found
within that file.

Usage:
  python 04_extract_RMSF_Master.py \
      --h5_dir /mnt/datasets/MD_CATH/data \
      --output_dir /home/s_felix/RMSF_data \
      --verbose
"""

import os
import argparse
import subprocess
from multiprocessing import Pool, cpu_count

def parse_master_arguments():
    parser = argparse.ArgumentParser(
        description="Parallelise RMSF extraction for multiple .h5 files."
    )
    parser.add_argument("--h5_dir", required=True, help="Directory containing .h5 files.")
    parser.add_argument(
        "--output_dir",
        default="rmsf_output",
        help="Directory for storing RMSF outputs for all .h5 files."
    )
    parser.add_argument("--verbose", action="store_true", help="Print detailed logs.")
    return parser.parse_args()

def discover_h5_files(h5_dir):
    """
    Return a list of absolute paths to .h5 files in 'h5_dir'.
    """
    files = []
    for filename in sorted(os.listdir(h5_dir)):
        if filename.endswith(".h5"):
            files.append(os.path.join(h5_dir, filename))
    return files

def worker(task):
    """
    Worker function that calls 03_extract_RMSF.py on a single .h5 file.
    """
    h5_file, output_dir, verbose = task
    script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "03_extract_RMSF.py")

    cmd = [
        "python", script_path,
        "--h5_file", h5_file,
        "--output_dir", output_dir
    ]
    if verbose:
        cmd.append("--verbose")
        print(f"[DEBUG] Executing command: {' '.join(cmd)}")

    # Run the single extraction script
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[ERROR] Failed to process {h5_file}. Error:\n{result.stderr}")
    else:
        if verbose:
            print(f"[INFO] Finished processing {h5_file}.\n{result.stdout}")

def main():
    args = parse_master_arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    h5_files = discover_h5_files(args.h5_dir)
    if not h5_files:
        print("[INFO] No .h5 files found in the specified directory.")
        return

    tasks = [(h5_file, args.output_dir, args.verbose) for h5_file in h5_files]
    total_tasks = len(tasks)
    if args.verbose:
        print(f"[INFO] Found {total_tasks} .h5 file(s) to process.")

    # Number of parallel workers: min(#files, #CPUs)
    num_processes = min(cpu_count(), total_tasks)
    if args.verbose:
        print(f"[INFO] Starting pool with {num_processes} parallel process(es).")

    with Pool(processes=num_processes) as pool:
        pool.map(worker, tasks)

    if args.verbose:
        print("[INFO] All tasks have been processed.")

if __name__ == "__main__":
    main()
