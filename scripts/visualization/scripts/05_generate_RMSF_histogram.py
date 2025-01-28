#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Your Name
Date: 27 January 2025
Environment: felix_env (conda)

Description:
    This script pools across multiple CPUs to parse only the average RMSF
    data (ignoring replicates), computing a single mean RMSF per domain
    for each temperature. It then plots a histogram for each temperature
    side-by-side using Seaborn's FacetGrid.

    British English is used throughout.

Steps:
    1. Identify domain directories of the form: rmsf_<domainname>_data.
    2. Within each domain directory, enter the subdirectory "average_RMSF".
    3. For each CSV matching:
         <domainid>_temperature_<TEMP>_average_rmsf.csv
       - Parse the column "rmsf_<TEMP>".
       - Compute a single mean RMSF across all residues in that file.
       - Append (domain_id, temperature, domain_mean_rmsf) to a list.
    4. Combine results across all domains.
    5. Plot distribution histograms for each temperature using FacetGrid.

Output:
    - Figure saved to /home/s_felix/visualize/RMSF_histogram/rmsf_distribution.png

Usage:
    1. conda activate felix_env
    2. conda install pandas matplotlib seaborn -y
    3. Place this script in /home/s_felix/drWiggle/scripts/
    4. python /home/s_felix/drWiggle/scripts/05_generate_RMSF_histogram.py
"""

import os
import sys
import glob
import logging
import multiprocessing

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Global paths
BASE_DIR = "/home/s_felix/RMSF_data"
OUTPUT_DIR = "/home/s_felix/visualize/RMSF_histogram"

# Logging config
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def parse_domain_average_rmsf(domain_dir):
    """
    Parse the average_RMSF directory in a single domain folder and gather
    domain-level mean RMSF per temperature.

    Returns a list of dicts with keys: domain, temperature, mean_rmsf
    """
    results = []

    domain_name = os.path.basename(domain_dir)  # e.g. "rmsf_2eefA01_data"
    # The domain_id might be "2eefA01" if we strip the "rmsf_" prefix and "_data" suffix
    # We'll attempt to parse that for clarity, but we can simply store domain_name.
    # However, let's try to do it neatly:
    if domain_name.startswith("rmsf_") and domain_name.endswith("_data"):
        domain_id = domain_name[len("rmsf_"):-len("_data")]
    else:
        domain_id = domain_name

    avg_dir = os.path.join(domain_dir, "average_RMSF")
    if not os.path.isdir(avg_dir):
        logging.warning(f"No average_RMSF directory found for {domain_dir}. Skipping.")
        return results

    # Find CSV files of the form <domainid>_temperature_<TEMP>_average_rmsf.csv
    avg_files = glob.glob(os.path.join(avg_dir, "*_temperature_*_average_rmsf.csv"))
    if not avg_files:
        logging.warning(f"No CSV files found in {avg_dir}. Skipping.")
        return results

    for csv_path in avg_files:
        filename = os.path.basename(csv_path)
        try:
            # Example filename: 2eefA01_temperature_379_average_rmsf.csv
            left_part, right_part = filename.split("_temperature_")  # ["2eefA01", "379_average_rmsf.csv"]
            file_domain_id = left_part  # e.g. "2eefA01"
            temp_part = right_part.split("_average_rmsf.csv")[0]  # e.g. "379"
            temperature_str = temp_part.strip()

            # Read CSV
            df = pd.read_csv(csv_path)
            # Expected column: f"rmsf_{temperature}"
            expected_col = f"rmsf_{temperature_str}"
            if expected_col not in df.columns:
                logging.warning(
                    f"Column {expected_col} not found in {csv_path}. Skipping this file."
                )
                continue

            # Compute domain-level mean for this temperature
            # Convert the column explicitly to float
            rmsf_values = pd.to_numeric(df[expected_col], errors='coerce').dropna()
            if rmsf_values.empty:
                continue

            mean_rmsf = rmsf_values.mean()

            # Convert temperature to float if possible
            try:
                temperature_val = float(temperature_str)
            except ValueError:
                # if not convertible, skip
                logging.warning(f"Temperature '{temperature_str}' not numeric in {filename}. Skipping.")
                continue

            results.append({
                "domain": file_domain_id,
                "temperature": temperature_val,
                "mean_rmsf": mean_rmsf
            })

        except Exception as e:
            logging.error(f"Error parsing {filename} in {domain_dir}: {e}")

    logging.info(f"Parsed {len(results)} average RMSF entries from {domain_dir}")
    return results

def main():
    """
    Main entry point: 
    1) Identify domain directories
    2) Parallel parse each domain's average_RMSF data
    3) Plot a histogram for each temperature side by side
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Locate domain directories
    domain_dirs = glob.glob(os.path.join(BASE_DIR, "rmsf_*_data"))
    if not domain_dirs:
        logging.warning(f"No domain directories found in {BASE_DIR}. Exiting.")
        sys.exit(0)

    # Use multiprocessing to parse each domain in parallel
    # Adjust processes to suit your system (e.g. processes=8)
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(parse_domain_average_rmsf, domain_dirs)
    pool.close()
    pool.join()

    # Flatten list of lists
    all_data = []
    for r in results:
        all_data.extend(r)

    if not all_data:
        logging.warning("No average RMSF data gathered. Exiting without plotting.")
        sys.exit(0)

    # Convert to DataFrame
    df = pd.DataFrame(all_data)  # columns: domain, temperature, mean_rmsf

    # Sort by temperature
    df.sort_values(by="temperature", inplace=True)

    # Create a FacetGrid histogram, one column per temperature
    # This shows the distribution of domain-level mean RMSF at each temperature
    g = sns.FacetGrid(df, col="temperature", col_wrap=3, sharex=False, sharey=False)
    g.map(sns.histplot, "mean_rmsf", kde=True, bins=100, color="steelblue")

    g.set_titles(col_template="Temperature: {col_name}")
    g.set_xlabels("Mean RMSF")
    g.set_ylabels("Count")

    plt.suptitle("Domain-Level Mean RMSF Distribution at Each Temperature", y=1.02)
    plt.tight_layout()

    # Save the figure
    out_path = os.path.join(OUTPUT_DIR, "rmsf_distribution.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()

    logging.info(f"Histogram saved to {out_path}")

if __name__ == "__main__":
    main()

# /// Split into uneven shaped bins, 
# /// So that each bin has the same number of entries 
# 
