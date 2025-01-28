#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Your Name
Date: 27 January 2025
Environment: felix_env (conda)

Description:
    This script extends the simple histogram of domain-level mean RMSF values
    by producing an interactive Plotly histogram. For each temperature, a subplot
    (facet) displays a histogram of domain-level mean RMSF values on the x-axis
    and the count on the y-axis. When you hover over a bin, you see:
      - The range of mean RMSF for that bin
      - How many domain entries are in that bin
      - Which domains are included (domain_id)
      - The mean RMSF values (if multiple entries fall in the same bin)
    Additionally, we annotate each subplot with the total count of domain entries
    for that temperature. A footer annotation summarises the overall data set.

    Includes extra logging to confirm whether the script is working correctly
    and to show precisely what files have been parsed and where the output is saved.

Steps:
    1) Parallel parse each domain directory, ignoring replicates. We only look at:
       /home/s_felix/RMSF_data/rmsf_<domain>_data/average_RMSF/
       CSVs named <domainid>_temperature_<TEMP>_average_rmsf.csv
    2) Compute a single mean RMSF for each domain & temperature (across all residues).
    3) Build a Plotly FacetGrid histogram of mean RMSF (x) vs count (y), one subplot
       per temperature. Provide rich hover data and custom annotations.

Output:
    /home/felix/visualize/RMSF_histogram/rmsf_distribution_interactive.html
    (an interactive HTML file containing the dynamic plot)
"""

import os
import sys
import glob
import logging
import multiprocessing

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

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
    Parse the average_RMSF directory in a single domain folder, gathering
    domain-level mean RMSF per temperature.

    Returns a list of dicts with keys: domain, temperature, mean_rmsf
    """
    results = []

    domain_name = os.path.basename(domain_dir)
    logging.info(f"Starting parse for domain directory: {domain_name}")

    # Attempt to parse an ID from that name
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
            left_part, right_part = filename.split("_temperature_")
            file_domain_id = left_part  # e.g. "2eefA01"
            temp_part = right_part.split("_average_rmsf.csv")[0]
            temperature_str = temp_part.strip()

            # Read CSV
            df = pd.read_csv(csv_path)
            expected_col = f"rmsf_{temperature_str}"
            if expected_col not in df.columns:
                logging.warning(
                    f"Column {expected_col} not found in {csv_path}. "
                    f"Skipping this file."
                )
                continue

            # Convert RMSF column to numeric and compute mean
            rmsf_values = pd.to_numeric(df[expected_col], errors='coerce').dropna()
            if rmsf_values.empty:
                logging.warning(f"No valid RMSF values found in {filename}.")
                continue

            mean_rmsf = rmsf_values.mean()

            # Convert temperature to float
            try:
                temperature_val = float(temperature_str)
            except ValueError:
                logging.warning(f"Temperature '{temperature_str}' not numeric in {filename}. Skipping.")
                continue

            results.append({
                "domain": file_domain_id,
                "temperature": temperature_val,
                "mean_rmsf": mean_rmsf
            })

        except Exception as e:
            logging.error(f"Error parsing {filename} in {domain_dir}: {e}")

    logging.info(f"Finished parse for {domain_name}; collected {len(results)} entries.")
    return results

def main():
    """
    Main workflow:
    1) Identify domain directories and parse them in parallel to get domain-level RMSF.
    2) Create an interactive Plotly histogram with subplots (facets) by temperature.
    3) Annotate each subplot with total domain count, and add a footer annotation.
    4) Save output to an HTML file, logging each key step.
    """
    logging.info(f"Checking if output directory exists: {OUTPUT_DIR}")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    logging.info(f"Output directory ready at: {OUTPUT_DIR}")

    # Identify domain directories
    logging.info(f"Searching for domain directories in {BASE_DIR}")
    domain_dirs = glob.glob(os.path.join(BASE_DIR, "rmsf_*_data"))
    if not domain_dirs:
        logging.warning(f"No domain directories found in {BASE_DIR}. Exiting.")
        sys.exit(0)
    logging.info(f"Found {len(domain_dirs)} domain directories to parse.")

    # Parallel parse each domain's average RMSF data
    logging.info("Beginning parallel parsing of domain directories...")
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    results = pool.map(parse_domain_average_rmsf, domain_dirs)
    pool.close()
    pool.join()

    # Flatten
    all_data = []
    for r in results:
        all_data.extend(r)

    if not all_data:
        logging.warning("No average RMSF data gathered. Exiting without plotting.")
        sys.exit(0)
    logging.info(f"Collected {len(all_data)} total domain-temperature entries.")

    df = pd.DataFrame(all_data)  # columns: domain, temperature, mean_rmsf
    df.sort_values(by="temperature", inplace=True)

    # Count how many data points (domains) per temperature
    temp_counts = df.groupby("temperature")["domain"].count().reset_index(name="count")

    # Generate a facet histogram
    logging.info("Generating interactive Plotly histogram.")
    fig = px.histogram(
        df,
        x="mean_rmsf",
        facet_col="temperature",
        facet_col_wrap=3,
        nbins=15,
        hover_data=["domain", "mean_rmsf"],
        labels={
            "mean_rmsf": "Mean RMSF",
            "count": "Count",
            "temperature": "Temperature"
        }
    )

    fig.update_layout(
        title_text="Interactive Histogram of Domain-Level Mean RMSF by Temperature",
        showlegend=False,
        margin=dict(t=100, r=100, l=60, b=150)  # Increased bottom margin
    )

    # Add annotation in the top-right corner for each subplot
    unique_temps = sorted(df["temperature"].unique())
    for i, temp in enumerate(unique_temps):
        row = i // 3 + 1
        col = i % 3 + 1
        count_val = temp_counts.loc[temp_counts["temperature"] == temp, "count"].values[0]
        annotation_text = f"Total count: {count_val}"

        # Calculate relative position for each facet
        fig.add_annotation(
            text=annotation_text,
            x=0.95,
            y=0.95,
            xref=f"x{row}{col} domain",
            yref=f"y{row}{col} domain",
            showarrow=False,
            font=dict(size=12, color="black"),
            align="right",
            bordercolor="black",
            borderwidth=1,
            borderpad=2,
            bgcolor="white",
            opacity=0.8
        )

    # Add footer annotation across the entire figure with multi-line text
    total_domains = df["domain"].nunique()
    total_entries = len(df)
    footer_text = (
        f"Extracted data from {total_domains} unique domain(s),\n"
        f"totalling {total_entries} domain-temperature entries.\n"
        "Hover over the bars to see which domain(s) contributed to each bin."
    )

    fig.add_annotation(
        text=footer_text,
        x=0.5,
        y=-0.15,  # Adjusted y-position to move footer below the plot
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(size=12, color="black"),
        align="center"
    )

    # Tweak axis labels
    fig.update_xaxes(title_text="Mean RMSF")
    fig.update_yaxes(title_text="Count")

    # Save to an interactive HTML file
    out_html = os.path.join(OUTPUT_DIR, "rmsf_distribution_interactive.html")
    logging.info(f"Saving interactive plot to: {out_html}")
    try:
        fig.write_html(out_html, include_plotlyjs="cdn")
        logging.info("Interactive histogram successfully saved.")
    except Exception as e:
        logging.error(f"Failed to save interactive histogram: {e}")

if __name__ == "__main__":
    main()
