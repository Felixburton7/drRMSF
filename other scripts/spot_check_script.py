#!/usr/bin/env python3

import h5py
import numpy as np

def main():
    # Define input parameters
    h5_file = "mdcath_dataset_1r9lA02.h5"  # Ensure this is the correct path
    domain_id = "1r9lA02"
    temperature = "379"  # Adjust as needed
    replica = "0"        # Adjust as needed

    # Open the HDF5 file
    try:
        with h5py.File(h5_file, "r") as f:
            # Navigate to the domain group
            if domain_id not in f:
                raise KeyError(f"Domain '{domain_id}' not found in the file.")

            domain_group = f[domain_id]

            # Navigate to the specified temperature and replica
            if temperature not in domain_group:
                raise KeyError(f"Temperature group '{temperature}' not found in domain '{domain_id}'.")

            temp_group = domain_group[temperature]

            if replica not in temp_group:
                raise KeyError(f"Replica '{replica}' not found in temperature group '{temperature}'.")

            replica_group = temp_group[replica]

            # Extract RMSF data
            if "rmsf" not in replica_group:
                raise KeyError(f"'rmsf' dataset not found in replica '{replica}'.")

            rmsf_data = replica_group["rmsf"][:]
            print(f"RMSF Data for Domain: {domain_id}, Temperature: {temperature}K, Replica: {replica}")
            print(rmsf_data)

    except FileNotFoundError:
        print(f"Error: The file '{h5_file}' was not found.")
    except KeyError as e:
        print(f"KeyError: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()
