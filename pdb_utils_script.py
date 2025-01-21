# PDB Files are annoying and not compatible with pandas or other analytical tools. pdbUtils reads a PDB file and vonverts it into a pandas DataFrame, where each row correspnods to an atom, and each column representas and attribute. 
# Built in module for interacting with the operating systme. It provides tools for file and dir operations. 
#It is a module in python


import os
from pdbUtils import pdbUtils

# Input and output file paths
input_pdb_file = "example_outputs/output_pdbs/1r9lA02.pdb"
output_dir = "./"
#os.path submodule within os handles common operations on files. 'join' function combines multiple paths into on string
output_pdb_file = os.path.join(output_dir, "1r9lA02_cleaned.pdb")

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load PDB file into a DataFrame
pdb_df = pdbUtils.pdb2df(input_pdb_file)

# (Optional) Perform any modifications on the DataFrame here

# Write the DataFrame back to a PDB file
pdbUtils.df2pdb(pdb_df, output_pdb_file)

print(f"Processed PDB file saved to {output_pdb_file}")