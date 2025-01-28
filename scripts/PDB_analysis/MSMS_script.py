

import os
import subprocess
import pandas as pd
import os.path as p

# Function to classify residues into core and exterior based on SES scores
def findCoreExterior(pdbFile, msmsDir, pdbDf, proteinName, outDir):
    # Change the working directory so MSMS can find all the files it needs
    os.chdir(msmsDir)

    # Define paths to the required executables
    pdb2xyzrExe = "./pdb_to_xyzr"  # Converts PDB to MSMS .xyzr file
    msmsExe = "./msms.x86_64Linux2.2.6.1"  # Computes SES using MSMS

    # Convert the PDB file to an MSMS-compatible XYZR file
    xyzrFile = p.join(outDir, f'{proteinName}.xyzr')  # Define output .xyzr file path
    command = f"{pdb2xyzrExe} {pdbFile} > {xyzrFile}"  # Command to run pdb_to_xyzr
    subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # Execute command

    # Use MSMS to calculate the SES and create an area file
    areaOut = p.join(outDir, proteinName)  # Base name for area output
    command = f"{msmsExe} -if {xyzrFile} -af {areaOut}"  # Command to run MSMS
    subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # Execute command

    # Define the path to the generated .area file
    areaFile = p.join(outDir, f"{proteinName}.area")

    # Convert the .area file into a DataFrame containing SES scores
    areaDf = area2df(areaFile=areaFile)

    # Merge the SES scores with the original PDB DataFrame
    pdbDf = pd.concat([pdbDf, areaDf], axis=1)

    # Group by residue ID and calculate the average SES score for each residue
    meanSesPerResidue = pdbDf.groupby('RES_ID')['SES'].mean()

    # Identify residues with average SES > 1 (considered as exterior residues)
    exteriorResiduesIndex = meanSesPerResidue[meanSesPerResidue > 1].index

    # Split the DataFrame into exterior and core residues
    exteriorDf = pdbDf[pdbDf['RES_ID'].isin(exteriorResiduesIndex)]  # Exterior residues
    coreDf = pdbDf[~pdbDf['RES_ID'].isin(exteriorResiduesIndex)]  # Core residues

    # Clean up temporary files created during processing
    os.remove(xyzrFile)  # Delete the .xyzr file
    os.remove(areaFile)  # Delete the .area file

    # Return DataFrames for exterior and core residues
    return exteriorDf, coreDf


# Helper function to convert the .area file into a pandas DataFrame
def area2df(areaFile):
    ses = []  # List to store SES scores

    # Read the .area file line by line
    with open(areaFile, "r") as file:
        for line in file:
            if "Atom" in line:  # Skip header or irrelevant lines
                continue
            cols = line.split()  # Split line into columns
            ses.append(float(cols[1]))  # Extract SES score (second column)

    # Create a DataFrame from the SES scores
    data = {"SES": ses}
    pdbDf = pd.DataFrame(data)

    # Return the DataFrame
    return pdbDf