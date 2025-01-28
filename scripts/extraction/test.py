

from pdbUtils import pdbUtils
from os import path as p; 
import os
import re
import h5py
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import subprocess
from multiprocessing import Pool, cpu_count
from typing import Tuple



class FilePath:
    pass
class DirPath:
    pass

def findCoreExterior(pdbFile: FilePath, msmsDir: DirPath, pdbDf: pd.DataFrame, proteinName: str, outDir: DirPath) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # change working directory so MSMS can find all the files it needs
    os.chdir(msmsDir)
    # find executables
    pdb2xyzrExe = "./pdb_to_xyzr"
    msmsExe = "./msms.x86_64Linux2.2.6.1"
    # convert pdb file to MSMS xyzr file
    xyzrFile = p.join(outDir, f'{proteinName}.xyzr')
    command = f"{pdb2xyzrExe} {pdbFile} > {xyzrFile}"
    subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # use MSMS to create an area file
    areaOut = p.join(outDir,proteinName)
    command = f"{msmsExe} -if {xyzrFile} -af {areaOut}"
    subprocess.run(command, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    areaFile=p.join(outDir,f"{proteinName}.area")
    # convert area file to dataframe, merge with main pdb dataframe
    areaDf = area2df(areaFile=areaFile)
    pdbDf = pd.concat([pdbDf,areaDf],axis=1)
 
    # Group by residue and calculate the average SES score
    meanSesPerResidue = pdbDf.groupby('RES_ID')['SES'].mean()
 
    # Get residue sequences with average SES > 1
    exteriorResiduesIndex = meanSesPerResidue[meanSesPerResidue > 1].index
 
    # Split the DataFrame based on average SES > 1
    exteriorDf = pdbDf[pdbDf['RES_ID'].isin(exteriorResiduesIndex)]
    coreDf = pdbDf[~pdbDf['RES_ID'].isin(exteriorResiduesIndex)]
 
    # clean up
    os.remove(xyzrFile)
    os.remove(areaFile)
 
    return exteriorDf, coreDf
 
def area2df(areaFile):
    ses =[]
    with open(areaFile,"r") as file:
        for line in file:
            if "Atom" in line:
                continue
            cols = line.split()
            ses.append(float(cols[1]))
    data = {"SES":ses}
    pdbDf = pd.DataFrame(data)
    return pdbDf