import os
import matplotlib
import h5py as h5
import numpy as np
import pandas as pd
from glob import glob
from os.path import join as opj
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Specify the case of interest
domain_id = '1r9lA02'
temperature = '348'
repl = '0'

# Define main function - This is what is called when python is called. 
def main():
    # Load the h5 file
    # data_dir = './'
    h5_file_path = opj(f"../mdcath_dataset_{domain_id}.h5") 
    
    # Check if the file exists
    if not os.path.exists(h5_file_path):
        print(f"Error: File '{h5_file_path}' not found.")
        return
    
    try:
        with h5.File(h5_file_path, 'r') as f:
            print('Successfully loaded the h5 file')
            print(f'dom: {list(f.keys())[0]}')
            print(f'layout: {f.attrs["layout"]}')
            print(f'molDatasets --> {list(f[domain_id].keys())}')
            print(f'molAttrs --> {list(f[domain_id].attrs.keys())}')
            print(domain_id)
            print(f"numChains --> {f[domain_id].attrs['numChains']}")
            print(f"numProteinAtoms --> {f[domain_id].attrs['numProteinAtoms']}")
            print(f"numResidues --> {f[domain_id].attrs['numResidues']}")
            # print(f"RESID ---- \n\n{f[domain_id]['resid'][()].decode('utf-8')}")

            
            print(f"z.shape --> {f[domain_id]['z'].shape}")
            print(f"z --> {f[domain_id]['z'][:10]}")

            # # Get CA idxs from the protein pdb file contained in the h5 file
            # pdbProteinAtoms = f[domain_id]['pdbProteinAtoms'][()].decode('utf-8').split('\n')[1:-3] # remove header and footer
            # atomtypes = [line.split()[2] for line in pdbProteinAtoms]
            # a_indices = np.where(np.array(atomtypes) == 'CA')[0]
            # print(f'Number of CA atoms: {len(ca_indices)}')
            # print(f"pdbProteinAtoms\n\n{f[domain_id]['pdbProteinAtoms'][()].decode('utf-8')}")
            print(f'available replicas ({temperature}K) --> {list(f[domain_id][temperature].keys())}')
            print(f'attrs ({temperature}K) --> {list(f[domain_id][temperature].attrs.keys())}')

            # for key, data in f[domain] AKA access domain and then [temperature] 
             # Access a specific nested group within the HDF5 file. So 

            for key, data in f[domain_id][temperature][str(repl)].items():
                print(f'prop {key} --> {data.shape}')
                for attr in data.attrs.keys():
                    print(f'{attr} --> {data.attrs[attr]}')
                print('' + 'TEST ')
                
            # Each replica contains the number of frames (numFrames) as attributes
            for replattr in f[domain_id][temperature][str(repl)].attrs.keys():
                print(f'{replattr} --> {f[domain_id][temperature][str(repl)].attrs[replattr]}')

            rmsd = f[domain_id][temperature][str(repl)]['rmsd'][:] # shape (numFrames,)
            rmsf = f[domain_id][temperature][str(repl)]['rmsf'][:] # shape (numResidues)
            gyration_radius = f[domain_id][temperature][str(repl)]['gyrationRadius'][:] # shape (numFrames,)
            print(f"rmsd.shape --> {rmsd.shape}")
            print(f"rmsf.shape --> {rmsf.shape}")
            print(f"gyration_radius.shape --> {gyration_radius.shape}")

            coords = f[domain_id][temperature][str(1)]["coords"]#[:,ca_indices,:]
            forces = f[domain_id][temperature][str(1)]["forces"]#[:,ca_indices,:]
            print(f'coords --> {coords.shape}, units: {f[domain_id][temperature]["0"]["coords"].attrs["unit"]}')
            print(f'forces --> {forces.shape}, units: {f[domain_id][temperature]["0"]["forces"].attrs["unit"]}')




            
    except Exception as e:
        print(f"An error occurred while reading the h5 file: {e}")

# Entry point
if __name__ == "__main__":
    main()