#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from quippy.potential import Potential
import numpy as np
from ase.io import read, write
from scipy.stats import spearmanr
from multiprocessing import Pool
from tqdm import tqdm
from utilities import CCC
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile
import subprocess
import os
from ase import Atoms
import json
from ase.io import read


# In[ ]:


def _load_ground_state_structures(ground_state_xyz_file):
    """
    Reads ground-state structures from an XYZ file and extracts DFT energies.
    Assumes that each structure has a unique tag (e.g., 'fcc', 'bcc') in the comment field.
    Returns:
    - A dictionary where keys are structure names (e.g., 'fcc') and values are {'atoms': ASE_Atoms, 'energy': DFT energy}.
    """
    ground_state_structures = read(ground_state_xyz_file, index=":")
    gs_atoms_dict = {}

    for atoms in ground_state_structures:
        structure_tag = atoms.info.get("config_type", None)  # Extract tag (e.g., 'fcc', 'bcc')
        dft_energy = atoms.get_potential_energy()/len(atoms)

        if structure_tag:
            gs_atoms_dict[structure_tag] = {"atoms": atoms, "energy": dft_energy}

    return gs_atoms_dict



def relax_structure( atoms, calculator, fmax=0.01, relax_cell=False):
    """
    Relax the structure using BFGS optimizer.
    
    Parameters:
    - atoms: ASE Atoms object
    - calculator: GAP potential calculator
    - fmax: Maximum force criterion for relaxation
    - relax_cell: If True, relax both atomic positions and unit cell. 
                  If False, only relax atomic positions (keep cell fixed).
    """
    atoms.set_calculator(calculator)

    if relax_cell:
        # Full relaxation (Atoms + Cell)
        ucf = UnitCellFilter(atoms)  # This allows cell to relax
        opt = BFGS(ucf, logfile='relax.log')
    else:
        # Atomic relaxation only (fix cell)
        opt = BFGS(atoms, logfile='relax.log')
    opt.run(fmax=fmax,steps=300)
    

def _get_lattice_parameter(atoms):
    cell = atoms.get_cell()
    
    # Lattice parameters
    a = np.linalg.norm(cell[0])  
    b = np.linalg.norm(cell[1])
    c = np.linalg.norm(cell[2])
    
    # Lattice angles (in degrees)
    alpha = np.degrees(np.arccos(np.dot(cell[1], cell[2]) / (b * c)))  # Angle between b and c
    beta = np.degrees(np.arccos(np.dot(cell[0], cell[2]) / (a * c)))   # Angle between a and c
    gamma = np.degrees(np.arccos(np.dot(cell[0], cell[1]) / (a * b)))  # Angle between a and b
    return [a, b, c, alpha, beta, gamma]
    

def compute_properties(potential_file,ground_state_xyz_file):
    """
    Computes relaxed lattice parameters and cohesive energies for ground-state structures
    using the specified GAP potential file.
    
    Parameters:
    - potential_file (str): Path to the GAP potential file.
    """
    results = {}  # Stores computed ground-state properties

    for structure, data in _load_ground_state_structures(ground_state_xyz_file).items():
        atoms = data['atoms'].copy()
        dft_energy_per_atom = data['energy']

        lattice_param_dft = _get_lattice_parameter(atoms)

        calculator = Potential(param_filename=potential_file)
        # Relax the structure (cell relaxation enabled)
        relax_structure(atoms, calculator, relax_cell=True)

        # Compute lattice parameter
        lattice_param_gap = _get_lattice_parameter(atoms)

        # Compute cohesive energy per atom
        atoms.set_calculator(calculator)
        gap_energy_per_atom = atoms.get_potential_energy() / len(atoms)

        # Store results
        results[structure] = {
            'lattice_param_gap': lattice_param_gap,
            'lattice_param_dft': lattice_param_dft,
            'gap_energy': gap_energy_per_atom,
            'dft_energy': dft_energy_per_atom
        }
    return results


# In[ ]:


data = compute_properties('gap_model.xml',
                          'ground_state.xyz')


# In[ ]:


with open('lattice_GS.json','w') as outfile:
    json.dump(data,outfile)


# In[ ]:





# In[ ]:





# In[ ]:


def ef_train_set(potential_file,train_data_file, modalities, save_e=None, save_f=None):
    
    E_modalities = {}
    F_modalities = {}

    # Read test structures
    structures = read(train_data_file, index=":")

    # Initialize modality dictionaries
    for modality in modalities:
        E_modalities[modality] = {"dft_energy": [], "gap_energy": []}
        F_modalities[modality] = {"dft_force": [], "gap_force": []}

    # Run QUIP to compute energies and forces
    temp_output_name = train_data_file.replace(".xyz", "_quip.xyz")
    quip_command = f"quip E=T F=T atoms_filename={train_data_file} param_filename={potential_file} | grep AT | sed 's/AT//' > {temp_output_name}"
    subprocess.run(quip_command, shell=True, check=True)

    # Read predicted structures from QUIP output
    atoms_list = read(temp_output_name, index=":")

    # Extract energies and forces for each structure
    for atomsi, atomsp in zip(structures, atoms_list):
        modality = atomsi.info['config_type']

        # Normalize by number of atoms
        E_modalities[modality]["dft_energy"].append(atomsi.get_potential_energy() / len(atomsi))
        E_modalities[modality]["gap_energy"].append(atomsp.get_potential_energy() / len(atomsp))

        # Collect forces
        F_modalities[modality]["dft_force"] += atomsi.get_forces().flatten().tolist()
        F_modalities[modality]["gap_force"] += atomsp.get_array('force').flatten().tolist()

    # Remove temporary output file
    os.system(f'rm {temp_output_name} *.idx')

    # Save energy and force data if specified
    if save_e:
        with open(save_e, 'w') as outfile:
            json.dump(E_modalities, outfile)

    if save_f:
        with open(save_f, 'w') as outfile:
            json.dump(F_modalities, outfile)



# In[ ]:


modalities = ['cluster', 'bulk']

ef_train_set('gap_model.xml',
             'train.xyz', 
             modalities, 
             save_e='E_train.json', 
             save_f='F_train.json')


# In[ ]:




