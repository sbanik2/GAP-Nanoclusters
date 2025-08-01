#!/usr/bin/env python
# coding: utf-8

import numpy as np
import subprocess
import os
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from quippy.potential import Potential
from ase.io import read, write
from utilities import CCC
import json

class PropertyCalculator:
    def __init__(self, ground_state_xyz_file):
        self.ground_state_xyz_file = ground_state_xyz_file
        self.ground_state_atoms = self._load_ground_state_structures()

    def _load_ground_state_structures(self):
        ground_state_structures = read(self.ground_state_xyz_file, index=":")
        gs_atoms_dict = {}

        for atoms in ground_state_structures:
            structure_tag = atoms.info.get("config_type", None)
            dft_energy = atoms.get_potential_energy() / len(atoms)

            if structure_tag:
                gs_atoms_dict[structure_tag] = {"atoms": atoms, "energy": dft_energy}

        return gs_atoms_dict

    def relax_structure(self, atoms, calculator, fmax=0.01, relax_cell=False):
        atoms.set_calculator(calculator)

        if relax_cell:
            ucf = UnitCellFilter(atoms)
            opt = BFGS(ucf, logfile='relax.log')
        else:
            opt = BFGS(atoms, logfile='relax.log')

        opt.run(fmax=fmax, steps=300)

    def _get_lattice_parameter(self, atoms):
        cell = atoms.get_cell()
        return [np.linalg.norm(cell[i]) for i in range(3)]

    def compute_properties(self, potential_file):
        results = {}

        for structure, data in self.ground_state_atoms.items():
            atoms = data['atoms'].copy()
            dft_energy_per_atom = data['energy']

            lattice_param_dft = self._get_lattice_parameter(atoms)
            calculator = Potential(param_filename=potential_file)

            self.relax_structure(atoms, calculator, relax_cell=True)
            lattice_param_gap = self._get_lattice_parameter(atoms)

            atoms.set_calculator(calculator)
            gap_energy_per_atom = atoms.get_potential_energy() / len(atoms)

            results[structure] = {
                'lattice_param_gap': lattice_param_gap,
                'lattice_param_dft': lattice_param_dft,
                'gap_energy': gap_energy_per_atom,
                'dft_energy': dft_energy_per_atom
            }

        return results

    def _compute_energy_ccc(self, potential_file, xyz_file):
        structures = read(xyz_file, index=":")
        gap_energies, dft_energies = [], []

        for atoms in structures:
            dft_energy = atoms.get_potential_energy()
            gap_energies.append(dft_energy)

            calculator = Potential(param_filename=potential_file)
            self.relax_structure(atoms, calculator, relax_cell=False)
            atoms.set_calculator(calculator)
            gap_energy = atoms.get_potential_energy()
            gap_energies.append(gap_energy)

        results_storage = {'gap_energy': gap_energies, 'dft_energy': dft_energies}
        return results_storage

    def compute_lattice_errors(self, potential_file, save_file=None):
        e_dft, e_gap, del_latt = [], [], []
        results = self.compute_properties(potential_file)

        if save_file:
            with open(save_file, 'w') as outfile:
                json.dump(results, outfile)

        for structure, data in results.items():
            gap_lattice = data['lattice_param_gap']
            e_gap.append(data['gap_energy'])
            dft_lattice = data['lattice_param_dft']
            e_dft.append(data['dft_energy'])
            del_latt.append(np.mean(np.abs(np.array(gap_lattice) - np.array(dft_lattice))))

        ordering_cohe = 1 if list(range(len(e_dft))) == list(range(len(e_gap))) else 0
        weights_normalized = np.abs(e_dft) / np.sum(np.abs(e_dft))

        e_dft, e_gap = np.array(e_dft) * weights_normalized, np.array(e_gap) * weights_normalized
        delE_cohe = 1 - CCC(e_dft - e_dft.min(), e_gap - e_gap.min())
        delLatt = np.sum(np.array(del_latt) * weights_normalized)

        return delE_cohe, 1 - ordering_cohe, delLatt

    def ef_test_set(self, potential_file, test_data_file, modalities, save_e=None, save_f=None):
        E_modalities = {modality: {"dft_energy": [], "gap_energy": []} for modality in modalities}
        F_modalities = {modality: {"dft_force": [], "gap_force": []} for modality in modalities}

        structures = read(test_data_file, index=":")
        temp_output_name = test_data_file.replace(".xyz", "_quip.xyz")
        quip_command = f"quip E=T F=T atoms_filename={test_data_file} param_filename={potential_file} | grep AT | sed 's/AT//' > {temp_output_name}"
        subprocess.run(quip_command, shell=True, check=True)

        atoms_list = read(temp_output_name, index=":")

        for atomsi, atomsp in zip(structures, atoms_list):
            modality = atomsi.info['config_type']
            E_modalities[modality]["dft_energy"].append(atomsi.get_potential_energy() / len(atomsi))
            E_modalities[modality]["gap_energy"].append(atomsp.get_potential_energy() / len(atomsp))

            F_modalities[modality]["dft_force"] += atomsi.get_forces().flatten().tolist()
            F_modalities[modality]["gap_force"] += atomsp.get_array('force').flatten().tolist()

        os.system(f'rm {temp_output_name} *.idx')

        if save_e:
            with open(save_e, 'w') as outfile:
                json.dump(E_modalities, outfile)

        if save_f:
            with open(save_f, 'w') as outfile:
                json.dump(F_modalities, outfile)

        ccc_errors = {"energy": {}, "force": {}}
        for modality, values in E_modalities.items():
            if values["dft_energy"] and values["gap_energy"]:
                ccc_errors["energy"][modality] = CCC(values["dft_energy"], values["gap_energy"])

        for modality, values in F_modalities.items():
            if values["dft_force"] and values["gap_force"]:
                ccc_errors["force"][modality] = CCC(values["dft_force"], values["gap_force"])

        return ccc_errors
