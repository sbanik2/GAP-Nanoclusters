# Gaussian Approximation Potentials for Elemental Nanoclusters

<p align="center">
  <img src="./index.png" alt="GAP Nanocluster Overview" width="800">
</p>

This repository contains **Gaussian Approximation Potential (GAP)** models for **54 elemental nanocluster systems** across the periodic table. For each element, the following are provided:

* GAP potential model files
* Tersoff potential files for comparison
* Performance metrics and dataset statistics
* Energy and force parity plots
* Results of dynamic stability tests
* Normal mode analysis
* Lattice parameter and cohesive energy ordering of ground state polytypes
* Comparison of energy and force errors with Tersoff models [Manna *et al.* (2022)](https://www.nature.com/articles/s41467-021-27849-6)
* Lattice parameter and cohesive energy prediction comparisons

This resource aims to support the development, benchmarking, and evaluation of machine-learning interatomic potentials against traditional models in low-dimensional nanocluster systems.

---

## Dataset

The training and test datasets are made available **upon request**.

The **ground state polymorphs** of the elemental nanoclusters and a portion of the **validation dataset** are publicly accessible via the [Quantum Cluster Database (QCD)](https://muellergroup.jhu.edu/qcd/).

---


## Usage

### Relaxing Structures with GAP

This example demonstrates how to relax a nanocluster structure using a trained **Gaussian Approximation Potential (GAP)** model via the **[QUIP](https://libatoms.github.io/QUIP/)** interface and **[ASE (Atomic Simulation Environment)](https://wiki.fysik.dtu.dk/ase/)**.

The script below performs a geometry optimization (relaxation) on a single `.xyz` structure using the GAP model:

```python
from ase.io import read, write
from ase.optimize import BFGS
from quippy.potential import Potential

# --- Settings ---
element = "Ag"  # Replace with your element symbol
structure_path = f"{element}_structure.xyz"  # Input structure file
gap_file = f"./gap_files/{element}.xml"      # Path to the GAP model file

# --- Read structure ---
atoms = read(structure_path)

# --- Assign GAP calculator ---
gap_calc = Potential(param_filename=gap_file)
atoms.calc = gap_calc

# --- Relaxation using BFGS optimizer ---
optimizer = BFGS(atoms, logfile=None)
optimizer.run(fmax=0.001, steps=500)

# --- Output results ---
relaxed_energy = atoms.get_potential_energy()
print(f"Relaxed GAP energy: {relaxed_energy:.6f} eV")

### 📌 Requirements

* `ase`
* `quippy` (install via QUIP or conda build of `libatoms`)


---


