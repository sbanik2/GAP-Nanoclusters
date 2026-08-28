# Tersoff-HyBOP Potentials

This folder contains the **Tersoff-HyBOP** bond-order potential files used as a comparison baseline against the GAP models in this repository (see [Manna *et al.* (2022)](https://www.nature.com/articles/s41467-021-27849-6) and the associated [Nature Communications paper](https://doi.org/10.1038/s41467-021-27885-0)).

## Model

The Tersoff-HyBOP model consists of **two components**:

* A **Tersoff-style (BOP-like)** three-body interaction potential
* A **Lennard-Jones (LJ) scaling** component for long-range dispersion interactions

All model parameters are stored in individual `.json` files per element, in this folder. Each file includes:

* A `tersoff` section with named Tersoff parameters
* A `scaling` section with `epsilon`, `sigma`, `k1`, `k2`, and `RcLR` used in `lj/cut/scaling`

## Converting to LAMMPS input

To convert these `.json` files into **LAMMPS-ready inputs**, use the Jupyter notebook:

> 📓 `tersoff_lammps_template.ipynb`

This notebook automates:

* Generating the `.tersoff` file used with `pair_style tersoff`
* Writing a matching LAMMPS input script (`.in`) that uses `hybrid/overlay` with `tersoff` and `lj/cut/scaling` styles

## LAMMPS setup example

```lammps
pair_style hybrid/overlay tersoff lj/cut/scaling 14
pair_coeff * * tersoff Ag.tersoff Ag
pair_coeff 1 1 lj/cut/scaling epsilon sigma k1 k2 0 RcLR
```

## Citation

If you use the Tersoff-HyBOP models, please cite:

```bibtex
@article{manna2022learning,
  title     = {Learning in continuous action space for developing high dimensional potential energy models},
  author    = {Manna, Sukriti and Loeffler, Troy D and Batra, Rohit and Banik, Suvo and Chan, Henry and Varughese, Bilvin and Sasikumar, Kiran and Sternberg, Michael and Peterka, Tom and Cherukara, Mathew J and others},
  journal   = {Nature Communications},
  volume    = {13},
  number    = {1},
  pages     = {368},
  year      = {2022},
  publisher = {Nature Publishing Group UK London},
  doi       = {10.1038/s41467-021-27885-0}
}
```
