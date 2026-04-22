# Soil Model for Enhanced Weathering (SMEW)

This folder contains all material related to the Soil Model for Enhanced Weathering (SMEW) published in Bertagni et al., 2025, JAMES
The material includes the numerical codes for the model, the experimental data used for the analyses, and the Jupyter notebooks for the model-experiment comparison.

Article: [https://doi.org/10.1029/2024MS004224](https://doi.org/10.1029/2024MS004224)  
Zenodo repository: [https://doi.org/10.5281/zenodo.14356660](https://doi.org/10.5281/zenodo.14356660)  
GitHub repository: [https://github.com/MatteoBertagni/SMEW](https://github.com/MatteoBertagni/SMEW)  
PyPi package: [https://pypi.org/project/smew](https://pypi.org/project/smew)

You can install and use SMEW simply by pulling the Python package as follows:
```bash
pip install smew
```
Be aware that this doesn't install the cminpack library which is required for the fast version of SMEW (see below for 
the full installation procedure). 

## Installation

The fast version of SMEW relies on cminpack, a C-library for solving nonlinear equations. 
Because of this compiled dependency, the installation process depends on your system setup.
If you choose to not install the cminpack library, SMEW will fall back to the scipy solver (no additional configuration 
is required to use the Scipy fallback).
If cminpack is not installed, SMEW automatically falls back to the scipy solver.
While scipy also uses minpack (written in Fortran), this implementation is slower (~ x2) due to Python function call 
overhead.
In contrast, the cminpack library is called directly from Cython, significantly reducing that overhead.

### With Conda (no root privilege required)

Create and activate the Conda environment for cminpack:

```bash
conda create -n smew_env python cminpack -c conda-forge
conda activate smew_env
```

Install SMEW:

```bash
pip install smew
```

### With your system package manager (root privilege required)

If you are installing smew on a personal machine or a VM where you have root access, you can install cminpack 
system-wide and use a standard Python virtual environment.

For Ubuntu:

```bash
sudo apt-get update
sudo apt-get install libcminpack-dev
```

If you don't already have a virtual environment create one and activate it:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install SMEW:

```bash
pip install smew
```


## Folders

- `pyEW`: contains the python codes for the SMEW numerical model
- `Exp_data`: contains the experimental data obtained from the various publications through a web plot digitizer (https://apps.automeris.io/wpd/)
- `pyeto`: contains the Python codes to estimate the potential evapotranspiration (Mark Richards, https://pyeto.readthedocs.io/en/latest/index.html)


## Jupyter notebooks

- `Example`: provides an example of simulation for an EW application 
- `Vials_Dietzen`: model-experiment comparisons with the experiments by Dietzen et al. (2018)
- `Bottles_tePas`: model-experiment comparisons with the experiments by tePas et al. (2023)
- `Mesocosm_Amann`: model-experiment comparisons with the experiments by Amann et al. (2020)
- `Mesocosm_Kelland`: model-experiment comparisons with the experiments by Kelland et al. (2020)


## Instructions

1. Download or pull the whole repository into a selected working directory.
2. Run the Juptyer notebook 'Example' to verify that the model components (pyEW) are correctly used within the notebooks.
3. Change the parameters in the file 'Example' to run specific simulations for different scenarios.
4. For the Jupyter notebooks of the model-experiment comparison, define the selected base directory in the first cell of each notebook.

## License

This project is licensed under the GNU Affero General Public License v3.0 (AGPL-3.0).

See the `license` file for details.

## Contact

You can contact me at @MatteoBertagni ([matteo.bertagni@polito.it](mailto:matteo.bertagni@polito.it)) for more information about the research.