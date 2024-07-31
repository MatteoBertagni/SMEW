# Soil Model for Enhanced Weathering (SMEW)
This folder contains all material related to the Soil Model for Enhanced Weathering (SMEW) published in Bertagni et al.
The material includes the numerical codes for the model, the experimental data used for the analyses, and the Jupyter notebooks for the model-experiment comparison.

# Folders

- 'pyEW': contains the python codes for the SMEW numerical model

- 'Exp_data': contains the experimental data obtained from the various publications through a web plot digitizer (https://apps.automeris.io/wpd/)

- 'pyeto': contains the Python codes to estimate the potential evapotranspiration (Mark Richards, https://pyeto.readthedocs.io/en/latest/index.html)

# Jupyter notebooks

- 'Example': provides an example of simulation for an EW application 

- 'Vials_Dietzen': model-experiment comparisons with the experiments by Dietzen et al. (2018)

- 'Bottles_tePas': model-experiment comparisons with the experiments by tePas et al. (2023)

- 'Mesocosm_Amann': model-experiment comparisons with the experiments by Amann et al. (2020)

- 'Mesocosm_Kelland': model-experiment comparisons with the experiments by Kelland et al. (2020)


# Instructions
1. Download or pull the whole repository into a selected working directory.
2. Run the Juptyer notebook 'Example' to verify that the model components (pyEW) are correctly used within the notebooks.
3. Change the parameters in the file 'Example' to run specific simulations for different scenarios.
4. For the Jupyter notebooks of the model-experiment comparison, define the selected base directory in the first cell of each notebook.

# Contact
You can contact me at @MatteoBertagni (matteo.bertagni@polito.it) for more information about the research.
