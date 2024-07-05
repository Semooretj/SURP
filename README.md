# SURP
This repository contains code and figures produced during the SURP 2024 program.
Supervised by Ioana Zelko and Joshua Speagle, the project aims to create improved versions of 3D Interstellar Medium Dust Temperature Maps from Zelko et al. 2022: https://arxiv.org/abs/2211.07667

Uses data files from dustfizz3d repository: https://github.com/ioanazelko/dustfizz3D/tree/main
## Overview
- Main project files for superpixel rotations and map shifts are located on the home page
- Files for CTA200 project contained in the *'project'* folder
- Extra files for both are located in the *'others'* folder
## Python files 
- **comparison_functions.py**

  To be used with the *compare_plank_rotations* notebook.
  
  To understand the effect of two rotation methods, this file has functions for generating residual, ratio and histogram plots of a rotated map. These functions were written specifically for evaluation with the Plank emission data.
  
- **superpixel_rotations.py**

  To be used with *generate_maps* and *testing_superpixels* notebooks

  This is the main file containing all the functions for working with rotating the emission and extinction data in preparation for the temperature fitting.

## Jupyter notebooks
- **compare_plank_rotations**

  This notebook uses functions from *comparison_functions.py* to generate residual, ratio and histogram plots of a rotated map. The notebook was written to be run from the dustfizz3d directory.
  
  *Note: requires a Dropbox access token to save generated figures. Info on how to generate one can be found here: https://help.displayr.com/hc/en-us/articles/360004116315-How-to-Create-an-Access-Token-for-Dropbox*
  
- **testing_superpixels**
  
  This notebook uses functions from *superpixel_rotations.py* to carry out the possible unique rotations on a test map, based on customizable nsides. Outputs hdf5 files of rotated maps in *rotated_maps* folder. Proceeds to read one of these rotated files to rotate back and check residuals.
  
- **generate_maps**

  This notebook uses functions from *superpixel_rotations.py* to rotate the plank emission maps. Outputs hdf5 files of rotated maps in *rotated_maps* folder. Also includes some tests for rotating back superpixels of temperature data
 
