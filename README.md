# cold_pool_variogram

## Overview
Repository for code and data used for paper "How variable are cold pools?" by Leah D. Grant, Bastian Kirsch, Jennie Bukowski, Nicholas M. Falk, Christine A. Neumaier, Mirjana Sakradzija, Susan C. van den Heever, and Felix Ament (2023), submitted to Geophysical Research Letters  

## File description
### Code
- **variogram_paper_obs_icon.py**: Main script for reading and analyzing FESSTVaL 2021 observations and ICON-LES data
- **variogram_paper_routines.py**: Routines for data reading and plotting used by variogram_paper_obs_icon.py
- **fesstval_routines.py**: Routines for reading, processing, and analyzing FESSTVaL observation data used by variogram_paper_obs_icon.py, variogram_paper_routines.py, and read_icon_data.ipynb
- **emprical_variogram.py**: Routines for calculation of empirical variograms
- **read_icon_data.ipynb**: Main script for processing ICON-LES data
- **icon_routines.py**: Routines for reading and processing ICON-LES data used by read_icon_data.ipynb
- **cp_spatial_analysis.py**: Routines needed for spatial interpolation of ICON data used by read_icon_data.ipynb
- **variogram_figures_2_3.m**: Matlab script for reading variogram netcdf files from observations and array of simulations, and creating Figures 2 and 3 of the manuscript
- **matlab_functions/**: sub-folder containing Matlab functions and colormap arrays used by variogram_figures_2_3.m script

### RAMS_Code
- **fig1_data_v04.ipynb**: Creates NetCDF file containing data plotted in figure 1.
- **get_network_Ts.py**: Reads in RAMS output, linerally interpolates temperature data to the FESSTVaL network, and saves this data as H5 files as an intermediate step. All of the other scripts in RAMS_code then use these H5 files.
- **Tref_dTmean_dTmin_stds.ipynb**: Calculates the reference temperature, mean and minimum temperatures, and mean and minimum standard devations. All information included in table 1. 
- **variogram_data_v02.ipynb**: Creates NetCDF file containing data plotted in figure 3.

### Data
- **ta_network_icon_dom0{x}_20210629.nc**: Simulated air temperature from ICON-LES case study simulations on 29 June 2021 interpolated to locations of FESSTVaL station network (created by read_icon_data.ipynb)
- **variogram_obs_20210629.nc**: Observed variograms of near-surface air temperature during FESSTVaL on 29 June 2021 at 1-min resolution (created by variogram_paper_obs_icon.py)
- **variogram_obs_20210517-20210827.nc**: Observed variograms for entire FESSTVaL period at 15-min resolution (created by variogram_paper_obs_icon.py)
- **variogram_icon_{xxx}m_20210629.nc**: Simulated variograms for ICON-LES case study simulations on 29 June 2021 at 15-min resolution (created by variogram_paper_obs_icon.py)
- **CS-TropOce-{xxxx}_variogram-data_v02.nc**: Simulated variograms for CS-TropOce simulations (used in Fig. 3)
- **variograms_IDEAL-DryBL-{xxx}.nc**: Simulated variograms for IDEAL-DryBL simulations (used in Fig. 3)
- **fval_uhh_apollowxt_l3_ta_v01_20210629.nc**: Spatially interpolated temperature observations during FESSTVaL on 29 June 2021 at 1-min resolution (used for Fig. 1, see https://doi.org/10.5281/ZENODO.8161036 for interpolation method)
- **ta_grid_icon_dom03_20210629.nc**: Simulated air temperature from ICON-LES case study simulations on 29 June 2021 interpolated to Cartesian grid (used for Fig. 1, created by read_icon_data.ipynb)
- **IDEAL-DryBL-50_fig-1-data.nc**: Simulated air temperature from IDEAL-DryBL-50m simulation (used for Fig. 1)

### Meta data
- **network_coordinates.txt**: x/y coordinates of FESSTVaL 2021 station network used by emprical_variogram.py
- **cold_pools_fesstval.txt**: List of cold pool events during FESSTVaL 2021 used by variogram_paper_obs_icon.py


## Contact
Bastian Kirsch (bastian.kirsch@uni-hamburg.de)<br>
Universität Hamburg, Hamburg, Germany

Leah D. Grant (leah.grant@colostate.edu)<br>
Colorado State University, Fort Collins, CO, United States

6 October 2023
