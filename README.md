# Oyster_salinity_var
Modeling oyster predator-prey dynamics under increasing salinity variability

Matlab code used in the Commander et al. manuscript describing oyster predator-prey dynamics under increasing salinity variability in Apalachicola Bay

Files:

# Setting up the different salinity scenarios:
TimeSeries_Resampling.m Analyzes salinity data from the NERR database to create a climatology and associated data needed for the resampling analysis
Salinity_runs.m is a wrapper function that loops over scenario options to generate the simulated timeseries
Create_Mock_Salinity.m uses the Denny et al. method to generate new simulated timeseries based on the Apalachicola climatology

# Running the oyster model
Oyster_run.m is a wrapper function that calls all of the needed functions for a given simulation scenario
Oyster_IPM_Disturbance.m performs one simulation of the model, given the inputs called by Oyster_run; it calls a number of functions listed below in #model files
Oyster_IPM_Disturbance_Analysis.m creates the figures in the manuscript

# model files
Oyster_Params.m creates a structure file with all parameters used in the model
mkkern.m creates IPM kernels of different types by calling kernmatSimp.m
kernmatSimp.m creates IPM kernels, given options
makeSimpVec.m creates a vector that will perform Simpson's integration 
predator_salinity_penalty.m calculates the appropriate predator attack rate, given a salinity and temperature

# Ancillary analyses
Apalach_data_comparison.m performs the comparison between historical data and model size distributions shown in Fig. S7
salinity_resid_analysis.m performs analyses on the time scales of variability in the historical salinity timeseries, including the FFT analysis
Oyster_global_sens.m performs the global sensitivity analysis
Oyster_Params_Sens.m creates the randomly generated parameter sets used in the sensitivity analysis
global_sens_randomforest.R performs the random forest analysis of for the global sensitivity results


